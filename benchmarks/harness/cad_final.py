"""Compile and measure the private rational final-stack prototype on one CPU.

From the repository root, with its Pixi environment active:
    pixi run --locked python -m benchmarks.harness.cad_final

The JSON report includes raw repetitions, medians, source/input hashes and
process memory. This measures a QQ stack, not multivariate CAD or QEPCAD B.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess

from . import measure

ROOT = Path(__file__).resolve().parents[2]
DRIVER = Path('benchmarks/drivers/poly_mojo/cad_final.mojo')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def output(command):
    return subprocess.check_output(command, cwd=ROOT, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repetitions', type=int, default=31)
    parser.add_argument('--cpu', type=int, default=min(os.sched_getaffinity(0)))
    parser.add_argument('--output', type=Path, default=ROOT / '.cache/revision-16/measurement.json')
    args = parser.parse_args()
    if args.repetitions < 3:
        parser.error('at least three repetitions are required')
    if args.cpu not in os.sched_getaffinity(0):
        parser.error('CPU is outside this process affinity')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cache = ROOT / '.cache/revision-16'
    cache.mkdir(parents=True, exist_ok=True)
    binary = cache / 'cad_final_bench'
    sources = sorted((ROOT / 'src').rglob('*.mojo')) + [ROOT / DRIVER, Path(__file__).resolve(),
              Path(measure.__file__).resolve(), ROOT / 'pixi.lock', ROOT / 'scripts/deps.py']
    source_hashes = {str(p.relative_to(ROOT)): digest(p) for p in sources}
    command = ['mojo', 'build', '--num-threads', '2', '-O3', '-I', 'src', '-I', 'build', str(DRIVER), '-o', str(binary)]
    subprocess.run(command, cwd=ROOT, check=True)
    stdout, stderr = cache / 'benchmark.jsonl', cache / 'benchmark.stderr'
    process = measure.run([str(binary), str(args.repetitions)], cwd=ROOT, env=os.environ.copy(),
                          timeout_s=120, memory_bytes=0, cpu=args.cpu, stdout=stdout, stderr=stderr)
    if process.returncode != 0 or process.timed_out:
        raise RuntimeError(f'benchmark failed: {process.as_dict()}\n{stderr.read_text()}')
    if source_hashes != {str(p.relative_to(ROOT)): digest(p) for p in sources}:
        raise RuntimeError('sources changed during measurement')
    records = [json.loads(line) for line in stdout.read_text().splitlines()]
    summaries = {}
    for name in sorted({row['case'] for row in records}):
        rows = [row for row in records if row['case'] == name]
        if len(rows) != args.repetitions or not all(row['check_ok'] for row in rows):
            raise RuntimeError(f'invalid repetitions for {name}')
        # Discard the first repetition of each case as warm-up; retain it above.
        measured = rows[1:]
        times = {key: statistics.median(row[key] for row in measured)
                 for key in measured[0] if key.endswith('_ns')}
        if 'input' in rows[0]:
            times['certified_scalar_ns'] = statistics.median(row['construct_and_assemble_ns'] + row['verify_scalar_ns'] for row in measured)
            times['certified_simd_ns'] = statistics.median(row['construct_and_assemble_ns'] + row['verify_simd_ns'] for row in measured)
            times['input_sha256'] = hashlib.sha256(rows[0]['input'].encode()).hexdigest()
        summaries[name] = times
    cpu_info = Path('/proc/cpuinfo').read_text()
    cpu_model = next((line.split(':', 1)[1].strip() for line in cpu_info.splitlines()
                      if line.startswith('model name')), platform.machine())
    report = dict(schema=1, date=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope='private complete one-variable QQ stacks; no projection; variable 0',
                  compiler=output(['mojo', '--version']), compile_command=command,
                  revision=output(['git', 'rev-parse', 'HEAD']),
                  dirty=bool(output(['git', 'status', '--porcelain'])),
                  source_sha256=source_hashes, binary_sha256=digest(binary),
                  dependency_sha256=digest(ROOT / 'build/apn_mojo.mojoc'),
                  cpu=args.cpu, cpu_model=cpu_model, platform=platform.platform(),
                  threads=1, repetitions=args.repetitions, warmup_repetitions=1,
                  process=process.as_dict(),
                  memory_note='maxrss_kb is for the whole benchmark process and can include the Python parent high-water mark; no per-case peak',
                  assembly_note='assembly is fused with construction; final_data_text_bytes is a deterministic text size, not a public archive',
                  medians=summaries, records=records)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(summaries, indent=2))
    print(f'Report: {args.output}')


if __name__ == '__main__':
    main()
