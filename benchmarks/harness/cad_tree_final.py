"""Build a frozen private CAD snapshot and time complete scalar/SIMD checks.

Run from the repository root with Pixi active:
    python -m benchmarks.harness.cad_tree_final --output .cache/cad-tree

The build manifest binds the executable to copied sources and APN binary.
--reuse-build verifies that manifest before measuring an existing executable.
Every successful repetition checks the entire decomposition twice; sign-set
comparison with --baseline-answers is an additional regression check.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess

from . import cad, measure
from .cad_revision15 import INITIAL, digest

ROOT = Path(__file__).resolve().parents[2]
DRIVER = Path('benchmarks/drivers/poly_mojo/cad_tree_final.mojo')


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def build(directory, reuse, *, driver=DRIVER, binary_name='cad-tree', extra_sources=()):
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / 'manifest.json'
    binary = directory / binary_name
    if reuse:
        record = json.loads(manifest.read_text())
        if digest(binary) != record['binary_sha256']:
            raise RuntimeError('executable does not match build manifest')
        for name, checksum in record['sources'].items():
            if digest(directory / name) != checksum:
                raise RuntimeError(f'frozen build input changed: {name}')
        return binary, record
    sources = sorted((ROOT / 'src').rglob('*.mojo'))
    sources += sorted((ROOT / 'tests/support').rglob('*.mojo'))
    sources += [ROOT / driver, ROOT / 'build/apn_mojo.mojoc', ROOT / 'build/apn_mojo.commit',
                ROOT / 'pixi.lock', ROOT / 'scripts/deps.py', Path(__file__).resolve(),
                ROOT / 'benchmarks/pixi.lock', Path(cad.__file__).resolve(),
                ROOT / 'benchmarks/harness/cad_revision15.py',
                ROOT / 'benchmarks/harness/polys.py',
                ROOT / 'benchmarks/harness/data/cad_example_bank_qepcad.txt',
                Path(measure.__file__).resolve(), *extra_sources]
    hashes = {}
    for source in sources:
        relative = source.relative_to(ROOT)
        target = directory / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        hashes[str(relative)] = digest(target)
    command = ['mojo', 'build', '--num-threads', '1', '-O3', '-I', str(directory / 'src'),
               '-I', str(directory / 'build'), '-I', str(directory / 'tests'),
               str(directory / driver), '-o', str(binary)]
    with (directory / 'compile.log').open('w') as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(command=command, sources=hashes, binary_sha256=digest(binary),
                  compiler=subprocess.check_output(['mojo', '--version'], text=True).strip(),
                  revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  dirty=bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True)))
    save(manifest, record)
    return binary, record


def signatures(path, legacy):
    lines = path.read_text().splitlines()
    return sorted({line.split(' ')[1] for line in lines[1:]}) if legacy else sorted(set(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--build', type=Path, help='Reuse a separately frozen build without overwriting its measurements')
    parser.add_argument('--cases', default=INITIAL)
    parser.add_argument('--modes', default='revision15,interval')
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--timeout', type=float, default=900)
    parser.add_argument('--cpu', type=int, default=min(os.sched_getaffinity(0)))
    parser.add_argument('--baseline-answers', type=Path)
    parser.add_argument('--reuse-build', action='store_true')
    parser.add_argument('--build-only', action='store_true')
    args = parser.parse_args()
    modes = args.modes.split(',')
    if args.repetitions < 1 or args.timeout <= 0 or args.cpu not in os.sched_getaffinity(0):
        parser.error('invalid repetitions, timeout, or CPU affinity')
    if not modes or any(mode not in ('lazard', 'revision15', 'interval') for mode in modes):
        parser.error('modes must be lazard, revision15, or interval')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    binary, manifest = build(args.build.resolve() if args.build else args.output / 'build', args.reuse_build or bool(args.build))
    if args.build_only:
        return
    cad.generate()
    cases = [case for case in cad.load() if case['kind'] != 'sat' and re.search(args.cases, case['id'])]
    if not cases:
        parser.error('no full-CAD cases selected')
    report = dict(schema=1, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope='private complete CAD; construction includes final-data assembly; terminal sections remain relative',
                  build=manifest, command=os.sys.argv, cpu=args.cpu, platform=platform.platform(),
                  repetitions=args.repetitions, warmup_repetitions=0,
                  cpu_model=next((line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                                  if line.startswith('model name')), 'unknown'),
                  note='One pinned execution CPU; shared-host contention is not controlled. No public consumer/materialization timing.',
                  results=[])
    failed = False
    for case in cases:
        work = args.output / case['id']
        work.mkdir(exist_ok=True)
        names = [f'x{i}' for i in range(len(case['variables']))]
        expressions = [cad.infix(poly, names, times='*') for poly in cad.polys_of(case)]
        expressions = [('-' + expr[2:]) if expr.startswith('- ') else expr for expr in expressions]
        inp = work / 'input.txt'
        inp.write_text(str(len(names)) + '\n' + '\n'.join(expressions) + '\n')
        first_answer = None
        for mode in modes:
            answer = work / f'{mode}.answer'
            stdout, stderr = work / f'{mode}.stdout', work / f'{mode}.stderr'
            resource = measure.run([str(binary), mode, str(inp), str(answer), str(args.repetitions)],
                                   cwd=ROOT, env=os.environ.copy(), timeout_s=args.timeout,
                                   memory_bytes=None, cpu=args.cpu, stdout=stdout, stderr=stderr).as_dict()
            item = dict(case=case['id'], mode=mode, variables=case['variables'], input_sha256=digest(inp), resource=resource)
            progress = Path(str(answer) + '.progress')
            if progress.exists():
                item['progress'] = []
                for line in progress.read_text().splitlines():
                    try:
                        item['progress'].append(json.loads(line))
                    except json.JSONDecodeError:
                        item['progress_truncated'] = True
            report['results'].append(item)
            if resource['returncode'] == 0 and not resource['timed_out']:
                rows = [json.loads(line) for line in stdout.read_text().splitlines()]
                item['records'] = rows
                if len(rows) != args.repetitions or not all(row['check_ok'] for row in rows):
                    raise RuntimeError('incomplete or unchecked measurement')
                item['medians'] = {key: statistics.median(row[key] for row in rows)
                                   for key in rows[0] if key.endswith('_ns')}
                for method in ('scalar', 'simd'):
                    item['medians'][f'certified_{method}_ns'] = statistics.median(
                        row['construct_ns'] + row[f'verify_{method}_ns'] for row in rows)
                item['answer_sha256'] = digest(answer)
                current = signatures(answer, False)
                if first_answer is not None:
                    item['mode_signatures_agree'] = current == first_answer
                    failed |= not item['mode_signatures_agree']
                first_answer = current
                if args.baseline_answers:
                    baseline = args.baseline_answers / f"{case['id']}-lazard.answer"
                    item['baseline_signatures_agree'] = current == signatures(baseline, True)
                    failed |= not item['baseline_signatures_agree']
            else:
                failed = True
                item['stderr'] = stderr.read_text(errors='replace')[-4000:]
            save(args.output / 'report.json', report)
            print(case['id'], mode, item.get('medians', resource), flush=True)
    if failed:
        raise SystemExit('one or more runs failed, timed out, or disagreed; see report.json')


if __name__ == '__main__':
    main()
