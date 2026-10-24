"""Reproduce the private full-CAD projection comparison (not SAT timings).

Build the driver against the sources under test first; --binary identifies that
executable. Every Mojo repetition checks its complete certificate. Independent
sample checks and realized sign sets supplement, but do not replace, that check.
QEPCAD's d-stat construction CPU total includes normalization, projection,
lifting and garbage collection. Its post-initialization CPU clock and matching
trivial-input overhead are also retained separately from process startup.
"""
import argparse
import datetime
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import re
import shlex
import subprocess

from . import cad, measure

ROOT = Path(__file__).resolve().parents[2]
INITIAL = ('full-bank-(mccallum-trivariate-random-polynomial|off-center-ellipse|'
           'simplified-putnum-example|x-axis-ellipse-problem)$')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command, work, name, args, env=None):
    return measure.run(command, cwd=work, env=env or os.environ.copy(),
                       timeout_s=args.timeout, memory_bytes=None, cpu=args.cpu,
                       stdout=work / f'{name}.stdout', stderr=work / f'{name}.stderr').as_dict()


def _sample_worker(polys, answer, limit, queue):
    try:
        checked, total, valid = cad.check_cells(polys, answer, limit)
        queue.put(dict(checked=checked, total=total, valid=valid))
    except Exception as error:
        queue.put(dict(error=str(error), valid=None))


def sample_check(polys, answer, args):
    queue = multiprocessing.Queue()
    worker = multiprocessing.Process(target=_sample_worker, args=(polys, answer, args.samples, queue))
    worker.start()
    worker.join(args.sample_timeout)
    if worker.is_alive():
        worker.terminate()
        worker.join()
        return dict(valid=None, timed_out=True, seconds=args.sample_timeout)
    result = queue.get(timeout=2) if worker.exitcode == 0 else dict(valid=None, exitcode=worker.exitcode)
    queue.close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--qepcad-bin', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--cases', default=INITIAL)
    parser.add_argument('--modes', default='revision15,revision15-finite')
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--samples', type=int, default=64)
    parser.add_argument('--sample-timeout', type=float, default=60)
    parser.add_argument('--resume', action='store_true', help='reuse complete raw repetitions with the same binary and input')
    parser.add_argument('--timeout', type=float, default=1200)
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--baseline-answers', type=Path)
    args = parser.parse_args()
    if args.repetitions < 1 or args.samples < 1:
        parser.error('repetitions and samples must be positive')
    if not args.binary and not args.qepcad_bin:
        parser.error('provide a Mojo binary or QEPCAD binary')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    previous = json.loads((args.output / 'report.json').read_text()) if args.resume and (args.output / 'report.json').exists() else {}
    cad.generate()
    cases = [c for c in cad.load() if c['kind'] != 'sat' and re.search(args.cases, c['id'])]
    if not cases:
        parser.error('no cases selected')
    cpu = next((line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                if line.startswith('model name')), 'unknown')
    report = dict(schema=1, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  harness_sha256=digest(Path(__file__)),
                  input_generator_sha256=digest(Path(cad.__file__)),
                  corpus_sha256=digest(cad.BANK),
                  benchmark_lock_sha256=digest(ROOT / 'benchmarks/pixi.lock'),
                  revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  cpu=cpu, pinned_cpu=args.cpu, threads=1, repetitions=args.repetitions,
                  command=os.sys.argv, results=[],
                  sources={str(p.relative_to(ROOT)): digest(p) for p in sorted((ROOT / 'src').rglob('*.mojo'))})
    if args.binary:
        args.binary = args.binary.resolve()
        report['binary_sha256'] = digest(args.binary)
    if args.qepcad_bin:
        args.qepcad_bin = args.qepcad_bin.resolve()
        qe = args.qepcad_bin.parent.parent
        report['qepcad'] = dict(binary_sha256=digest(args.qepcad_bin), heap_cells=100000000,
                               revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=qe, text=True).strip())
    for case in cases:
        work = args.output / case['id']
        work.mkdir(exist_ok=True)
        names = [f'x{i}' for i in range(len(case['variables']))]
        expressions = [cad.infix(p, names, times='*') for p in cad.polys_of(case)]
        expressions = [('-' + p[2:]) if p.startswith('- ') else p for p in expressions]
        inp = work / 'input.txt'
        inp.write_text(str(len(names)) + '\n' + '\n'.join(expressions) + '\n')
        if args.binary:
            for mode in args.modes.split(','):
                answer = work / f'{mode}.answer'
                saved = next((r for r in previous.get('results', []) if r['case'] == case['id'] and r['mode'] == mode), {})
                rows = []
                if saved.get('input_sha256') == digest(inp) and previous.get('binary_sha256') == report['binary_sha256'] and answer.exists():
                    try:
                        rows = [json.loads(s) for s in (work / f'{mode}.stdout').read_text().splitlines()]
                    except (OSError, ValueError):
                        rows = []
                if len(rows) == args.repetitions and all(r.get('check_ok') for r in rows):
                    resource = saved.get('resource', {})
                    if resource.get('returncode') != 0:
                        resource = dict(returncode=0, resumed=True, resource_measurement_unavailable=True)
                else:
                    resource = run([str(args.binary), mode, str(inp), str(answer), str(args.repetitions)], work, mode, args)
                item = dict(case=case['id'], mode=mode, variables=case['variables'], input_sha256=digest(inp), resource=resource)
                report['results'].append(item)
                if resource['returncode'] == 0:
                    rows = [json.loads(s) for s in (work / f'{mode}.stdout').read_text().splitlines()]
                    item.update(repetitions=rows, answer_sha256=digest(answer))
                    if len(rows) != args.repetitions or not all(r['check_ok'] for r in rows):
                        raise RuntimeError('missing or rejected certificate repetition')
                    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
                    item['independent_samples'] = sample_check(cad.polys_of(case), answer, args)
                    if item['independent_samples'].get('valid') is False:
                        raise RuntimeError('independent sample sign mismatch')
                    if args.baseline_answers:
                        baseline = args.baseline_answers / f"{case['id']}-lazard.answer"
                        if baseline.exists():
                            signatures = lambda p: sorted({line.split(' ')[1] for line in p.read_text().splitlines()[1:]})
                            item['baseline_signatures_agree'] = signatures(answer) == signatures(baseline)
                            if not item['baseline_signatures_agree']:
                                raise RuntimeError('realized input sign set changed')
                (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
                print(case['id'], mode, resource['returncode'], flush=True)
        if args.qepcad_bin:
            for trivial in (True, False):
                label = 'qepcad-trivial' if trivial else 'qepcad'
                script = work / f'{label}.input'
                script.write_text(cad.Qepcad('fullcad').script(case, trivial=trivial))
                for repetition in range(args.repetitions):
                    name = f'{label}-{repetition}'
                    command = ['sh', '-c', f'exec {shlex.quote(str(args.qepcad_bin))} +N100000000 < {shlex.quote(str(script))}']
                    resource = run(command, work, name, args, {**os.environ, 'qe': str(qe)})
                    text = (work / f'{name}.stdout').read_text(errors='replace')
                    elapsed = re.search(r'System time after the initialization: (\d+) milliseconds', text)
                    item = dict(case=case['id'], mode=label, repetition=repetition, resource=resource,
                                input_sha256=digest(script), variables=case['variables'], command=command)
                    if resource['returncode'] == 0 and elapsed:
                        item['cpu_after_initialization_ns'] = int(elapsed[1]) * 1000000
                        # d-stat is emitted before solution output. Its total
                        # includes normalization, projection, lifting and GC,
                        # without timing initialization or printed factors.
                        for key, pattern in [
                            ('construction_cpu_ns', r'^Total time with garbage collection\s*:\s*(\d+)'),
                            ('garbage_collection_cpu_ns', r'^Total time for garbage collection\s*:\s*(\d+)'),
                            ('normalization_cpu_ns', r'^Time taken for Normalization Phase\s*:\s*(\d+)')]:
                            matched = re.search(pattern, text, re.M)
                            if matched:
                                item[key] = int(matched[1]) * 1000000
                        for key, phase_label in [('projection_cpu_ns', 'PROJECT'), ('lifting_cpu_ns', 'TICAD')]:
                            matched = re.search(r'^' + phase_label + r'((?:\s+\d+)+)\s*$', text, re.M)
                            if matched:
                                item[key] = int(matched[1].split()[-1]) * 1000000
                        for key, pattern in [('cells_by_level', 'Cells'), ('factors_by_level', 'Proj Fact')]:
                            matched = re.search(r'^' + pattern + r'((?:\s+\d+)+)\s*$', text, re.M)
                            if matched:
                                item[key] = [int(v) for v in matched[1].split()[:-1]]
                    else:
                        item['failure_stdout'] = text[-5000:]
                        item['failure_stderr'] = (work / f'{name}.stderr').read_text(errors='replace')[-2000:]
                    report['results'].append(item)
                    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
                print(case['id'], label, flush=True)


if __name__ == '__main__':
    main()
