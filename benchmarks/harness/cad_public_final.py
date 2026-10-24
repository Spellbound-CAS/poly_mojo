"""Measure accepted public CAD calls from a frozen O3 source snapshot.

Use cad_tree_final on the same sources for separate producer/checker phases,
and cad_revision15 --qepcad-bin for the matching reference full-CAD workload.
This driver times construction plus mandatory acceptance exactly once. The
extra scalar recheck and sign-set aggregation are outside that timed region.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import platform
import re
import statistics

from . import cad, measure
from .cad_revision15 import INITIAL, digest
from .cad_tree_final import build, save

ROOT = Path(__file__).resolve().parents[2]
DRIVER = Path('benchmarks/drivers/poly_mojo/cad_public_final.mojo')


def historical_signatures(directory, case):
    """Accept either legacy absolute-output layout or private sign-only output."""
    candidates = [(directory / f'{case}-lazard.answer', True),
                  (directory / case / 'lazard.answer', True),
                  (directory / case / 'interval.answer', False)]
    for path, legacy in candidates:
        if path.exists():
            lines = path.read_text().splitlines()
            return path, sorted({line.split(' ')[1] for line in lines[1:]} if legacy else set(lines))
    return None, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--build', type=Path, help='Reuse a separately frozen build directory')
    parser.add_argument('--build-only', action='store_true')
    parser.add_argument('--cases', default=INITIAL)
    parser.add_argument('--modes', default='final')
    parser.add_argument('--repetitions', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=1800)
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--baseline-answers', action='append', type=Path, default=[])
    args = parser.parse_args()
    if args.repetitions < 1 or args.timeout <= 0 or args.cpu not in os.sched_getaffinity(0):
        parser.error('invalid repetitions, timeout, or CPU affinity')
    modes = args.modes.split(',')
    if not modes or any(mode not in ('final', 'lazard', 'mccallum') for mode in modes):
        parser.error('unknown mode')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    build_directory = args.build.resolve() if args.build else args.output / 'build'
    binary, manifest = build(build_directory, bool(args.build), driver=DRIVER,
                             binary_name='cad-public', extra_sources=[Path(__file__).resolve()])
    if args.build_only:
        return
    cad.generate()
    cases = [case for case in cad.load() if case['kind'] != 'sat' and re.search(args.cases, case['id'])]
    if not cases:
        parser.error('no full-CAD cases selected')
    report = dict(schema=1, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope='public cad_final: construction + assembly + mandatory complete SIMD verification; relative terminal sections',
                  build=manifest, command=os.sys.argv, cpu=args.cpu, platform=platform.platform(),
                  repetitions=args.repetitions, warmup_repetitions=0,
                  note='One pinned CPU; no other task benchmarks or compilation overlap. Shared host contention remains possible. '
                       'The first repetition also runs an extra scalar check, reported separately and excluded from accepted_ns. '
                       'No absolute output materialization or archive writing is included.', results=[])
    failed = False
    for case in cases:
        work = args.output / case['id']
        work.mkdir(exist_ok=True)
        names = [f'x{i}' for i in range(len(case['variables']))]
        expressions = [cad.infix(poly, names, times='*') for poly in cad.polys_of(case)]
        expressions = [('-' + value[2:]) if value.startswith('- ') else value for value in expressions]
        inp = work / 'input.txt'
        inp.write_text(str(len(names)) + '\n' + '\n'.join(expressions) + '\n')
        for mode in modes:
            answer = work / f'{mode}.answer'
            stdout, stderr = work / f'{mode}.stdout', work / f'{mode}.stderr'
            resource = measure.run([str(binary), mode, str(inp), str(answer), str(args.repetitions)],
                                   cwd=ROOT, env=os.environ.copy(), timeout_s=args.timeout,
                                   memory_bytes=None, cpu=args.cpu, stdout=stdout, stderr=stderr).as_dict()
            item = dict(case=case['id'], mode=mode, variables=case['variables'], input_sha256=digest(inp), resource=resource)
            report['results'].append(item)
            if resource['returncode'] == 0 and not resource['timed_out']:
                rows = [json.loads(line) for line in stdout.read_text().splitlines()]
                if len(rows) != args.repetitions or not all(row['check_ok'] for row in rows):
                    raise RuntimeError('incomplete or unchecked measurement')
                item['records'] = rows
                item['medians'] = dict(accepted_ns=statistics.median(row['accepted_ns'] for row in rows))
                item['range_ns'] = [min(row['accepted_ns'] for row in rows), max(row['accepted_ns'] for row in rows)]
                item['answer_sha256'] = digest(answer)
                signatures = sorted(set(answer.read_text().splitlines()))
                for directory in args.baseline_answers:
                    path, expected = historical_signatures(directory, case['id'])
                    if path:
                        item['baseline_answer'] = str(path)
                        item['baseline_answer_sha256'] = digest(path)
                        item['baseline_signatures_agree'] = signatures == expected
                        failed |= not item['baseline_signatures_agree']
                        break
            else:
                failed = True
                item['stderr'] = stderr.read_text(errors='replace')[-4000:]
                progress = Path(str(answer) + '.progress')
                if progress.exists():
                    item['progress'] = progress.read_text()
            save(args.output / 'report.json', report)
            print(case['id'], mode, item.get('medians', resource), flush=True)
    report['utc_finished'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    report['passed'] = not failed
    save(args.output / 'report.json', report)
    if failed:
        raise SystemExit('a call failed, timed out, or disagreed; see report.json')


if __name__ == '__main__':
    main()
