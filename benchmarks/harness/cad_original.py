"""Measure the frozen pre-revision full CAD and its complete checker.

The binary is built with cad_original.mojo against an untouched source archive.
Its manifest must list every frozen input and the binary hash. This runner does
not compile or use the current library. Each repetition is a separate process;
construction and checking timers exclude startup and answer aggregation.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import re

from . import cad, measure
from .cad_revision15 import INITIAL, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--answers', required=True, type=Path,
                        help='New public results whose final.answer sign sets must agree')
    parser.add_argument('--cases', default=INITIAL)
    parser.add_argument('--repetitions', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=1800)
    parser.add_argument('--cpu', type=int, default=0)
    args = parser.parse_args()
    if args.repetitions < 1 or args.timeout <= 0 or args.cpu not in os.sched_getaffinity(0):
        parser.error('invalid repetitions, timeout, or CPU affinity')
    binary = args.binary.resolve()
    manifest = json.loads(args.manifest.read_text())
    if digest(binary) != manifest['binary_sha256']:
        raise RuntimeError('frozen binary changed')
    for name, expected in manifest['sources'].items():
        if digest(binary.parent / name) != expected:
            raise RuntimeError('frozen build input changed: ' + name)
    cad.generate()
    cases = [case for case in cad.load()
             if case['kind'] != 'sat' and re.search(args.cases, case['id'])]
    if not cases:
        parser.error('no full-CAD cases selected')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = dict(schema=1, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope='original plain Lazard construction plus complete legacy verification',
                  build=manifest, harness_sha256=digest(Path(__file__)), command=os.sys.argv,
                  cpu=args.cpu, repetitions=args.repetitions, warmup_repetitions=0,
                  note='One pinned CPU with no overlapping task compiler or benchmark. '
                       'Shared-host contention remains possible. The legacy producer builds '
                       'absolute sample coordinates and intermediate proof histories.', results=[])
    failed = False
    for case in cases:
        work = output / case['id']
        work.mkdir(exist_ok=True)
        names = [f'x{i}' for i in range(len(case['variables']))]
        expressions = [cad.infix(p, names, times='*') for p in cad.polys_of(case)]
        expressions = [('-' + p[2:]) if p.startswith('- ') else p for p in expressions]
        inp = work / 'input.txt'
        inp.write_text(str(len(names)) + '\n' + '\n'.join(expressions) + '\n')
        expected_path = args.answers / case['id'] / 'final.answer'
        expected = set(expected_path.read_text().splitlines())
        for repetition in range(args.repetitions):
            answer = work / f'{repetition}.answer'
            stdout, stderr = work / f'{repetition}.stdout', work / f'{repetition}.stderr'
            command = [str(binary), str(inp), str(answer)]
            resource = measure.run(command, cwd=work, env=os.environ.copy(),
                                   timeout_s=args.timeout, memory_bytes=None, cpu=args.cpu,
                                   stdout=stdout, stderr=stderr).as_dict()
            item = dict(case=case['id'], repetition=repetition, variables=case['variables'],
                        input_sha256=digest(inp), resource=resource, command=command)
            if resource['returncode'] == 0 and not resource['timed_out']:
                row = json.loads(stdout.read_text())
                item.update(record=row, answer_sha256=digest(answer),
                            new_public_signatures_agree=set(answer.read_text().splitlines()) == expected)
                failed |= not row['check_ok'] or not item['new_public_signatures_agree']
            else:
                failed = True
                item.update(stderr=stderr.read_text(errors='replace')[-4000:],
                            progress=Path(str(answer) + '.progress').read_text())
            report['results'].append(item)
            (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
            print(case['id'], repetition, item.get('record', resource), flush=True)
    report['utc_finished'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    report['passed'] = not failed
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    if failed:
        raise SystemExit('original CAD failed or disagreed; see report.json')


if __name__ == '__main__':
    main()
