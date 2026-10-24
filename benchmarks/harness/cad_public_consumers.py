"""Measure facade output costs and unchanged selective solving separately.

Compile the existing CAD driver at O3 with NVARS=3. Pass that binary here; this
runner does not compile anything while performance measurements are active.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import re
import statistics

from . import cad, measure
from .cad_revision15 import digest
from .cad_tree_final import save
from .polys import write_records

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--cases', default='^(full-bank-(off-center-ellipse|mccallum-trivariate-random-polynomial)|sat-bank-(ball-and-circular-cylinder|collins-and-johnson))$')
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--timeout', type=float, default=1800)
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    cad.generate()
    selected = [case for case in cad.load() if re.search(args.cases, case['id'])]
    if not selected or any(len(case['variables']) != 3 for case in selected):
        parser.error('select three-variable cases for this driver')
    report = dict(schema=1, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  binary_sha256=digest(args.binary), cpu=args.cpu, command=os.sys.argv,
                  scope='facade decompose plus explicit output costs; selective solve retains its legacy backend',
                  notes=['Decomposition times_ns already includes mandatory verification; check_ns is an extra recheck.',
                         'cells_ns includes absolute materialization and facade description together. absolute_materialize_ns measures core geometry independently; do not sum these overlapping output requests.',
                         'archive_write_ns and archive_read_ns include their required acceptance checks. Formatting and in-memory iteration are timed separately.',
                         'Each case runs in a separate process, with four call repetitions; output operations run once.'], results=[])
    failed = False
    for case in selected:
        work = args.output / case['id']
        work.mkdir(exist_ok=True)
        inp = work / 'input.txt'
        inp.write_text(write_records(cad.polys_of(case), 3, None))
        mode = 'solve' if case['kind'] == 'sat' else 'consumers'
        masks = ','.join(map(str, case['masks'])) if mode == 'solve' else '-'
        result, answer = work / 'result.json', work / 'answer.txt'
        # Negative MIN_SECONDS selects exactly MIN_REPS calls.
        command = [str(args.binary), mode, str(inp), str(result), str(answer), masks, '-1', '4']
        resource = measure.run(command, cwd=ROOT, env=os.environ.copy(), timeout_s=args.timeout,
                               memory_bytes=None, cpu=args.cpu, stdout=work / 'stdout', stderr=work / 'stderr').as_dict()
        item = dict(case=case['id'], mode=mode, input_sha256=digest(inp), variables=case['variables'], resource=resource)
        report['results'].append(item)
        if resource['returncode'] == 0 and not resource['timed_out']:
            raw = json.loads(result.read_text())
            item['measurement'] = raw
            item['median_call_ns'] = statistics.median(raw['times_ns'])
            item['answer_sha256'] = digest(answer)
            failed |= not raw['check_ok'] or len(raw['times_ns']) != 4
            if mode == 'solve' and raw['sat']:
                item['independent_sample_ok'] = cad.check_point(cad.polys_of(case), case['masks'], answer.read_text().splitlines()[1].split(' '))
                failed |= not item['independent_sample_ok']
        else:
            failed = True
            item['stderr'] = (work / 'stderr').read_text(errors='replace')[-4000:]
        save(args.output / 'report.json', report)
        print(case['id'], mode, item.get('median_call_ns', resource), flush=True)
    report['passed'] = not failed
    save(args.output / 'report.json', report)
    if failed:
        raise SystemExit('consumer measurement failed; see report.json')


if __name__ == '__main__':
    main()
