"""Compile invalid programs separately and require their expected diagnostics."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from build_state import prepared_build
from packages import ROOT
from suite import _first_line, _report
from test_executable import run


def check(path, snapshot, timeout):
    expected = _first_line(path, '# expect:')
    if expected is None:
        return dict(name=str(path.relative_to(ROOT)), passed=False, detail='first line must be "# expect: ..."', seconds=0)
    output = ROOT / '.cache/compile_fail'
    output.mkdir(parents=True, exist_ok=True)
    built = run(['mojo', 'build', '-I', str(snapshot.path), '-I', str(ROOT / 'tests'),
                 str(path), '-o', str(output / path.stem)], ROOT, timeout)
    diagnostic = built['stdout'] + built['stderr']
    passed = built['returncode'] not in (0, None) and expected in diagnostic
    detail = '' if passed else ('compiled, but must not' if built['returncode'] == 0 else f'error does not contain {expected!r}:\n{diagnostic}')
    return dict(name=str(path.relative_to(ROOT)), passed=passed, detail=detail, seconds=built['seconds'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', default='build')
    parser.add_argument('--report', type=Path, default=ROOT / '.cache/tests/rejections.json')
    parser.add_argument('--jobs', type=int, default=2, help='Concurrent rejection compiles')
    parser.add_argument('--timeout', type=float, default=600)
    parser.add_argument('pattern', nargs='?')
    args = parser.parse_args()
    paths = [p for p in sorted((ROOT / 'tests/compile_fail').glob('*.mojo')) if not args.pattern or args.pattern in p.name]
    if not paths:
        parser.error('No compiler-rejection tests match')
    records = []
    with prepared_build(args.build) as snapshot:
        print(f'Using library build {snapshot.identity}: {snapshot.path}', flush=True)
        _report(args.report, snapshot, records, complete=False, total=len(paths))
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            pending = {pool.submit(check, path, snapshot, args.timeout): path for path in paths}
            for future in as_completed(pending):
                try:
                    record = future.result()
                except Exception as error:
                    record = dict(name=str(pending[future].relative_to(ROOT)), passed=False, detail=str(error), seconds=None)
                records.append(record)
                print(f"{'PASS' if record['passed'] else 'FAIL'} {record['name']} [{len(records)}/{len(paths)}]", flush=True)
                if not record['passed']:
                    print(record['detail'][-4000:], flush=True)
                _report(args.report, snapshot, records, complete=False, total=len(paths))
    _report(args.report, snapshot, records, complete=True, total=len(paths))
    return not all(r['passed'] for r in records)


if __name__ == '__main__':
    raise SystemExit(main())
