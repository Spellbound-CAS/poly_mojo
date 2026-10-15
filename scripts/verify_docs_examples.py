"""Compile and execute documentation programs against the shared library build.

One executable per optimization setting covers the selected programs, with --Werror.
Its standard output must match the documented result. Reports identify the
prepared library; this runner never compiles packages or imports src/.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

from build_state import prepared_build
from test_executable import prepare, CompilationFailed
from docs_contract import ROOT, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, help='A fresh report directory')
    parser.add_argument('--build', default='build', help='Prepared library directory')
    parser.add_argument('--only', help='Run only programs whose path contains this text')
    parser.add_argument('--compiler-threads', type=int, default=2, help='Threads in the single compiler process (default: 2)')
    parser.add_argument('--compile-timeout', type=float, help='Optional timeout per suite compilation')
    parser.add_argument('--jobs', type=int, default=2, help='Concurrent example processes')
    args = parser.parse_args()
    with prepared_build(args.build) as snapshot:
        records, executables = check(args, snapshot, parser)
    passed = all(r['passed'] for r in records)
    (Path(args.output).resolve() / 'report.json').write_text(json.dumps(
        dict(passed=passed, library_build=snapshot.report(), executables=executables, records=records), indent=2))
    print('DOCUMENTATION EXAMPLES PASS' if passed else 'FAILED', flush=True)
    return not passed


def check(args, snapshot, parser):
    declared = data('examples.json')['examples']
    selected = {p: s for p, s in declared.items() if not args.only or args.only in p}
    if not selected:
        parser.error('No documentation programs match --only')
    folder = Path(args.output).resolve()
    folder.mkdir(parents=True, exist_ok=False)
    print(f'Using library build {snapshot.identity}: {snapshot.path}', flush=True)
    records, executables = [], []
    for level in (0, 3):
        try:
            executable = prepare([ROOT / p for p in selected], snapshot, ROOT,
                                 flags=('--Werror', f'-O{level}'), timeout=args.compile_timeout, compiler_threads=args.compiler_threads)
        except CompilationFailed as error:
            records.append(dict(source='suite', optimization=level, passed=False, compile=error.record['compile'], run=None))
            break
        executables.append(executable.report())

        def execute(path):
            executed = executable.execute(ROOT / path, ROOT, 120)
            passed = executed['returncode'] == 0 and executed['stdout'] == selected[path]['stdout']
            print('PASS' if passed else 'FAIL', f'{path}-package-O{level}', flush=True)
            return dict(source=path, route='package', optimization=level, passed=passed,
                        compile=executable.record['compile'], run=executed)

        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
            records.extend(executor.map(execute, selected))
    return records, executables


if __name__ == '__main__':
    raise SystemExit(main())
