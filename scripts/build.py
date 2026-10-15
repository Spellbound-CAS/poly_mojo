"""Prepare one library build, including the diagnostics reused by test gates."""
import argparse
import json
import subprocess
import uuid

from packages import ROOT, present
from build_state import MANIFEST, SCHEMA, artifacts, build_lock, compiler_identity, inputs
from compiler_diagnostics import classify


def build(output, extra=()):
    with build_lock(output, exclusive=True) as output:
        manifest = output / MANIFEST
        manifest.unlink(missing_ok=True)
        if not (output / 'apn_mojo.mojoc').is_file():
            raise SystemExit('Missing apn_mojo. Prepare dependencies with `pixi run --locked deps`.')
        before = inputs()
        compiler = compiler_identity()
        docstrings = []
        for package in present(ROOT):
            # Collect doc diagnostics in the only precompile pass. Werror cannot
            # exclude compiler-generated docstrings, so reject other warnings below.
            cmd = ['mojo', 'precompile', '--diagnose-missing-doc-strings', '-I', str(output), *extra,
                   f'src/{package}', '-o', str(output / f'{package}.mojoc')]
            print(' '.join(cmd), flush=True)
            result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
            docs, warnings = classify(result.stdout + result.stderr, ROOT)
            docstrings.extend(docs)
            if result.returncode or warnings:
                raise SystemExit(f'Compilation failed for {package}:\n' + result.stdout + result.stderr)
        if inputs() != before or compiler_identity() != compiler:
            raise SystemExit('Build inputs changed during compilation. Run the build again.')
        record = dict(schema=SCHEMA, build_id=uuid.uuid4().hex, inputs=before,
                      compiler=compiler, artifacts=artifacts(output), docstrings=docstrings)
        manifest.write_text(json.dumps(record, indent=2) + '\n')
        print(f"Library build {record['build_id']}: {output}", flush=True)
        print(f'{len(docstrings)} missing or partial docstrings recorded', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build', help='Directory for the shared .mojoc files')
    build(parser.parse_args().output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
