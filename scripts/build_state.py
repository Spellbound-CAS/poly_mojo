"""Identify and lock the single library build consumed by all test runners."""
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess

from packages import ROOT, present, test_dependencies

MANIFEST = 'library-build.json'
SCHEMA = 1


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def compiler_identity():
    return subprocess.run(['mojo', '--version'], cwd=ROOT, check=True,
                          capture_output=True, text=True).stdout.strip()


def inputs():
    """Inputs that require a new library build; test and prose edits do not."""
    paths = sorted((ROOT / 'src').rglob('*.mojo'))
    paths += [ROOT / name for name in (
        'pixi.lock', 'scripts/build.py', 'scripts/build_state.py',
        'scripts/compiler_diagnostics.py', 'scripts/packages.py', 'scripts/deps.py')]
    return {str(path.relative_to(ROOT)): digest(path) for path in paths if path.is_file()}


def artifacts(output):
    names = ['apn_mojo.mojoc', 'apn_mojo.commit']
    names += [f'{name}.mojoc' for name in present(ROOT)]
    return {name: digest(output / name) for name in names}


def resolve(output):
    return (ROOT / output).resolve()


@contextmanager
def build_lock(output, *, exclusive=False):
    """Readers share packages; a build cannot replace them during a test run."""
    output = resolve(output)
    if exclusive:
        output.mkdir(parents=True, exist_ok=True)
    if not output.is_dir():
        raise SystemExit('No library build. Run `pixi run --locked build` first.')
    with (output / '.library.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
        try:
            yield output
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


@dataclass(frozen=True)
class LibraryBuild:
    path: Path
    manifest: dict

    @property
    def identity(self):
        return self.manifest['build_id']

    def report(self):
        return {'directory': str(self.path), 'build_id': self.identity,
                'compiler': self.manifest['compiler'], 'artifacts': self.manifest['artifacts']}


def validate(output):
    """Reject missing, stale, or replaced packages without compiling anything."""
    try:
        record = json.loads((output / MANIFEST).read_text())
        if record['schema'] != SCHEMA or not isinstance(record['build_id'], str):
            raise ValueError('unsupported build manifest')
        if record['inputs'] != inputs():
            raise ValueError('library sources or build inputs changed')
        if record['compiler'] != compiler_identity():
            raise ValueError('the Mojo toolchain changed')
        if record['artifacts'] != artifacts(output):
            raise ValueError('a library package or dependency changed')
        if not isinstance(record['docstrings'], list):
            raise ValueError('invalid docstring diagnostics')
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise SystemExit(f'No current library build in {output}: {error}. '
                         'Run `pixi run --locked build` once before testing.') from None
    return LibraryBuild(output, record)


@contextmanager
def prepared_build(output='build', *, kernel=False):
    with build_lock(output) as folder:
        snapshot = validate(folder)
        if kernel:
            stamp = folder / 'kernel/library-build-id'
            if not stamp.is_file() or stamp.read_text().strip() != snapshot.identity:
                raise SystemExit('Kernel dependencies do not match this library build. '
                                 'Run `pixi run --locked python scripts/deps.py --output build --tests-only` first.')
            names = [Path(p).name for spec in test_dependencies().values() for p in spec['packages']]
            for name in names:
                if not (folder / 'kernel' / f'{name}.mojoc').is_file():
                    raise SystemExit('Missing kernel packages. Run `pixi run --locked python scripts/deps.py --output build --tests-only` first.')
        yield snapshot
        # Also catch source edits and changes made outside the managed build commands.
        if validate(folder).identity != snapshot.identity:
            raise SystemExit('The library build changed during testing; results are invalid.')
