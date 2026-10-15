"""Compile a suite once and dispatch individual programs from the cached executable."""
from dataclasses import dataclass
import fcntl
import hashlib
import io
import tokenize
import json
from pathlib import Path
import subprocess
import time

from build_state import digest


def run(command, root, timeout=None):
    started = time.monotonic()
    try:
        result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=timeout)
        return dict(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr,
                    seconds=round(time.monotonic() - started, 3))
    except subprocess.TimeoutExpired:
        return dict(returncode=None, stdout='', stderr=f'timeout after {timeout} seconds',
                    seconds=round(time.monotonic() - started, 3))


def source_inputs(root, programs):
    """Include sibling modules and shared helpers used by the selected programs."""
    paths = set(programs)
    for parent in {path.parent for path in programs}:
        paths.update(parent.rglob('*.mojo'))
    paths.update((root / 'tests/support').rglob('*.mojo'))
    return {str(path.relative_to(root)): digest(path) for path in sorted(paths)}


def render(programs, root):
    """Import each main under a distinct name; dispatch without rebuilding."""
    lines = ['from std.sys import argv']
    for i, path in enumerate(programs):
        module = '.'.join(path.relative_to(root).with_suffix('').parts)
        lines.append(f'from _test_cases.{module} import _suite_entry_point as _case_{i}')
    lines += ['', 'def main() raises:', '    var args = argv()',
              '    if len(args) != 2:', '        raise Error("expected a test index")',
              '    var selected = Int(args[1])']
    for i in range(len(programs)):
        lines += [f'    if selected == {i}:', f'        _case_{i}()', '        return']
    lines.append('    raise Error("unknown test index")')
    return '\n'.join(lines) + '\n'


def module_source(source):
    """Rename the actual top-level entry point, leaving comments/strings intact."""
    tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    entries = [tokens[i + 1] for i, token in enumerate(tokens[:-1])
               if token.type == tokenize.NAME and token.string == 'def' and token.start[1] == 0
               and tokens[i + 1].type == tokenize.NAME and tokens[i + 1].string == 'main']
    if len(entries) != 1:
        raise ValueError('A runtime program must define exactly one top-level main')
    entry = entries[0]
    lines = source.splitlines(keepends=True)
    row, column = entry.start
    lines[row - 1] = lines[row - 1][:column] + '_suite_entry_point' + lines[row - 1][entry.end[1]:]
    return ''.join(lines)


@dataclass
class Executable:
    binary: Path
    programs: tuple
    record: dict
    reused: bool

    def report(self):
        return {**self.record, 'binary': str(self.binary), 'reused': self.reused}

    def execute(self, path, root, timeout=None):
        return run([str(self.binary), str(self.programs.index(path))], root, timeout)


class CompilationFailed(Exception):
    def __init__(self, record):
        self.record = record
        super().__init__(record['compile']['stderr'] or
                         f"suite compiler exited with status {record['compile']['returncode']}")


def prepare(programs, snapshot, root, *, flags=(), cache=None, kernel=False, timeout=None, compiler_threads=2):
    """Produce at most one executable for this complete source/configuration snapshot."""
    if compiler_threads < 1:
        raise ValueError('compiler_threads must be positive')
    root = root.resolve()
    programs = tuple(sorted(set(path.resolve() for path in programs)))
    if not programs:
        raise ValueError('Cannot compile an empty suite')
    includes = [snapshot.path]
    if kernel:
        includes.append(snapshot.path / 'kernel')
    inputs = {
        'schema': 1, 'library_build': snapshot.report(), 'kernel': kernel,
        'programs': [str(path.relative_to(root)) for path in programs],
        'sources': source_inputs(root, programs), 'flags': list(flags),
        'compiler_threads': compiler_threads,
        'kernel_artifacts': {str(p): digest(p) for p in sorted((snapshot.path / 'kernel').glob('*.mojoc'))} if kernel else {},
        'generator': digest(Path(__file__)),
    }
    key = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
    folder = (cache or root / '.cache/tests/executables') / key
    folder.mkdir(parents=True, exist_ok=True)
    binary = folder / 'runtime'
    manifest = folder / 'manifest.json'
    with (folder / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            record = json.loads(manifest.read_text())
            if record['inputs'] == inputs and record['binary_sha256'] == digest(binary):
                return Executable(binary, programs, record, True)
        except (OSError, ValueError, KeyError):
            pass
        manifest.unlink(missing_ok=True)
        source = folder / 'suite.mojo'
        source.write_text(render(programs, root))
        package = folder / '_test_cases'
        for path in programs:
            target = package / path.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            parent = target.parent
            while parent != folder:
                (parent / '__init__.mojo').write_text('')
                parent = parent.parent
            target.write_text(module_source(path.read_text()))
        includes += [folder, root / 'tests', *sorted({p.parent for p in programs})]
        command = ['mojo', 'build', '--num-threads', str(compiler_threads), *flags]
        for include in dict.fromkeys(includes):
            command += ['-I', str(include)]
        pending = folder / 'runtime.tmp'
        command += [str(source), '-o', str(pending)]
        print(f'Compiling one executable for {len(programs)} runtime programs', flush=True)
        built = run(command, root, timeout)
        record = {'inputs': inputs, 'compile': built, 'command': command}
        if built['returncode'] != 0:
            pending.unlink(missing_ok=True)
            raise CompilationFailed(record)
        if source_inputs(root, programs) != inputs['sources']:
            pending.unlink(missing_ok=True)
            raise RuntimeError('Test sources changed during compilation; result discarded')
        pending.replace(binary)
        record['binary_sha256'] = digest(binary)
        manifest.write_text(json.dumps(record, indent=2) + '\n')
        return Executable(binary, programs, record, False)
