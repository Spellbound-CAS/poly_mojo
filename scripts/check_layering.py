"""Fail if a package imports something the package graph does not allow, or calls foreign code.

The graph is the one of the architecture note (section 27.4) and of
docs/content/architecture/packages.md. `A: {B, ...}` lists what package A may
import directly, besides `std` and itself. Packages not yet in `src/` are
skipped, so the graph can be written down before the code exists.

The packages are pure Mojo (ground rule PG7): no `external_call`,
no dynamic loading, and nothing from `std.ffi` except `_get_global`, which
makes a process-wide value without calling foreign code.
"""
import io
import re
import tokenize

from packages import GRAPH, ROOT, exports, import_source, present

FOREIGN = frozenset({'external_call', 'DLHandle', 'OwnedDLHandle', 'dlopen', 'dlsym'})
FFI_MODULES = frozenset({'std.ffi', 'sys.ffi', 'ffi'})
FFI_ALLOWED = frozenset({'_get_global'})


def foreign_names(text):
    """The names in `text`, outside strings and comments, that call foreign code."""
    found = []
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type == tokenize.NAME and token.string in FOREIGN:
            found.append((token.start[0], token.string))
    return found


IMPORT = re.compile(r'^\s*(?:from\s+([A-Za-z_.][\w.]*)\s+import\s+(?:\(([^)]*)\)|([^\n]+))'
                    r'|import\s+([^\n]+))', re.M)


def check(root=ROOT):
    """Report forbidden package edges, private imports, and relative escapes."""
    problems = []
    unknown = sorted(p.name for p in (root / 'src').iterdir() if p.is_dir() and p.name not in GRAPH)
    problems += [f'src/{name}: not in the package graph of scripts/packages.py' for name in unknown]
    public = {p: exports(p, root) for p in present(root)}
    for package in present(root):
        allowed = GRAPH[package] | {'std', package}
        for path in sorted((root / 'src' / package).rglob('*.mojo')):
            try:
                source = import_source(path.read_text())
                foreign = foreign_names(path.read_text())
            except (tokenize.TokenError, IndentationError) as error:
                problems.append(f'{path.relative_to(root)}: cannot scan imports: {error}')
                continue
            for line, name in foreign:
                problems.append(f'{path.relative_to(root)}:{line}: {name} calls foreign code; '
                                f'{package} is pure Mojo (PG7)')
            for match in IMPORT.finditer(source):
                modules = [match[1]] if match[1] else [s.strip().split()[0] for s in match[4].split('#')[0].split(',') if s.strip()]
                for module in modules:
                    prefix = f'{path.relative_to(root)}: {package} may not import {module}'
                    if module.startswith('.'):
                        dots = len(module) - len(module.lstrip('.'))
                        if dots > len(path.relative_to(root / 'src').parent.parts):
                            problems.append(prefix + ' (relative import leaves its package)')
                        continue
                    if module in FFI_MODULES:
                        if match[1]:
                            names = ' '.join(line.split('#')[0] for line in (match[2] or match[3]).splitlines())
                            imported = {item.strip().split()[0] for item in names.split(',') if item.strip()}
                        else:
                            imported = {'*'}
                        for name in sorted(imported - FFI_ALLOWED):
                            problems.append(prefix + f': {name} (pure Mojo, PG7)')
                        continue
                    dependency = module.split('.')[0]
                    if dependency not in allowed:
                        problems.append(prefix)
                    elif dependency in GRAPH and dependency != package:
                        if module != dependency:
                            problems.append(prefix + ' (use the package exports)')
                        elif match[1] and dependency in public:
                            names = ' '.join(line.split('#')[0] for line in (match[2] or match[3]).splitlines())
                            for item in names.split(','):
                                name = item.strip().split()[0] if item.strip() else ''
                                if name and name not in public[dependency]:
                                    problems.append(prefix + f': {name} is not exported')
    return problems


def main():
    problems = check()
    if problems:
        print('\n'.join(problems))
        return 1
    print('layering ok:', ', '.join(f'{p} -> {", ".join(sorted(GRAPH[p])) or "std"}' for p in present()))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
