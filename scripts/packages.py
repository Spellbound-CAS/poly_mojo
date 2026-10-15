"""Repository package metadata, shared by build, import, and documentation gates.

Dependencies are direct imports, as in architecture/overview.md; external
packages (apn_mojo, and any test-only packages) are fetched by scripts/deps.py.
Repository packages expose only the names exported by their __init__.mojo;
internal modules remain accessible to their owning package.
"""
from pathlib import Path
import io
import json
import re
import tokenize

ROOT = Path(__file__).resolve().parent.parent
GRAPH = {
    'polynomials': frozenset({'apn_mojo'}),
    'poly': frozenset({'polynomials', 'apn_mojo'}),
}
# Insertion order is build order; every repository dependency precedes its user.
PACKAGES = tuple(GRAPH)
# An optional kernel suite in tests/kernel/ declares its own test-only
# dependencies here, in the format of scripts/deps.py's DEPS entries.
TEST_DEPS = Path('tests/kernel/deps.json')
_EXPORT = re.compile(r'^\s*from \.(\w+) import (?:\(([^)]*)\)|([^(\n]+)$)', re.M)


def import_source(text):
    """Canonicalize import statements from Mojo's lexical tokens.

    Ignore strings and comments. Dotted paths, explicit continuations, and
    multiple statements per line all get the same spelling before scanning.
    Tokenization understands these lexical forms without parsing Mojo types.
    """
    tokens = list(tokenize.generate_tokens(io.StringIO(text).readline))
    statements = []
    i = 0
    ignored = {tokenize.COMMENT, tokenize.STRING, tokenize.NL, tokenize.INDENT, tokenize.DEDENT}
    while i < len(tokens):
        token = tokens[i]
        if token.type != tokenize.NAME or token.string not in {'from', 'import'}:
            i += 1
            continue
        kind = token.string
        i += 1
        parts = []
        while i < len(tokens) and tokens[i].type not in {tokenize.NEWLINE, tokenize.ENDMARKER} and tokens[i].string != ';':
            if tokens[i].type not in ignored:
                parts.append(tokens[i].string)
            i += 1
        if kind == 'from' and 'import' in parts:
            split = parts.index('import')
            statements.append('from ' + ''.join(parts[:split]) + ' import ' + ' '.join(parts[split + 1:]))
        elif kind == 'import':
            modules = []
            for group in ' '.join(parts).split(','):
                module, separator, alias = group.partition(' as ')
                modules.append(module.replace(' ', '') + (separator + alias.strip() if separator else ''))
            statements.append('import ' + ', '.join(modules))
    return '\n'.join(statements) + '\n'


def test_dependencies(root=ROOT):
    """Test-only dependencies declared by the optional kernel suite; none without it."""
    path = root / TEST_DEPS
    if not path.is_file():
        return {}
    return {name: dict(spec, test_only=True) for name, spec in json.loads(path.read_text()).items()}


def present(root=ROOT):
    """Packages with an entry point, in build order, without requiring Mojo."""
    return [p for p in PACKAGES if (root / 'src' / p / '__init__.mojo').is_file()]


def exports(package, root=ROOT):
    """Exported names and their defining source paths, read from the entry point."""
    text = import_source((root / 'src' / package / '__init__.mojo').read_text())
    out = public_modules(package, root)
    for module, grouped, single in _EXPORT.findall(text):
        for item in (grouped or single).replace('\n', ' ').split(','):
            name = item.strip()
            if name:
                out[name] = f'src/{package}/{module}.mojo'
    return out


def public_modules(package, root=ROOT):
    """Submodules explicitly exported as namespaces from a package root."""
    text = import_source((root / 'src' / package / '__init__.mojo').read_text())
    out = {}
    for statement in re.findall(r'^from \. import (.*)$', text, re.M):
        for item in statement.strip('() ').split(','):
            name = item.strip()
            if name:
                out[name] = f'src/{package}/{name}.mojo'
    return out
