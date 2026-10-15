"""Read the public API of Mojo sources: structs, traits, functions, and constants, with their docstrings.

The documentation's `api` directive (docs_contract.py) renders reference
pages from what this module reads. It understands the subset of Mojo that
the packages use: top-level `struct`, `trait`, `def`, and `comptime` declarations,
preceded by decorators; inside a struct, `var` fields, `comptime` members,
and `def` methods; and a docstring right after each declaration, in the
Google style of the sources (a summary, paragraphs, and `Args:`,
`Parameters:`, `Returns:`, and `Raises:` sections). Names that start with
`_` are private, and so are dunder methods other than `__init__`.
"""
from dataclasses import dataclass, field
from pathlib import Path
import re
import textwrap

SECTIONS = ('Args', 'Parameters', 'Returns', 'Raises')
_DECL = re.compile(r'^(struct|trait|def|comptime|var)\s+([A-Za-z_][A-Za-z0-9_]*)')


@dataclass
class Item:
    """A declaration and its docstring."""

    kind: str
    name: str
    signature: str
    doc: str
    line: int
    members: list = field(default_factory=list)


def _docstring(lines, i):
    """The docstring that starts at line `i`, if any, and the line after it."""
    if i >= len(lines) or not lines[i].strip().startswith('"""'):
        return '', i
    first = lines[i].strip()[3:]
    if '"""' in first:
        return first[:first.index('"""')].strip(), i + 1
    parts = [first]
    k = i + 1
    while k < len(lines) and '"""' not in lines[k]:
        parts.append(lines[k])
        k += 1
    if k < len(lines):
        parts.append(lines[k][:lines[k].index('"""')])
    head = parts[0]
    rest = textwrap.dedent('\n'.join(parts[1:]))
    return (head + '\n' + rest).strip(), k + 1


def _depth(text):
    """The bracket depth at the end of a text, ignoring strings."""
    depth = 0
    quote = None
    for ch in text:
        if quote:
            if ch == quote:
                quote = None
        elif ch in '"\'':
            quote = ch
        elif ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
    return depth


def _signature(lines, i, kind):
    """A declaration, joined across lines up to the colon that ends it (or the line, for `var` and `comptime`)."""
    text = lines[i].strip()
    k = i
    if kind in ('var', 'comptime'):
        while _depth(text) > 0 and k + 1 < len(lines):
            k += 1
            text += ' ' + lines[k].strip()
        return text, k + 1
    while not (_depth(text) == 0 and text.rstrip().endswith(':')) and k + 1 < len(lines):
        k += 1
        text += ' ' + lines[k].strip()
    return re.sub(r'\s+', ' ', text.rstrip()[:-1]).replace('( ', '(').replace(' )', ')').replace(', )', ')'), k + 1


def _public(name):
    """Whether a name is part of the public API."""
    if name == '__init__':
        return True
    return not name.startswith('_')


def parse(path):
    """The module docstring and the public top-level items of a Mojo source file."""
    lines = Path(path).read_text().splitlines()
    module_doc, i = _docstring(lines, 0)
    items = []
    decorators = []
    current = None
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())
        if stripped.startswith('@'):
            decorators.append(stripped)
            i += 1
            continue
        m = _DECL.match(stripped)
        if not m or indent not in (0, 4):
            decorators = []
            i += 1
            continue
        kind, name = m.groups()
        if indent == 0:
            current = None
        elif current is None:
            decorators = []
            i += 1
            continue
        signature, i = _signature(lines, i, kind)
        doc, i = _docstring(lines, i)
        item = Item(kind, name, '\n'.join(decorators + [signature]), doc, i)
        decorators = []
        if indent == 0:
            if kind == 'var' or not _public(name):
                continue
            items.append(item)
            if kind in ('struct', 'trait'):
                current = item
        elif _public(name):
            current.members.append(item)
    return module_doc, items


def find(path, name):
    """A public top-level item of a file, by name."""
    for item in parse(path)[1]:
        if item.name == name:
            return item
    return None


def split_doc(doc):
    """Split a docstring into its prose and its sections, each a list of (name, text) entries or one text."""
    prose = []
    sections = {}
    section = None
    for line in doc.splitlines():
        head = line.strip().rstrip(':')
        if line.strip().endswith(':') and head in SECTIONS and not line.startswith(' '):
            section = head
            sections[section] = []
            continue
        if section is None:
            prose.append(line)
        else:
            sections[section].append(line)
    out = {}
    for name, body in sections.items():
        entries = []
        for line in body:
            if not line.strip():
                continue
            if line.startswith('    ') and not line.startswith('        '):
                text = line.strip()
                m = re.match(r'([A-Za-z_][A-Za-z0-9_]*|\*?[a-z_]+):\s*(.*)', text)
                if m and name in ('Args', 'Parameters'):
                    entries.append([m[1], m[2]])
                elif m and name == 'Raises' and m[1][0].isupper():
                    entries.append([m[1], m[2]])
                elif entries:
                    entries[-1][1] += ' ' + text
                else:
                    entries.append(['', text])
            elif entries:
                entries[-1][1] += ' ' + line.strip()
            else:
                entries.append(['', line.strip()])
        out[name] = entries
    return '\n'.join(prose).strip(), out
