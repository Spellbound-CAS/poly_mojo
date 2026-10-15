"""Documentation inventories and rendering; no MkDocs or CLI dependency.

Three directives: `example` runs a maintained program, `source` shows one,
and `api` renders a public declaration of a package, with its docstring,
from the sources (`mojo_api.py`). The audit requires every name a package
exports from its `__init__.mojo` to have an `api` directive on its
reference page, and every rendered declaration and member a docstring.
"""
import json
import posixpath
import re
from pathlib import Path

from packages import ROOT, exports, present, public_modules
CONTENT = ROOT / 'docs/content'
EXAMPLES = ROOT / 'docs/examples'
DIRECTIVE = re.compile(r'<!--\s*(example|source|api):\s*([^>]*?)\s*-->')
FENCE = re.compile(r'(^```[^\n]*\n.*?^```[ \t]*$)', re.M | re.S)


def data(name):
    return json.loads((ROOT / 'docs' / name).read_text())


def relative(page, destination):
    return posixpath.relpath(destination, posixpath.dirname(page) or '.')


def directives(text):
    """The directives of a page, ignoring any inside fenced code."""
    return DIRECTIVE.findall(''.join(FENCE.split(text)[::2]))


def _example_bar(key, page):
    """The header line of an example: its file name and a download link."""
    source = relative(page, 'downloads/' + key)
    return (f'<div class="example-bar" markdown="span"><span class="example-file">{posixpath.basename(key)}</span> '
            f'[Download]({source}){{ .example-download }}</div>')


def example_markdown(key, page):
    """A complete program as one unit: its file, the code, how to run it and what it prints."""
    spec = data('examples.json')['examples'][key]
    code = (ROOT / key).read_text().rstrip()
    return (f'<div class="example" markdown="1">\n{_example_bar(key, page)}\n\n'
            f'```mojo\n{code}\n```\n\n'
            f'<div class="example-run" markdown="1"><span class="example-label">Run from the repository root</span> '
            f'`pixi run mojo run -I build {key}`</div>\n\n'
            f'<div class="example-output" markdown="1"><span class="example-label">Output</span>\n\n'
            f'```text\n{spec["stdout"].rstrip()}\n```\n\n</div>\n</div>\n')


def source_markdown(key, page):
    """Show a maintained example file that is not run on its own."""
    code = (ROOT / key).read_text().rstrip()
    return f'<div class="example" markdown="1">\n{_example_bar(key, page)}\n\n```mojo\n{code}\n```\n\n</div>\n'


def _api_target(value):
    """The file and the declaration an `api` directive names."""
    parts = value.split()
    assert len(parts) == 2, f'An api directive names a file and a declaration: {value}'
    return parts[0], parts[1]


def _sections_markdown(doc):
    """A docstring as Markdown: its prose, then its sections as lists."""
    from mojo_api import split_doc
    prose, sections = split_doc(doc)
    out = [prose] if prose else []
    labels = dict(Args='Arguments', Parameters='Parameters', Returns='Returns', Raises='Raises')
    for name in ('Parameters', 'Args', 'Returns', 'Raises'):
        if name not in sections:
            continue
        entries = sections[name]
        if len(entries) == 1 and not entries[0][0]:
            out.append(f'**{labels[name]}.** {entries[0][1]}')
            continue
        lines = [f'**{labels[name]}.**', '']
        for key, text in entries:
            lines.append(f'- `{key}`: {text}' if key else f'- {text}')
        out.append('\n'.join(lines))
    return '\n\n'.join(out)


def api_markdown(value, page):
    """A public declaration: its signature and docstring, and, for a struct, each public member."""
    from mojo_api import find
    path, name = _api_target(value)
    item = find(ROOT / path, name)
    assert item is not None, f'api directive names no public declaration: {value}'
    out = [f'### `{item.name}`', '', f'```mojo\n{item.signature}\n```', '', _sections_markdown(item.doc), '',
           f'<p class="api-source">Defined in <code>{path}</code>.</p>', '']
    fields = [m for m in item.members if m.kind == 'var']
    if fields:
        out += ['| Field | Type | |', '|---|---|---|']
        for m in fields:
            sig = m.signature.split(':', 1)[1].strip() if ':' in m.signature else ''
            out.append(f'| `{m.name}` | `{sig}` | {" ".join(m.doc.split())} |')
        out.append('')
    for m in item.members:
        if m.kind == 'var':
            continue
        out += [f'#### `{item.name}.{m.name}`', '', f'```mojo\n{m.signature}\n```', '', _sections_markdown(m.doc), '']
    return '\n'.join(out) + '\n'


def render_markdown(text, page):
    def replace(match):
        kind, value = match.groups()
        if kind == 'example':
            return example_markdown(value, page)
        if kind == 'api':
            return api_markdown(value, page)
        return source_markdown(value, page)
    return ''.join(part if index % 2 else DIRECTIVE.sub(replace, part)
                   for index, part in enumerate(FENCE.split(text)))


def audit_docs():
    import sys
    sys.path.insert(0, str(ROOT / 'scripts'))
    from mojo_api import find, parse
    examples = set()
    documented = set()
    pages = sorted(CONTENT.rglob('*.md'))
    for page in pages:
        text = page.read_text()
        assert text.startswith('# '), f'Missing page title: {page}'
        # Executable Mojo comes from maintained files, never drifting copies.
        assert not re.search(r'^```mojo\s*$', text, re.M), f'Use an example directive: {page}'
        for kind, value in directives(text):
            if kind == 'example':
                examples.add(value)
            elif kind == 'api':
                path, name = _api_target(value)
                item = find(ROOT / path, name)
                assert item is not None, f'api directive names no public declaration: {page}: {value}'
                assert item.doc, f'Public declaration without a docstring: {value}'
                for m in item.members:
                    assert m.doc, f'Public member without a docstring: {value}.{m.name}'
                documented.add((path, name))
            else:
                path = ROOT / value
                assert path.is_file() and path.suffix == '.mojo' and path.is_relative_to(EXAMPLES), \
                    f'Source directive must name a file in docs/examples: {page}: {value}'
    for package in present():
        modules = public_modules(package)
        for name, path in exports(package).items():
            names = {item.name for item in parse(ROOT / path)[1]} if name in modules else {name}
            for declaration in names:
                assert (path, declaration) in documented, f'Exported name without an api directive: {package}.{name}.{declaration} ({path})'
    manifest = data('examples.json')['examples']
    assert examples == set(manifest), 'Unused or missing executable example'
    for name, spec in manifest.items():
        assert (ROOT / name).is_file() and name.endswith('.mojo'), f'Missing example: {name}'
        assert isinstance(spec['stdout'], str) and spec['stdout'].endswith('\n'), f'Bad expected output: {name}'
    return dict(pages=len(pages), examples=len(examples), api=len(documented))
