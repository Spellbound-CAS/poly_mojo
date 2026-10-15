"""Filter compiler-generated docstrings while preserving authored diagnostics."""
import re

GENERATED = {'__init__', '__deinit__', '__copyinit__', '__moveinit__'}
_DIAGNOSTIC = re.compile(r"^(?P<path>[^:\s]+\.mojo):(?P<line>\d+):\d+: (?:warning|error): (?P<message>.*)$")
_SYMBOL = re.compile(r"public symbol '([^']+)'")


def classify(text, root):
    """Return authored docstring problems and all other compiler warnings."""
    docs, warnings = [], []
    for line in text.splitlines():
        match = _DIAGNOSTIC.match(line)
        if match and 'doc string' in match['message']:
            symbol = _SYMBOL.search(match['message'])
            path = root / match['path']
            lines = path.read_text().splitlines() if path.is_file() else []
            index = int(match['line']) - 1
            source = lines[index] if 0 <= index < len(lines) else ''
            if symbol and symbol[1] in GENERATED and not re.match(r'\s*def\s', source):
                continue
            docs.append(f"{match['path']}:{match['line']}: {match['message']}")
        elif re.search(r'\bwarning:', line):
            warnings.append(line)
    return docs, warnings
