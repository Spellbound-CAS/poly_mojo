"""Interval gate for CAD answers: a new answer may not be coarser than a baseline.

Usage: python -m harness.cad_gate BASELINE NEW [CASE ...]

BASELINE and NEW are answer files written by the Poly Mojo CAD driver (one
`dim signs coords` line per terminal cell) or directories of `<case>.txt`
answers, in which case every case present in both is compared. A pair passes
when both list the same cells in the same order, with the same dimensions,
signs, and every algebraic coordinate's isolating interval in NEW is no wider
than the baseline's; an exact rational is never wider, while a baseline
rational that becomes an interval fails. Sample values may differ:
canonical sector samples move the fibres above them.
"""
from fractions import Fraction
from pathlib import Path
import sys


def coordinate(text):
    if text.startswith('Q:'):
        return ('Q', None)
    _, _, low, high = text.split(':')
    return ('A', Fraction(high) - Fraction(low))


def cells(path):
    lines = Path(path).read_text().splitlines()
    rows = []
    for line in lines[1:]:
        if line:
            parts = line.split(' ')
            rows.append((parts[0], parts[1], [coordinate(t) for t in parts[2:]]))
    return int(lines[0]), rows


def compare(baseline, new):
    """Failures as strings; empty when NEW passes against BASELINE."""
    count, old_rows = cells(baseline)
    new_count, new_rows = cells(new)
    if count != new_count or len(old_rows) != len(new_rows):
        return [f'cell count {count} != {new_count}']
    failures = []
    for index, (old, cur) in enumerate(zip(old_rows, new_rows)):
        if old[0] != cur[0] or old[1] != cur[1]:
            failures.append(f'cell {index}: dimension or signs differ')
            continue
        if len(old[2]) != len(cur[2]):
            failures.append(f'cell {index}: coordinate count differs')
            continue
        for axis, (a, b) in enumerate(zip(old[2], cur[2])):
            if a[0] == 'Q' and b[0] == 'A':
                failures.append(f'cell {index} axis {axis}: exact baseline coordinate is now an interval')
            elif a[0] == 'A' and b[0] == 'A' and b[1] > a[1]:
                failures.append(f'cell {index} axis {axis}: interval wider than baseline')
    return failures


def main(argv):
    baseline, new = Path(argv[0]), Path(argv[1])
    if baseline.is_dir():
        names = argv[2:] or sorted(p.stem for p in baseline.glob('*.txt') if (new / p.name).exists())
        pairs = [(name, baseline / f'{name}.txt', new / f'{name}.txt') for name in names]
    else:
        pairs = [(new.stem, baseline, new)]
    passed = True
    for name, old, cur in pairs:
        failures = compare(old, cur)
        passed = passed and not failures
        print(f'{name}: {"pass" if not failures else "FAIL " + "; ".join(failures[:5])}')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
