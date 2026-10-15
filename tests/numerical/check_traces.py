"""Independent standard-library verifier for trace_oracle.mojo output.

Run the Mojo fixture, then pass its stdout file as this program's argument.
No imports from the polynomial implementation or apn_mojo are used.
"""
from fractions import Fraction as F
from pathlib import Path
import sys


def contains(a, b):
    return all(a[k] <= b[k] <= b[k + 1] <= a[k + 1] for k in (0, 2))


def product(a, b):
    candidates = [x * y for x in a for y in b]
    return min(candidates), max(candidates)


def mul(a, b):
    rr, ii = product(a[:2], b[:2]), product(a[2:], b[2:])
    ri, ir = product(a[:2], b[2:]), product(a[2:], b[:2])
    return rr[0] - ii[1], rr[1] - ii[0], ri[0] + ir[0], ri[1] + ir[1]


def verify(lines):
    rows = iter(lines)
    checked = steps = 0
    # Independently specified source polynomials and exactly dyadic arguments.
    fixtures = [
        ([F(1, 3), F(-2), F(1)], (F(11, 8), F(13, 8), F(0), F(0))),
        ([F(1, 3), F(-2), F(1)], (F(-1, 16), F(1, 16), F(15, 16), F(17, 16))),
        ([F(1), F(-4), F(6), F(-4), F(1)], (F(1023, 1024), F(1025, 1024), F(0), F(0))),
    ]
    for header in rows:
        tag, ni, ns, final = header.strip().split(',')
        assert tag == 'trace'
        inputs = [tuple(map(F, next(rows).split(','))) for _ in range(int(ni))]
        coefficients, argument = fixtures[checked]
        expected_inputs = [(x, x, F(0), F(0)) for x in coefficients]
        expected_inputs += [argument, (F(0),) * 4]
        assert inputs == expected_inputs
        instructions = [(0, i, 0) for i in range(len(expected_inputs))]
        result = len(coefficients) + 1
        for coefficient in reversed(range(len(coefficients))):
            instructions.append((2, result, len(coefficients)))
            instructions.append((1, len(instructions) - 1, coefficient))
            result = len(instructions) - 1
        assert (int(ns), int(final)) == (len(instructions), result)
        bounds = []
        for expected_instruction in instructions:
            op, left, right = map(int, next(rows).split(','))
            assert (op, left, right) == expected_instruction
            box = tuple(map(F, next(rows).split(',')))
            expected = (inputs[left] if op == 0 else
                        tuple(x + y for x, y in zip(bounds[left], bounds[right])) if op == 1 else
                        mul(bounds[left], bounds[right]))
            assert contains(box, expected), (op, box, expected)
            bounds.append(box)
            steps += 1
        value = tuple(map(F, next(rows).split(',')))
        assert contains(value, bounds[int(final)])
        checked += 1
    return checked, steps


if __name__ == '__main__':
    lines = Path(sys.argv[1]).read_text().splitlines()
    checked, steps = verify(lines)
    assert checked == 3
    # A deliberately impossible final enclosure must fail this independent path.
    changed = lines.copy()
    changed[-1] = '100,101,0,0'
    try:
        verify(changed)
    except AssertionError:
        pass
    else:
        raise AssertionError('corrupted trace accepted')
    print(f'{checked} certificates / {steps} steps verified with Python Fraction; corruption rejected')
