"""Reference implementation of the kernel's exact literal kinds.

It is deliberately simple and independent of apn_mojo: Python integers and
fractions. The kernel's literal operations are tested against it.

Each operation takes the kind, `INTEGER` or `RATIONAL` (the rational kind by
default), and returns the result, or None when the kernel must decline. An
operation evaluates exactly where the literal laws determine its value, and
declines everywhere else. Addition and
multiplication never decline: when the actual result's encoding exceeds
`MAX_LITERAL_BYTES`, they raise `CapacityError`, as the kernel raises
`CAPACITY`. Run it from the repository root:

    python3 tests/reference/kernel/literal_reference.py           # print the encodings, check the vector file
    python3 tests/reference/kernel/literal_reference.py --write   # regenerate the vector file

The vector file, literal_vectors.txt next to this script, drives the kernel's
differential tests: values drawn from a fixed seed by SplitMix64 (each
draw adds GOLDEN to the state and returns its finalizer), with small values,
values around word boundaries, signs, zero, one, values up to 75 bytes,
and the declining cases. It is plain text, one record per line, fields
separated by one space; a value is written in decimal, `N` or `N/D`:

    literal-vectors SEED
    enc KIND VALUE HEX               the encoding of VALUE
    bad KIND HEX                     an encoding decoding must reject (ENCODING)
    op KIND OP A B RESULT            OP is add, mul, lt, le, eq, or pow
    op KIND OP A - RESULT            OP is neg or inv
    fold KIND OP N V1 ... VN RESULT  OP is add or mul

RESULT is a value, `True` or `False`, `DECLINE`, or `CAPACITY`. KIND is
`integer` or `rational`.
"""
import sys
from fractions import Fraction
from pathlib import Path


def int_bytes(n: int) -> bytes:
    """Sign byte (0 non-negative, 1 negative), then |n| little-endian, no trailing zero byte."""
    magnitude = abs(n)
    body = magnitude.to_bytes((magnitude.bit_length() + 7) // 8, 'little') if magnitude else b''
    return bytes([1 if n < 0 else 0]) + body


def int_from_bytes(data: bytes) -> int:
    if not data or data[0] not in (0, 1) or (len(data) > 1 and data[-1] == 0):
        raise ValueError('not a canonical integer encoding')
    if data[0] == 1 and len(data) == 1:
        raise ValueError('negative zero')
    magnitude = int.from_bytes(data[1:], 'little')
    return -magnitude if data[0] else magnitude


def rational_bytes(q: Fraction) -> bytes:
    """4-byte little-endian length of the numerator's encoding, the numerator, then the denominator's magnitude."""
    num = int_bytes(q.numerator)
    den = q.denominator.to_bytes((q.denominator.bit_length() + 7) // 8, 'little')
    return len(num).to_bytes(4, 'little') + num + den


def rational_from_bytes(data: bytes) -> Fraction:
    if len(data) < 5:
        raise ValueError('too short')
    length = int.from_bytes(data[:4], 'little')
    num = int_from_bytes(data[4:4 + length])
    den_bytes = data[4 + length:]
    if not den_bytes or den_bytes[-1] == 0:
        raise ValueError('not a canonical denominator')
    den = int.from_bytes(den_bytes, 'little')
    q = Fraction(num, den)
    if q.numerator != num or q.denominator != den:
        raise ValueError('not reduced')
    return q


# The two literal kinds of release 1 (EXACT_INTEGER and EXACT_RATIONAL).
INTEGER = 'integer'
RATIONAL = 'rational'

MAX_LITERAL_BYTES = 2 ** 20  # the ceiling on a literal's encoding
POWER_LIMIT = 4096           # the largest exponent magnitude the kernel evaluates


class CapacityError(Exception):
    """A result's encoding exceeds MAX_LITERAL_BYTES: the kernel raises CAPACITY."""


def encode(value, kind: str) -> bytes:
    """The payload of a literal of the kind (the values 0 and 1 are never literal atoms)."""
    value = _value(value, kind)
    return int_bytes(int(value)) if kind == INTEGER else rational_bytes(value)


def encoded_length(value, kind: str) -> int:
    """The length the kernel reads for an operand: 1 for 0 and 1, which have no pool entry."""
    value = _value(value, kind)
    if value in (0, 1):
        return 1
    return len(encode(value, kind))


def _value(a, kind: str) -> Fraction:
    """Validate that a is a value of the kind."""
    if kind not in (INTEGER, RATIONAL):
        raise ValueError(f'unknown literal kind {kind!r}')
    q = Fraction(a)
    if kind == INTEGER and q.denominator != 1:
        raise ValueError(f'{a} is not a value of the integer kind')
    return q


def _result(q: Fraction, kind: str):
    """A value in the kind's Python type: int for the integer kind, Fraction for the rational kind."""
    return int(q) if kind == INTEGER else q


def _within_ceiling(q: Fraction, kind: str):
    if encoded_length(q, kind) > MAX_LITERAL_BYTES:
        raise CapacityError(f'result of {encoded_length(q, kind)} bytes exceeds {MAX_LITERAL_BYTES}')
    return _result(q, kind)


def _operand(a, kind: str) -> Fraction:
    """An operand is a literal, so its encoding is within the ceiling (mk_lit raises CAPACITY otherwise)."""
    q = _value(a, kind)
    if encoded_length(q, kind) > MAX_LITERAL_BYTES:
        raise CapacityError('operand exceeds the ceiling')
    return q


def add(a, b, kind: str = RATIONAL):
    """a + b; never declines, raises CapacityError when the actual result is too large."""
    return _within_ceiling(_operand(a, kind) + _operand(b, kind), kind)


def mul(a, b, kind: str = RATIONAL):
    """a * b; never declines, raises CapacityError when the actual result is too large."""
    return _within_ceiling(_operand(a, kind) * _operand(b, kind), kind)


def neg(a, kind: str = RATIONAL):
    """-a; its encoding is as long as a's, so it cannot exceed the ceiling."""
    return _result(-_operand(a, kind), kind)


def inv(a, kind: str = RATIONAL):
    """1/a for the rational kind; declines on 0, where the field laws are silent."""
    if kind != RATIONAL:
        raise ValueError('INV is bound only by the rational kind')
    q = _operand(a, kind)
    return None if q == 0 else 1 / q


def lt(a, b, kind: str = RATIONAL): return _operand(a, kind) < _operand(b, kind)
def le(a, b, kind: str = RATIONAL): return _operand(a, kind) <= _operand(b, kind)
def eq(a, b, kind: str = RATIONAL): return _operand(a, kind) == _operand(b, kind)


def power(a, n, kind: str = RATIONAL):
    """a**n, or None where the kernel declines.

    It evaluates for an integer n with |n| <= POWER_LIMIT when a != 0, or
    a == 0 and n > 0, and, for the integer kind, only n >= 0. It also
    declines, before computing, when |n| times the base's encoded length
    exceeds the ceiling; the bound is never below the result's length.
    """
    base = _operand(a, kind)
    exponent = _operand(n, kind)
    if exponent.denominator != 1 or abs(exponent) > POWER_LIMIT:
        return None
    n = int(exponent)
    if base == 0 and n <= 0:
        return None
    if kind == INTEGER and n < 0:
        return None
    if n != 0 and abs(n) * encoded_length(base, kind) > MAX_LITERAL_BYTES:
        return None
    result = base ** n
    assert encoded_length(result, kind) <= MAX_LITERAL_BYTES
    return _result(result, kind)


ADD, MUL = 'ADD', 'MUL'


def fold(op: str, values, kind: str = RATIONAL):
    """Combine the literal arguments of an ADD or MUL spine, in the order given.

    Canonical construction folds them in the spine's canonical order; each intermediate result is checked against the ceiling,
    so the first one that exceeds it raises CapacityError. With no
    values, the result is the unit.
    """
    if op == ADD:
        step, acc = add, Fraction(0)
    elif op == MUL:
        step, acc = mul, Fraction(1)
    else:
        raise ValueError(f'only ADD and MUL spines fold, not {op!r}')
    for v in values:
        acc = Fraction(step(acc, v, kind))
    return _result(acc, kind)


def _raises_capacity(f, *args) -> bool:
    try:
        f(*args)
    except CapacityError:
        return True
    return False


def self_test() -> None:
    """The declining cases of powers and the encoding ceiling; silent when they hold."""
    for kind in (INTEGER, RATIONAL):
        for v in (2, -1, -3, 255, 256, -(2 ** 64), 2 ** 64 - 1):
            data = encode(v, kind)
            back = int_from_bytes(data) if kind == INTEGER else rational_from_bytes(data)
            assert back == v
        assert add(2, 3, kind) == 5 and mul(-2, 3, kind) == -6 and neg(4, kind) == -4
        assert lt(2, 3, kind) and le(3, 3, kind) and not eq(2, 3, kind)
        # Powers.
        assert power(0, 0, kind) is None
        assert power(0, -1, kind) is None
        assert power(0, 3, kind) == 0
        assert power(5, 0, kind) == 1
        assert power(-2, 3, kind) == -8
        assert power(2, POWER_LIMIT, kind) == 2 ** POWER_LIMIT
        assert power(2, POWER_LIMIT + 1, kind) is None
        assert power(2, -(POWER_LIMIT + 1), kind) is None
    assert power(2, -1, INTEGER) is None
    assert power(2, -3, RATIONAL) == Fraction(1, 8)
    assert power(Fraction(-2, 3), -2, RATIONAL) == Fraction(9, 4)
    assert power(4, Fraction(1, 2), RATIONAL) is None
    assert inv(0) is None and inv(Fraction(-3, 4)) == Fraction(-4, 3)
    try:
        inv(2, INTEGER)
        raise AssertionError('INV accepted for the integer kind')
    except ValueError:
        pass
    try:
        add(Fraction(1, 2), 1, INTEGER)
        raise AssertionError('a rational accepted by the integer kind')
    except ValueError:
        pass
    # The power bound: |n| times the base's length, checked before computing.
    base = 2 ** (8 * 1024 - 1)                    # 1,025 bytes as an integer literal
    assert encoded_length(base, INTEGER) == 1025
    assert power(base, 1023, INTEGER) is not None   # 1,048,575 bytes of bound
    assert power(base, 1024, INTEGER) is None       # 1,049,600 bytes of bound
    # Addition and multiplication check the actual result.
    a = 2 ** (8 * 600 * 1024)                     # a 600 KB literal
    assert add(a, 2, INTEGER) == a + 2 and add(a, 2, RATIONAL) == a + 2
    b = 2 ** (8 * 900 * 1024)                     # above 256 KB, and half the ceiling
    assert add(b, -b, INTEGER) == 0 and fold(ADD, [b, -b], RATIONAL) == 0
    top = 2 ** (8 * (MAX_LITERAL_BYTES - 1)) - 1  # an integer literal exactly at the ceiling
    assert encoded_length(top, INTEGER) == MAX_LITERAL_BYTES
    assert add(top, -1, INTEGER) == top - 1
    assert _raises_capacity(add, top, 1, INTEGER)
    half = 2 ** (8 * 524288 - 1)                  # its square has 2^20 magnitude bytes
    assert _raises_capacity(mul, half, half, INTEGER)
    assert mul(half, Fraction(1, half), RATIONAL) == 1
    # A fold checks every intermediate result, so the order matters.
    assert _raises_capacity(fold, ADD, [top, top, -top], INTEGER)
    assert fold(ADD, [top, -top, top], INTEGER) == top
    assert fold(MUL, [], RATIONAL) == 1 and fold(ADD, [], INTEGER) == 0
    assert fold(MUL, [2, Fraction(1, 2), 3], RATIONAL) == 3


VECTOR_FILE = Path(__file__).with_name('literal_vectors.txt')
SEED = 0x6C69746572616C73  # 'literals'
MASK = (1 << 64) - 1
GOLDEN = 0x9E3779B97F4A7C15


class SplitMix64:
    """Vigna's SplitMix64: each draw adds GOLDEN to the state and returns its finalizer."""

    def __init__(self, seed: int):
        self.state = seed & MASK

    def next(self) -> int:
        self.state = (self.state + GOLDEN) & MASK
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
        return z ^ (z >> 31)

    def below(self, n: int) -> int:
        return self.next() % n

    def bits(self, n: int) -> int:
        """A random integer of at most n bits."""
        value = 0
        for _ in range((n + 63) // 64):
            value = (value << 64) | self.next()
        return value & ((1 << n) - 1)


def _text(q) -> str:
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else f'{q.numerator}/{q.denominator}'


def _special_integers():
    edges = []
    for k in (7, 8, 15, 16, 31, 32, 63, 64, 127, 128, 255, 256):
        edges += [2 ** k - 1, 2 ** k, 2 ** k + 1]
    return [0, 1, -1, 2, -2, 3, -3] + edges + [-e for e in edges]


def _random_integer(rng: SplitMix64) -> int:
    size = (1, 8, 63, 64, 65, 128, 129, 600)[rng.below(8)]
    value = rng.bits(1 + rng.below(size))
    return -value if rng.below(2) else value


def _values(rng: SplitMix64, kind: str, count: int):
    values = list(_special_integers())
    if kind == RATIONAL:
        values += [Fraction(1, 2), Fraction(-1, 2), Fraction(-3, 4), Fraction(2 ** 64, 3), Fraction(-1, 2 ** 64 + 1)]
    while len(values) < count:
        n = _random_integer(rng)
        if kind == RATIONAL and rng.below(3):
            d = 0
            while d == 0:
                d = abs(_random_integer(rng))
            values.append(Fraction(n, d))
        else:
            values.append(n)
    return [Fraction(v) for v in values]


def _outcome(f, *args) -> str:
    try:
        r = f(*args)
    except CapacityError:
        return 'CAPACITY'
    if r is None:
        return 'DECLINE'
    if isinstance(r, bool):
        return 'True' if r else 'False'
    return _text(r)


def _bad_encodings(kind: str):
    """Encodings that are not canonical: each must fail to decode."""
    if kind == INTEGER:
        # empty; sign byte 2; negative zero; trailing zero bytes
        return ['', '02', '01', '0000', '000100', '010500']
    # too short; numerator length 0; no denominator; trailing zero in the
    # denominator; not reduced (2/4 and 0/2); numerator length past the end;
    # a negative-zero numerator
    return ['', '01000000', '0000000000', '0100000000', '010000000000', '0200000000020400',
            '02000000000204', '010000000002', '0500000000020102', '010000000101']


def vector_lines():
    rng = SplitMix64(SEED)
    lines = [f'literal-vectors {SEED:#x}']
    for kind in (INTEGER, RATIONAL):
        values = _values(rng, kind, 120)
        for v in values:
            lines.append(f'enc {kind} {_text(v)} {encode(v, kind).hex()}')
        for b in _bad_encodings(kind):
            data = bytes.fromhex(b)
            try:
                int_from_bytes(data) if kind == INTEGER else rational_from_bytes(data)
                raise AssertionError(f'{b} decodes')
            except ValueError:
                pass
            except ZeroDivisionError:
                pass
            lines.append(f'bad {kind} {b if b else "-"}')
        for _ in range(300):
            a, b = values[rng.below(len(values))], values[rng.below(len(values))]
            for name, f in (('add', add), ('mul', mul), ('lt', lt), ('le', le), ('eq', eq)):
                lines.append(f'op {kind} {name} {_text(a)} {_text(b)} {_outcome(f, a, b, kind)}')
        for v in values:
            lines.append(f'op {kind} neg {_text(v)} - {_outcome(neg, v, kind)}')
            if kind == RATIONAL:
                lines.append(f'op {kind} inv {_text(v)} - {_outcome(inv, v, kind)}')
        exponents = [0, 1, 2, 3, 5, -1, -2, -3, POWER_LIMIT, POWER_LIMIT + 1, -POWER_LIMIT, -(POWER_LIMIT + 1)]
        if kind == RATIONAL:
            exponents += [Fraction(1, 2), Fraction(-3, 2)]
        bases = [Fraction(b) for b in (0, 1, -1, 2, -2, 3, 255, 2 ** 64 + 1)]
        if kind == RATIONAL:
            bases += [Fraction(1, 2), Fraction(-2, 3)]
        bases += values[20:40]
        for a in bases:
            for n in exponents:
                if abs(n) >= POWER_LIMIT and a not in (0, 1, -1, 2, Fraction(1, 2)):
                    continue
                if abs(n) > 5 and abs(a) > 2 ** 64:
                    continue
                lines.append(f'op {kind} pow {_text(a)} {_text(n)} {_outcome(power, a, n, kind)}')
        for op in (ADD, MUL):
            for _ in range(60):
                count = rng.below(6)
                vs = [values[rng.below(len(values))] for _ in range(count)]
                text = ' '.join(_text(v) for v in vs)
                lines.append(f'fold {kind} {op.lower()} {count}{" " if vs else ""}{text} {_outcome(fold, op, vs, kind)}')
    return lines


def main(argv) -> int:
    for q in (Fraction(0), Fraction(2), Fraction(-3, 4), Fraction(256)):
        print(f'{str(q):6} {rational_bytes(q).hex(" ")}')
    self_test()
    text = '\n'.join(vector_lines()) + '\n'
    if argv == ['--write']:
        VECTOR_FILE.write_text(text, encoding='utf-8')
        print(f'wrote {VECTOR_FILE.name}')
        return 0
    if argv:
        print('usage: literal_reference.py [--write]', file=sys.stderr)
        return 2
    if not VECTOR_FILE.is_file() or VECTOR_FILE.read_text(encoding='utf-8') != text:
        print(f'{VECTOR_FILE.name} is missing or stale: run with --write', file=sys.stderr)
        return 1
    print(f'{VECTOR_FILE.name} is up to date')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
