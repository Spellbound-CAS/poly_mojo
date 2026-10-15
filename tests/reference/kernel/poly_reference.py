"""Reference implementation of the kernel's FORMAL_POLYNOMIAL literal kind.

It is deliberately simple and independent of the kernel and of the
`polynomials` package: Python integers and fractions, and a polynomial as a
dict from monomials to nonzero Fractions, a monomial being a tuple of
(index, exponent) pairs with indices strictly increasing and exponents at
least 1. The kernel's polynomial operations, and the `polynomials`
package's conversion to and from the kernel's encoding, are tested against
it.

Each operation returns the result, or None when the kernel must decline
Run it from the repository root to print the encodings
of the specification's table:

    python3 tests/reference/kernel/poly_reference.py
"""
import struct
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from literal_reference import MAX_LITERAL_BYTES, POWER_LIMIT, CapacityError, rational_bytes, rational_from_bytes  # noqa: E402

INDEX_LIMIT = 2 ** 63  # INDET evaluates on integers i with |i| < 2**63


def monomial_key(monomial):
    """Monomials compare pair by pair, a proper prefix first: Python's tuple order."""
    return monomial


def canonical(poly):
    """The canonical term list: monomials strictly increasing, no zero coefficient."""
    return [(m, c) for m, c in sorted(poly.items(), key=lambda item: monomial_key(item[0])) if c != 0]


def encode(poly) -> bytes:
    """The kernel encoding; the values 0 and 1 have none as atoms."""
    terms = canonical(poly)
    if not terms or terms == [((), Fraction(1))]:
        raise ValueError('0 and 1 are the bound ZERO and ONE, not literal atoms')
    out = bytearray(struct.pack('<I', len(terms)))
    for monomial, coefficient in terms:
        out += struct.pack('<I', len(monomial))
        for index, exponent in monomial:
            out += struct.pack('<qI', index, exponent)
        encoded = rational_bytes(Fraction(coefficient))
        out += struct.pack('<I', len(encoded)) + encoded
    return bytes(out)


def decode(data: bytes):
    """Decode a canonical encoding; raise ValueError on anything else (ENCODING)."""
    position = 0

    def take(fmt):
        nonlocal position
        size = struct.calcsize(fmt)
        if position + size > len(data):
            raise ValueError('truncated')
        values = struct.unpack_from(fmt, data, position)
        position += size
        return values

    (count,) = take('<I')
    if count == 0:
        raise ValueError('no terms')
    poly, previous = {}, None
    for _ in range(count):
        (pairs,) = take('<I')
        monomial, last = [], None
        for _ in range(pairs):
            index, exponent = take('<qI')
            if exponent == 0 or (last is not None and index <= last):
                raise ValueError('not a canonical monomial')
            monomial.append((index, exponent))
            last = index
        monomial = tuple(monomial)
        (length,) = take('<I')
        if position + length > len(data):
            raise ValueError('truncated coefficient')
        coefficient = rational_from_bytes(data[position:position + length])
        position += length
        if coefficient == 0 or (previous is not None and monomial <= previous):
            raise ValueError('not canonical')
        poly[monomial] = coefficient
        previous = monomial
    if position != len(data):
        raise ValueError('trailing bytes')
    if list(poly.items()) == [((), Fraction(1))]:
        raise ValueError('1 is the bound ONE')
    return poly


def _clean(poly):
    return {m: c for m, c in poly.items() if c != 0}


def _checked(poly):
    if _clean(poly) and _clean(poly) != {(): Fraction(1)} and len(encode(_clean(poly))) > MAX_LITERAL_BYTES:
        raise CapacityError()
    return _clean(poly)


def add(p, q):
    out = dict(p)
    for m, c in q.items():
        out[m] = out.get(m, Fraction(0)) + c
    return _checked(out)


def neg(p):
    return {m: -c for m, c in p.items()}


def _mul_monomials(a, b):
    exponents = dict(a)
    for index, exponent in b:
        exponents[index] = exponents.get(index, 0) + exponent
    return tuple(sorted(exponents.items()))


def mul(p, q):
    out = {}
    for ma, ca in p.items():
        for mb, cb in q.items():
            m = _mul_monomials(ma, mb)
            out[m] = out.get(m, Fraction(0)) + ca * cb
    return _checked(out)


def power(p, n):
    """p**n for a natural number n <= POWER_LIMIT; None (decline) otherwise."""
    n = Fraction(n)
    if n.denominator != 1 or n < 0 or n > POWER_LIMIT:
        return None
    result, base, e = {(): Fraction(1)}, dict(p), int(n)
    while e:
        if e & 1:
            result = mul(result, base)
        e >>= 1
        if e:
            base = mul(base, base)
    return result


def indet(i):
    """X_i for an integer i with |i| < 2**63; None (decline) otherwise."""
    i = Fraction(i)
    if i.denominator != 1 or abs(i) >= INDEX_LIMIT:
        return None
    return {((int(i), 1),): Fraction(1)}


def coeff(c):
    """The constant polynomial c, for every rational c."""
    c = Fraction(c)
    return {(): c} if c != 0 else {}


def eq(p, q):
    """True for the same literal; None (decline) otherwise: EQ never evaluates to False."""
    return True if _clean(p) == _clean(q) else None


def table():
    """The rows of the kernel's encoding table."""
    return [
        ('2', coeff(2)),
        ('X_0', indet(0)),
        ('X_{-1}^2 - 1/2', add(mul(indet(-1), indet(-1)), coeff(Fraction(-1, 2)))),
    ]


def self_test():
    for _, poly in table():
        assert decode(encode(poly)) == poly
    x, y = indet(0), indet(1)
    left = power(add(x, y), 2)
    right = add(add(mul(x, x), mul(coeff(2), mul(x, y))), mul(y, y))
    assert left == right
    assert add(x, neg(x)) == {}
    assert eq(x, y) is None and eq(x, indet(0)) is True
    assert power(x, -1) is None and power(x, Fraction(1, 2)) is None and indet(Fraction(1, 2)) is None
    for bad in (b'', encode(x) + b'\x00', struct.pack('<I', 0)):
        try:
            decode(bad)
        except ValueError:
            continue
        raise AssertionError('accepted a bad encoding')


if __name__ == '__main__':
    self_test()
    for name, poly in table():
        print(f'{name:18} {encode(poly).hex(" ")}')
