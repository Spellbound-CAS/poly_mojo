"""Outward-rounded dyadic intervals for certified sign filters.

An interval is [low * 2^exponent, high * 2^exponent] with Integer mantissas
and one shared exponent, as in libpoly's dyadic intervals. Every operation
rounds outward, so the result contains every value of the exact operation on
the operand intervals: a sign decided here is the sign of the exact value.
Precision only bounds mantissa growth; it never affects soundness.
"""
from apn_mojo import Integer
from .errors import Budget, PolyError, PolyErrorKind, _charge_size
from .coefficients import QQ
from .dense import DensePoly
from .algebraic_types import RealAlgebraic


@fieldwise_init
struct _Dyadic(Copyable, Movable):
    var low: Integer
    var high: Integer
    var exponent: Int

    def sign(self) -> Int:
        """1 or -1 when every point has that sign, 0 when undecided."""
        if self.low.sign() > 0:
            return 1
        if self.high.sign() < 0:
            return -1
        return 0


def _bits(x: Integer) -> Int:
    return x.magnitude_bit_length()


def _floor_quotient(n: Integer, d: Integer) raises PolyError -> Integer:
    """floor(n / d) for d > 0, independent of the division's rounding convention."""
    try:
        var q = n // d
        while q * d > n:
            q = q - Integer(1)
        while (q + Integer(1)) * d <= n:
            q = q + Integer(1)
        return q^
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _ceil_shift(x: Integer, s: Int) raises PolyError -> Integer:
    """ceil(x / 2^s)."""
    try:
        return -((-x) >> s)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _floor_shift(x: Integer, s: Int) raises PolyError -> Integer:
    try:
        return x >> s
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _left(x: Integer, s: Int) raises PolyError -> Integer:
    if s == 0 or not x:
        return x.copy()
    try:
        return x << s
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _rounded(var x: _Dyadic, precision: Int, mut budget: Budget) raises PolyError -> _Dyadic:
    """Drop low mantissa bits outward until both mantissas fit the precision."""
    var bits = max(_bits(x.low), _bits(x.high))
    _charge_size(bits // 64 + 1, budget)
    if bits <= precision:
        return x^
    var s = bits - precision
    return _Dyadic(_floor_shift(x.low, s), _ceil_shift(x.high, s), x.exponent + s)


def _bounds(low: QQ, high: QQ, precision: Int, mut budget: Budget) raises PolyError -> _Dyadic:
    """A dyadic interval containing the rational interval [low, high]."""
    var a = low.value.numerator()
    var b = low.value.denominator()
    var c = high.value.numerator()
    var d = high.value.denominator()
    _charge_size((_bits(a) + _bits(b) + _bits(c) + _bits(d)) // 64 + 1, budget)
    # The exponent puts about `precision` bits in the larger mantissa.
    var top = max(max(_bits(a) - _bits(b), _bits(c) - _bits(d)), 0)
    var exponent = top - precision
    var lo: Integer
    var hi: Integer
    if exponent >= 0:
        lo = _floor_quotient(a, _left(b, exponent))
        hi = -_floor_quotient(-c, _left(d, exponent))
    else:
        lo = _floor_quotient(_left(a, -exponent), b)
        hi = -_floor_quotient(-_left(c, -exponent), d)
    return _Dyadic(lo^, hi^, exponent)


def _point(q: QQ, precision: Int, mut budget: Budget) raises PolyError -> _Dyadic:
    return _bounds(q, q, precision, budget)


def _of_real(x: RealAlgebraic, precision: Int, mut budget: Budget) raises PolyError -> _Dyadic:
    return _bounds(x.low, x.high, precision, budget)


def _aligned(a: _Dyadic, b: _Dyadic) raises PolyError -> Tuple[Integer, Integer, Integer, Integer, Int]:
    var e = min(a.exponent, b.exponent)
    return (_left(a.low, a.exponent - e), _left(a.high, a.exponent - e), _left(b.low, b.exponent - e), _left(b.high, b.exponent - e), e)


def _add(a: _Dyadic, b: _Dyadic, precision: Int, mut budget: Budget) raises PolyError -> _Dyadic:
    if not a.low and not a.high:
        return b.copy()
    if not b.low and not b.high:
        return a.copy()
    var parts = _aligned(a, b)
    _charge_size((_bits(parts[0]) + _bits(parts[2])) // 64 + 1, budget)
    var result: _Dyadic
    try:
        result = _Dyadic(parts[0] + parts[2], parts[1] + parts[3], parts[4])
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    # Rounding can exhaust the caller's work allowance. Keep that typed
    # failure outside the handler for arbitrary-precision allocation errors.
    return _rounded(result^, precision, budget)


def _mul(a: _Dyadic, b: _Dyadic, precision: Int, mut budget: Budget) raises PolyError -> _Dyadic:
    _charge_size((_bits(a.low) + _bits(a.high) + _bits(b.low) + _bits(b.high)) // 32 + 1, budget)
    var result: _Dyadic
    try:
        var p1 = a.low * b.low
        var p2 = a.low * b.high
        var p3 = a.high * b.low
        var p4 = a.high * b.high
        var low = p1.copy()
        var high = p1.copy()
        for p in [p2^, p3^, p4^]:
            if p < low:
                low = p.copy()
            if p > high:
                high = p.copy()
        result = _Dyadic(low^, high^, a.exponent + b.exponent)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return _rounded(result^, precision, budget)


def _zero() -> _Dyadic:
    return _Dyadic(Integer(0), Integer(0), 0)


def _dense_value(p: DensePoly[QQ], x: _Dyadic, precision: Int, mut budget: Budget) raises PolyError -> _Dyadic:
    """An enclosure of p over the interval x, by Horner's rule."""
    var value = _zero()
    for i in range(len(p.coefficients) - 1, -1, -1):
        budget.charge(1)
        value = _add(_mul(value, x, precision, budget), _point(p.coefficients[i], precision, budget), precision, budget)
    return value^


def _precision(width_low: QQ, width_high: QQ) -> Int:
    """Enough mantissa bits to resolve an interval of the given rational endpoints."""
    var w = width_high.value.denominator().magnitude_bit_length() + width_low.value.denominator().magnitude_bit_length()
    return max(64, w + 32)
