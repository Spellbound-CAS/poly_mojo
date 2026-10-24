"""Exact scalar comparisons and sign-variation checks, independent of producers."""
from std.bit import count_trailing_zeros
from apn_mojo import Integer, Rational, lcm as integer_lcm
from ._coefficient_ops import _copy
from .coefficients import QQ, NoContext, _work
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _sum, _product
from .evaluate import _derivative_dense


def _compare(a: QQ, b: QQ, mut budget: Budget) raises PolyError -> Int:
    # Exact ordering without forming the normalized difference.
    budget.charge(1)
    try:
        return _order(a.value, b.value)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _order(a: Rational, b: Rational) raises -> Int:
    # sign(a - b) by cross-multiplication: two products, no gcd, and no
    # continued-fraction expansion of nearly equal operands.
    var sa = a.sign()
    var sb = b.sign()
    if sa != sb:
        return -1 if sa < sb else 1
    var da = a.denominator()
    var db = b.denominator()
    if da == db:
        var difference = a.numerator() - b.numerator()
        return difference.sign()
    var ea = _power_of_two(da)
    var eb = _power_of_two(db)
    if ea >= 0 and eb >= 0:
        # Dyadic operands compare after one shift instead of two products.
        if ea < eb:
            return ((a.numerator() << (eb - ea)) - b.numerator()).sign()
        return (a.numerator() - (b.numerator() << (ea - eb))).sign()
    return (a.numerator() * db - b.numerator() * da).sign()


def _twos(x: Integer) -> Int:
    """The number of trailing zero bits of a nonzero magnitude, from its words."""
    for i in range(x._word_count()):
        var word = x._word(i)
        if word != 0:
            return 32 * i + Int(count_trailing_zeros(word))
    return 0


def _power_of_two(x: Integer) -> Int:
    """k when |x| = 2^k, else -1."""
    var bits = x.magnitude_bit_length()
    if bits == 0:
        return -1
    var k = _twos(x)
    return k if k == bits - 1 else -1


def _midpoint(a: QQ, b: QQ, mut budget: Budget) raises PolyError -> QQ:
    var two = QQ.one(NoContext()).add(QQ.one(NoContext()), budget)
    return a.add(b, budget).mul(two.inverse(budget), budget)


def _homogeneous(p: DensePoly[QQ], x: QQ) raises -> Tuple[Integer, Integer]:
    # For x = a/b with b > 0 and the positive common denominator L of the
    # coefficients, returns (L b^n p(a/b), L b^n) as integers.
    var n = len(p.coefficients)
    var common = Integer(1)
    for i in range(n):
        ref value = p.coefficients[i].value
        if value.denominator() != 1:
            common = integer_lcm(common, value.denominator())
    var a = x.value.numerator()
    var b = x.value.denominator()
    var unit = b == 1
    var power = Integer(1)
    ref top = p.coefficients[n - 1].value
    var result = top.numerator() * (common // top.denominator())
    for i in range(n - 2, -1, -1):
        ref value = p.coefficients[i].value
        var term = value.numerator()
        if value.denominator() != common:
            term = term * (common // value.denominator())
        if not unit:
            power = power * b
            term = term * power
        result = result * a + term
    return (result^, common * power)


def _charge_evaluation(p: DensePoly[QQ], x: QQ, mut budget: Budget) raises PolyError:
    var n = len(p.coefficients)
    var coefficient_bytes = 0
    for i in range(n):
        budget.charge(1)
        coefficient_bytes = _sum(coefficient_bytes, p.coefficients[i].size())
    var x_bytes = x.size()
    for i in range(n - 1):
        budget.charge(1)
        _work(_sum(coefficient_bytes, _product(i, x_bytes)), x_bytes, True, budget)


def _value_at(p: DensePoly[QQ], x: QQ, mut budget: Budget) raises PolyError -> QQ:
    """The exact value p(x), normalized once instead of at every Horner step."""
    if len(p.coefficients) == 0:
        return QQ.zero(NoContext())
    _charge_evaluation(p, x, budget)
    try:
        var parts = _homogeneous(p, x)
        _work(parts[0].magnitude_bit_length() // 8 + 1, parts[1].magnitude_bit_length() // 8 + 1, True, budget)
        return QQ(Rational(parts[0], parts[1]))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _sign_at(p: DensePoly[QQ], x: QQ, mut budget: Budget) raises PolyError -> Int:
    """The exact sign of p(x), without rational normalization.

    For x = a/b with b > 0 and L the positive common denominator of the
    coefficients, L b^n p(a/b) is an integer of the same sign; homogeneous
    Horner computes it with integer products and sums only. The charge is
    a size bound fixed before the evaluation.
    """
    if len(p.coefficients) == 0:
        return 0
    _charge_evaluation(p, x, budget)
    try:
        return _homogeneous(p, x)[0].sign()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _variations(polynomials: List[DensePoly[QQ]], a: QQ, mut budget: Budget, side: Int = 0) raises PolyError -> Int:
    var previous = 0
    var count = 0
    for p in polynomials:
        budget.charge(1)
        var sign = _one_sided(p, _sign_at(p, a, budget), a, side, budget)
        if sign != 0:
            if previous != 0 and previous != sign:
                count += 1
            previous = sign
    return count


def _chain_variations(polynomials: List[DensePoly[QQ]], quotients: List[DensePoly[QQ]], a: QQ, mut budget: Budget, side: Int = 0) raises PolyError -> Int:
    # Equal to _variations for a chain whose recurrences
    # polynomials[i - 1] == quotients[i - 1] * polynomials[i] - polynomials[i + 1]
    # hold: callers either built the chain or have verified it. Two Horner
    # evaluations and one short quotient per entry give the same exact values.
    var n = len(polynomials)
    if n < 3 or len(quotients) + 1 < n:
        return _variations(polynomials, a, budget, side)
    var previous = 0
    var count = 0
    var before = QQ.zero(NoContext())
    var current = QQ.zero(NoContext())
    for i in range(n):
        budget.charge(1)
        if i < 2:
            before = current^
            current = _value_at(polynomials[i], a, budget)
        else:
            var following = _value_at(quotients[i - 2], a, budget).mul(current, budget).sub(before, budget)
            before = current^
            current = following^
        var sign = _one_sided(polynomials[i], current.value.sign(), a, side, budget)
        if sign != 0:
            if previous != 0 and previous != sign:
                count += 1
            previous = sign
    return count


def _one_sided(p: DensePoly[QQ], value_sign: Int, a: QQ, side: Int, mut budget: Budget) raises PolyError -> Int:
    var sign = value_sign
    if sign == 0 and side != 0 and not p.is_zero():
        # At an endpoint root, the first nonzero Taylor coefficient
        # gives the exact one-sided sign, including repeated roots.
        var derivative = _derivative_dense(p, budget)
        var parity = 1
        while True:
            budget.charge(1)
            sign = _sign_at(derivative, a, budget)
            if sign != 0:
                if side < 0 and parity:
                    sign = -sign
                break
            derivative = _derivative_dense(derivative, budget)
            parity = 1 - parity
    return sign


# Certified signs of integer polynomials at rationals, by fixed-point Horner
# with a proved error bound and an exact fallback (used by _real_refine and
# by short Sturm chains).

def _bits(x: Integer) -> Int:
    return x.magnitude_bit_length()


def _power_shift(q: Integer) raises -> Int:
    return _power_of_two(q)


def _error_bits(n: Int, p: Integer, q: Integer) -> Int:
    # Floor rounding in each of the n Horner steps errs by less than one unit,
    # and later steps scale it by |x|, so the total is below
    # sum_{j<n} |x|^j <= n max(1, |x|)^(n-1) units.
    var lx = max(0, _bits(p) - _bits(q) + 1)
    var log_n = 0
    while (1 << log_n) < n:
        log_n += 1
    return log_n + (n - 1) * lx + 1


def _fixed(c: List[Integer], p: Integer, q: Integer, shift: Int, scale: Int) raises -> Integer:
    # 2^scale * c(p/q) within the error of _error_bits, by floor-rounded Horner.
    # Narrow values at points of at most two words take the limb path,
    # which saves an allocation per step; wider ones gain from Integer's
    # multiplication.
    if shift >= 0 and p.magnitude_bit_length() <= 128 and _fixed_width(c, p, shift, scale) <= 40:
        return _fixed_limbs(c, p, shift, scale)
    var n = len(c) - 1
    var t = c[n] << scale
    for i in range(n - 1, -1, -1):
        var product = t * p
        if shift >= 0:
            t = product >> shift
        else:
            t = product // q
        t = t + (c[i] << scale)
    return t^


def _magnitude_limbs(x: Integer) -> List[UInt64]:
    # |x| as 64-bit words, lowest first.
    var count = x._word_count()
    var m = List[UInt64](length=(count + 1) // 2, fill=0)
    for j in range(count):
        m[j // 2] |= UInt64(x._word(j)) << UInt64(32 * (j % 2))
    return m^


struct _Words(Movable):
    """The magnitudes of integer coefficients as 64-bit words in one buffer."""
    var data: List[UInt64]
    var start: List[Int]
    """Coefficient i takes data[start[i]:start[i + 1]], lowest word first."""
    var negative: List[Bool]
    var top: Int
    """The largest bit length."""

    def __init__(out self, c: List[Integer]):
        self.start = List[Int](capacity=len(c) + 1)
        self.negative = List[Bool](capacity=len(c))
        self.top = 0
        var total = 0
        for i in range(len(c)):
            self.start.append(total)
            total += (c[i]._word_count() + 1) // 2
        self.start.append(total)
        self.data = List[UInt64](length=total, fill=0)
        for i in range(len(c)):
            ref x = c[i]
            var base = self.start[i]
            for j in range(x._word_count()):
                self.data[base + j // 2] |= UInt64(x._word(j)) << UInt64(32 * (j % 2))
            self.negative.append(x.sign() < 0)
            self.top = max(self.top, x.magnitude_bit_length())


def _add_shifted(p: Pointer[mut=True, UInt64, _], w: Int, words: _Words, index: Int, shift: Int):
    # t +- c_index 2^shift modulo 2^(64 w), on two's complement limbs.
    var subtract = words.negative[index]
    var first = words.start[index]
    var length = words.start[index + 1] - first
    var m = words.data.unsafe_ptr().unsafe_offset(first)
    var offset = shift // 64
    var bits = shift % 64
    var carry = UInt64(0)
    var k = offset
    var i = 0
    var previous = UInt64(0)
    while k < w and (i <= length or carry != 0):
        var current = m.unsafe_offset(i).unsafe_load() if i < length else UInt64(0)
        var part = (current << UInt64(bits)) | (previous >> UInt64(64 - bits) if bits else UInt64(0)) if i <= length else UInt64(0)
        previous = current
        var x = p.unsafe_offset(k).unsafe_load()
        if subtract:
            var d = x - part
            var b1 = UInt64(1) if x < part else UInt64(0)
            var e = d - carry
            carry = b1 | (UInt64(1) if d < carry else UInt64(0))
            p.unsafe_offset(k).unsafe_store(e)
        else:
            var sum = x + part
            var c1 = UInt64(1) if sum < x else UInt64(0)
            var e = sum + carry
            carry = c1 | (UInt64(1) if e < sum else UInt64(0))
            p.unsafe_offset(k).unsafe_store(e)
        k += 1
        i += 1


def _fixed_width(c: List[Integer], p: Integer, shift: Int, scale: Int) -> Int:
    # 64-bit words that bound |t| in _fixed: below 2^scale (n + 1) max|c|
    # max(1, |x|)^n plus the error.
    var n = len(c) - 1
    var top = 0
    for i in range(n + 1):
        top = max(top, c[i].magnitude_bit_length())
    var lx = max(0, p.magnitude_bit_length() - shift + 1)
    return (scale + top + n * lx + 2 * _bits(Integer(n + 1)) + 8) // 64 + 2


def _fixed_limbs(c: List[Integer], p: Integer, shift: Int, scale: Int) raises -> Integer:
    """_fixed for q = 2^shift, with the same floor roundings, on two's complement limbs.

    Each step multiplies by |p|, negates for a negative p, and shifts right
    arithmetically, which floors as Integer's >> does: the result is the
    same integer, without an allocation per step.
    """
    var n = len(c) - 1
    var words = _Words(c)
    var top = words.top
    var factor = _magnitude_limbs(p)
    var negative = p.sign() < 0
    var r = len(factor)
    var w = _fixed_width(c, p, shift, scale)
    var t = List[UInt64](length=w + r + 2, fill=0)
    var a = t.unsafe_ptr()
    var product = List[UInt64](length=w + r + 2, fill=0)
    var b = product.unsafe_ptr()
    # Only the first `used` limbs of t carry its value; the rest repeat its
    # sign and are written only when a step reads them.
    var used = (scale + top) // 64 + 2
    _add_shifted(a, used, words, n, scale)
    var whole = shift // 64
    var part = shift % 64
    for i in range(n - 1, -1, -1):
        # product = t |p| over used + r limbs, t sign-extended.
        var fill = UInt64.MAX if (a.unsafe_offset(used - 1).unsafe_load() >> 63) != 0 else UInt64(0)
        var span = used + r
        for k in range(used, span):
            a.unsafe_offset(k).unsafe_store(fill)
        for k in range(span):
            b.unsafe_offset(k).unsafe_store(0)
        for j in range(r):
            var f = UInt128(factor[j])
            if f == 0:
                continue
            var carry = UInt128(0)
            for k in range(span - j):
                var x = UInt128(a.unsafe_offset(k).unsafe_load()) * f + UInt128(b.unsafe_offset(j + k).unsafe_load()) + carry
                b.unsafe_offset(j + k).unsafe_store(UInt64(x & UInt128(UInt64.MAX)))
                carry = x >> 64
        if negative:
            var borrow = UInt64(1)
            for k in range(span):
                var x = ~b.unsafe_offset(k).unsafe_load() + borrow
                borrow = UInt64(1) if borrow == 1 and x == 0 else UInt64(0)
                b.unsafe_offset(k).unsafe_store(x)
        # t = floor(product / 2^shift): an arithmetic shift right.
        var sign_fill = UInt64.MAX if (b.unsafe_offset(span - 1).unsafe_load() >> 63) != 0 else UInt64(0)
        var kept = max(1, span - whole)
        for k in range(kept):
            var lo = b.unsafe_offset(k + whole).unsafe_load() if k + whole < span else sign_fill
            var hi = b.unsafe_offset(k + whole + 1).unsafe_load() if k + whole + 1 < span else sign_fill
            a.unsafe_offset(k).unsafe_store((lo >> UInt64(part)) | (hi << UInt64(64 - part)) if part else lo)
        # Room for the sum with the next coefficient, sign-extended.
        var bound = max(kept, (scale + c[i].magnitude_bit_length()) // 64 + 1) + 1
        for k in range(kept, bound):
            a.unsafe_offset(k).unsafe_store(sign_fill)
        _add_shifted(a, bound, words, i, scale)
        used = bound
        while used > 1:
            var high = a.unsafe_offset(used - 1).unsafe_load()
            var below = a.unsafe_offset(used - 2).unsafe_load()
            if high != (UInt64.MAX if (below >> 63) != 0 else UInt64(0)):
                break
            used -= 1
    w = used
    # Back to an Integer: magnitude words and sign.
    var below = (a.unsafe_offset(w - 1).unsafe_load() >> 63) != 0
    var out = List[UInt32](length=2 * w, fill=0)
    var borrow = UInt64(1)
    for k in range(w):
        var x = a.unsafe_offset(k).unsafe_load()
        if below:
            x = ~x + borrow
            borrow = UInt64(1) if borrow == 1 and x == 0 else UInt64(0)
        out[2 * k] = UInt32(x & 0xFFFFFFFF)
        out[2 * k + 1] = UInt32(x >> 32)
    return Integer._from_words(out^, below)


def _exact_sign(c: List[Integer], p: Integer, q: Integer) raises -> Int:
    # sign(q^n c(p/q)) by homogeneous Horner; q is positive.
    var n = len(c) - 1
    var t = c[n].copy()
    var power = Integer(1)
    for i in range(n - 1, -1, -1):
        power = power * q
        t = t * p + c[i] * power
    return t.sign()


def _sign(c: List[Integer], x: Rational) raises -> Int:
    var n = len(c) - 1
    if n <= 0:
        return c[0].sign() if n == 0 else 0
    var p = x.numerator()
    var q = x.denominator()
    if not p:
        return c[0].sign()
    var shift = _power_shift(q)
    var error = _error_bits(n, p, q)
    var exact_bits = n * (_bits(p) + _bits(q))
    # Near a root, |c(x)| is about |c'| times the grid step of x, so start
    # with as many guard bits as x has fractional bits.
    var extra = 64 + _bits(q)
    while error + extra < exact_bits:
        var t = _fixed(c, p, q, shift, error + extra)
        if t and _bits(t) - 1 > error + 1:
            return t.sign()
        extra *= 4
    return _exact_sign(c, p, q)
