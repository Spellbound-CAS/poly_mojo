"""Rational series by integer arithmetic: composition, inverses and logarithms.

Inverse and logarithm. With f = F / d and c = F_0, U(x) = F(c x) / c is an
integer series with constant term one, so G = 1 / U is integral and
(1 / f)_k = d G_k / c^(k+1). Newton steps for G run on integers through the
low and middle products of `_ntt`, and only the outputs are reduced. For
the logarithm f_0 = 1 forces c = d, and f'/f has coefficients
(F^' G)_k / c^(k+1) with the integer series F^'_j = F'_j c^j.

Composition.

With f = F / d_f and g = G / d_g for integer polynomials F and G,
f(g) = H / (d_f d_g^(m-1)) where H = sum_i F_i d_g^(m-1-i) G^i over the m
coefficients of f. H is computed with integer arithmetic only: the baby steps
G^0 ... G^k and the giant step G^k go through the transform products of
`_ntt`, the block sums are integer multiply-adds, and one division per
output coefficient reduces the fractions. No rational gcd runs inside the
loops.
"""
from apn_mojo import Integer, Rational, lcm
from .coefficient_buffer import _CoefficientView, _view
from .coefficients import ZZ, QQ
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _product, _sum
from ._ntt import _ntt_low_product, _ntt_middle_product, _magnitude, _integer


def _common(values: _CoefficientView[QQ, _]) raises -> Integer:
    """The lcm of the denominators."""
    var result = Integer(1)
    for i in range(len(values)):
        var d = values[i].value.denominator()
        if d != Integer(1):
            result = lcm(result, d)
    return result^


def _scaled(values: _CoefficientView[QQ, _], scale: Integer) raises -> List[ZZ]:
    var result = List[ZZ](capacity=len(values))
    for i in range(len(values)):
        ref v = values[i].value
        result.append(ZZ(v.numerator() * (scale // v.denominator())))
    return result^


def _powers(base: Integer, count: Int) raises -> List[Integer]:
    """base^0 ... base^(count - 1)."""
    var result = List[Integer](capacity=count)
    var power = Integer(1)
    for _ in range(count):
        result.append(power)
        power = power * base
    return result^


def _fractions(values: List[ZZ], denominator: Integer, length: Int) raises -> List[QQ]:
    var result = List[QQ](capacity=length)
    for t in range(length):
        result.append(QQ(Rational(values[t].value, denominator)) if t < len(values) else QQ(Rational(0)))
    return result^


def _low(a: List[ZZ], b: List[ZZ], length: Int, mut budget: Budget) raises PolyError -> Optional[List[ZZ]]:
    return _ntt_low_product(_view(a), _view(b), length, budget)


def _compose_rational(outer: _CoefficientView[QQ, _], inner: _CoefficientView[QQ, _], length: Int, mut budget: Budget) raises PolyError -> Optional[List[QQ]]:
    """outer(inner) modulo x^length, or None when a transform product does not apply."""
    var m = len(outer)
    if m == 0 or length <= 0:
        return List[QQ]()
    var d_f: Integer
    var d_g: Integer
    var f: List[ZZ]
    var g: List[ZZ]
    var scales: List[Integer]
    try:
        d_f = _common(outer)
        d_g = _common(inner)
        f = _scaled(outer, d_f)
        g = _scaled(inner, d_g)
        scales = _powers(d_g, m) if d_g != Integer(1) else List[Integer]()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    if len(scales):
        # F_i d_g^(m-1-i).
        for i in range(m):
            f[i] = f[i].mul(ZZ(scales[m - 1 - i]), budget)
    var k = 1
    while k * k < m:
        k += 1
    var blocks = (m + k - 1) // k
    var valuation = len(g) == 0 or g[0].is_zero()
    _allocation[List[ZZ]](k + 1, budget)
    var powers = List[List[ZZ]](capacity=k + 1)
    powers.append([ZZ(Integer(1))])
    for i in range(1, k + 1):
        var next = _low(powers[i - 1], g, length, budget)
        if not next:
            return None
        powers.append(next.take())
    # With word-sized F_i, each row entry is one run of word-by-limbs
    # multiply-adds over the baby steps' limbs instead of k integer products.
    var words = True
    for i in range(m):
        if f[i].value.magnitude_bit_length() > 63:
            words = False
            break
    var weights = List[Int64]()
    var limbs = List[List[List[UInt64]]]()
    var signs = List[List[Bool]]()
    var widest = 0
    if words:
        _allocation[Int64](m, budget)
        try:
            for i in range(m):
                weights.append(f[i].value.to_native_exact[DType.int64]())
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        for i in range(k):
            ref p = powers[i]
            var mags = List[List[UInt64]](capacity=len(p))
            var negs = List[Bool](capacity=len(p))
            for t in range(len(p)):
                var mag = List[UInt64]()
                _magnitude(p[t].value, mag)
                widest = max(widest, len(mag))
                mags.append(mag^)
                negs.append(p[t].value.sign() < 0)
            limbs.append(mags^)
            signs.append(negs^)
    var bytes = List[UInt8](length=8 * (widest + 2), fill=0)
    var result = List[ZZ]()
    for j in range(blocks - 1, -1, -1):
        var limit = max(0, length - j * k) if valuation else length
        _charge_size(_product(limit, k) + 1, budget)
        _allocation[ZZ](limit, budget)
        var row = List[ZZ](length=limit, fill=ZZ(Integer(0)))
        if words:
            _charge_size(_product(_product(limit, k), widest + 2), budget)
            var acc = List[UInt64](length=widest + 2, fill=0)
            for t in range(limit):
                for w in range(len(acc)):
                    acc[w] = 0
                var any = False
                for i in range(k):
                    var index = j * k + i
                    if index >= m or weights[index] == 0 or t >= len(limbs[i]) or len(limbs[i][t]) == 0:
                        continue
                    var x = weights[index]
                    var magnitude = UInt64(x) if x > 0 else UInt64(0) - UInt64(x)
                    if (x < 0) != signs[i][t]:
                        _submul(acc, limbs[i][t], magnitude)
                    else:
                        _addmul(acc, limbs[i][t], magnitude)
                    any = True
                if not any:
                    continue
                var below = (acc[len(acc) - 1] >> 63) == 1
                if below:
                    var carry = UInt64(1)
                    for w in range(len(acc)):
                        var v = ~acc[w] + carry
                        carry = UInt64(1) if (carry == 1 and v == 0) else UInt64(0)
                        acc[w] = v
                var used = len(acc)
                while used and acc[used - 1] == 0:
                    used -= 1
                if used:
                    row[t] = ZZ(_integer(acc, used, below, bytes))
        else:
            for i in range(k):
                var index = j * k + i
                if index >= m or f[index].is_zero():
                    continue
                ref p = powers[i]
                for t in range(min(limit, len(p))):
                    if not p[t].is_zero():
                        row[t] = row[t].add(f[index].mul(p[t], budget), budget)
        if len(result):
            var shifted = _low(result, powers[k], limit, budget)
            if not shifted:
                return None
            var low = shifted.take()
            for t in range(min(limit, len(low))):
                row[t] = row[t].add(low[t], budget)
        result = row^
    try:
        var denominator = d_f * (d_g ** (m - 1)) if len(scales) else d_f
        return _fractions(result, denominator, min(length, len(result)))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _middle(a: List[ZZ], b: List[ZZ], lo: Int, hi: Int, mut budget: Budget) raises PolyError -> Optional[List[ZZ]]:
    return _ntt_middle_product(_view(a), _view(b), lo, hi, budget)


comptime _INVERSE_QUADRATIC = 512
"""Orders below which inverses of word-sized unit series run the quadratic recurrence (measured)."""


def _inverse_words(u: List[ZZ], n: Int, mut budget: Budget) raises PolyError -> Optional[List[ZZ]]:
    """1 / u modulo x^n for an integer series with u_0 = 1, or None unless u_1 ... u_(n-1) are signed words.

    G_j = -sum_(i>=1) u_i G_(j-i): each coefficient is one run of
    word-by-limbs multiply-adds into a two's-complement accumulator, with no
    division at all.
    """
    var count = min(len(u), n)
    _allocation[Int64](count, budget)
    var weights = List[Int64](length=max(count, 1), fill=0)
    try:
        for i in range(1, count):
            if u[i].value.magnitude_bit_length() > 63:
                return None
            weights[i] = u[i].value.to_native_exact[DType.int64]()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var magnitudes = List[List[UInt64]](capacity=n)
    var negative = List[Bool](capacity=n)
    magnitudes.append([UInt64(1)])
    negative.append(False)
    var widest = 1
    for j in range(1, n):
        var terms = min(j, count - 1)
        var width = widest + 2
        _charge_size(_product(max(terms, 1), width), budget)
        var acc = List[UInt64](length=width, fill=0)
        for i in range(1, terms + 1):
            var w = weights[i]
            if w == 0 or len(magnitudes[j - i]) == 0:
                continue
            var m = UInt64(w) if w > 0 else UInt64(0) - UInt64(w)
            # Accumulate -sum u_i G_(j-i) directly.
            if (w < 0) != negative[j - i]:
                _addmul(acc, magnitudes[j - i], m)
            else:
                _submul(acc, magnitudes[j - i], m)
        var below = (acc[width - 1] >> 63) == 1
        if below:
            var carry = UInt64(1)
            for k in range(width):
                var v = ~acc[k] + carry
                carry = UInt64(1) if (carry == 1 and v == 0) else UInt64(0)
                acc[k] = v
        while len(acc) and acc[len(acc) - 1] == 0:
            _ = acc.pop()
        widest = max(widest, len(acc))
        magnitudes.append(acc^)
        negative.append(below)
    _allocation[ZZ](n, budget)
    var result = List[ZZ](capacity=n)
    var bytes = List[UInt8](length=8 * widest, fill=0)
    for j in range(n):
        ref limbs = magnitudes[j]
        result.append(ZZ(_integer(limbs, len(limbs), negative[j], bytes)) if len(limbs) else ZZ(Integer(0)))
    return result^


def _unit_inverse(u: List[ZZ], n: Int, mut budget: Budget) raises PolyError -> Optional[List[ZZ]]:
    """1 / u modulo x^n for an integer series with constant term one, by Newton steps (the word recurrence below _INVERSE_QUADRATIC)."""
    if n < _INVERSE_QUADRATIC:
        var words = _inverse_words(u, n, budget)
        if words:
            return words^
    var g: List[ZZ] = [ZZ(Integer(1))]
    var known = 1
    while known < n:
        budget.charge(1)
        var next = min(2 * known, n)
        var prefix = List[ZZ](capacity=min(len(u), next))
        for i in range(min(len(u), next)):
            prefix.append(u[i].copy())
        # u g = 1 + x^known e + O(x^next); the next coefficients are -g e.
        var e = _middle(prefix, g, known, next, budget)
        if not e:
            return None
        var t = _low(g, e.take(), next - known, budget)
        if not t:
            return None
        var correction = t.take()
        for i in range(next - known):
            g.append(correction[i].neg(budget) if i < len(correction) else ZZ(Integer(0)))
        known = next
    return g^


@fieldwise_init
struct _Integral(Movable):
    """f = F / d with c = F_0 and U(x) = F(c x) / c."""
    var f: List[ZZ]
    var d: Integer
    var c: Integer
    var powers: List[Integer]
    """c^0 ... c^n."""
    var u: List[ZZ]


def _integral(values: _CoefficientView[QQ, _], n: Int) raises -> _Integral:
    var d = _common(values)
    var f = _scaled(values, d)
    var c = f[0].value
    var powers = _powers(c, n + 2) if c != Integer(1) else List[Integer]()
    var u = List[ZZ](capacity=min(len(f), n))
    u.append(ZZ(Integer(1)))
    for k in range(1, min(len(f), n)):
        u.append(ZZ(f[k].value * powers[k - 1]) if len(powers) else f[k].copy())
    return _Integral(f^, d^, c^, powers^, u^)


def _ratio(value: Integer, scale: Integer, power: Integer) raises -> QQ:
    return QQ(Rational(value * scale, power))


def _inverse_rational(values: _CoefficientView[QQ, _], n: Int, mut budget: Budget) raises PolyError -> Optional[List[QQ]]:
    """1 / f modulo x^n, or None when f(0) = 0 or a product does not apply."""
    if n < 1 or len(values) == 0 or values[0].is_zero():
        return None
    _charge_size(n, budget)
    var data: _Integral
    try:
        data = _integral(values, n)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var g = _unit_inverse(data.u, n, budget)
    if not g:
        return None
    var inverse = g.take()
    _allocation[QQ](n, budget)
    var result = List[QQ](capacity=n)
    try:
        for k in range(n):
            var value = inverse[k].value if k < len(inverse) else Integer(0)
            result.append(_ratio(value, data.d, data.powers[k + 1] if len(data.powers) else Integer(1)))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return result^


def _log_rational(values: _CoefficientView[QQ, _], n: Int, mut budget: Budget) raises PolyError -> Optional[List[QQ]]:
    """log f modulo x^n for f(0) = 1, or None for other inputs."""
    if n < 1 or len(values) == 0 or values[0].value != Rational(1):
        return None
    _charge_size(n, budget)
    var data: _Integral
    try:
        data = _integral(values, n)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var result = List[QQ](capacity=n)
    result.append(QQ(Rational(0)))
    if n == 1:
        return result^
    var g = _unit_inverse(data.u, n - 1, budget)
    if not g:
        return None
    # F^'_j = (j + 1) F_(j+1) c^j.
    var derivative = List[ZZ](capacity=min(len(data.f), n) - 1)
    for j in range(min(len(data.f), n) - 1):
        var term = data.f[j + 1].mul(ZZ(Integer(j + 1)), budget)
        derivative.append(term.mul(ZZ(data.powers[j]), budget) if len(data.powers) else term^)
    var slope = List[ZZ]()
    if len(derivative):
        var product = _low(derivative, g.take(), n - 1, budget)
        if not product:
            return None
        slope = product.take()
    _allocation[QQ](n, budget)
    try:
        for k in range(1, n):
            var value = slope[k - 1].value if k - 1 < len(slope) else Integer(0)
            var power = data.powers[k] if len(data.powers) else Integer(1)
            result.append(QQ(Rational(value, power * Integer(k))))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return result^


comptime _EXP_BASE = 16
"""Coefficients of the exponential from the quadratic recurrence before Newton steps."""
comptime _EXP_DENOMINATOR_BITS = 64
"""Common denominators above this size go to the field routine: L_m carries d^(m - 1),
which outgrows the reduced coefficients when the input denominators vary widely."""


def _at(values: List[ZZ], i: Int) -> Integer:
    return values[i].value if i < len(values) else Integer(0)


def _at_zz(values: List[ZZ], i: Int) -> ZZ:
    return values[i].copy() if i < len(values) else ZZ(Integer(0))


def _exp_recurrence(f: List[ZZ], d: Integer, sign: Int, k: Int, scale: Integer) raises -> List[ZZ]:
    """Numerators of exp(sign f) modulo x^k over scale, from j N_j = sum_i i f_i N_(j-i)."""
    var result = List[ZZ](capacity=k)
    result.append(ZZ(scale))
    for j in range(1, k):
        var s = Integer(0)
        for i in range(1, min(j, len(f) - 1) + 1):
            if f[i].value.sign() != 0:
                s = s + f[i].value * Integer(i) * result[j - i].value
        var q = s // (d * Integer(j))
        result.append(ZZ(q if sign > 0 else -q))
    return result^


struct _ExactDivisor(Movable):
    """Exact quotients by one divisor, 2-adically: x / (2^e o) is
    ((x >> e) o^-1) modulo 2^bits in the symmetric range, for quotients below
    2^(bits - 1) in magnitude."""
    var shift: Int
    var mask: Integer
    var modulus: Integer
    var half: Integer
    var inverse: Integer

    def __init__(out self, divisor: Integer, bits: Int) raises:
        self.shift = (divisor & (-divisor)).magnitude_bit_length() - 1
        var odd = divisor >> self.shift
        var size = max(bits, 2)
        self.modulus = Integer(1) << size
        self.mask = self.modulus - Integer(1)
        self.half = Integer(1) << (size - 1)
        # Newton: y <- y (2 - o y) doubles the 2-adic precision of 1 / o.
        var y = Integer(1)
        var precision = 1
        while precision < size:
            precision = min(2 * precision, size)
            var window = (Integer(1) << precision) - Integer(1)
            y = (y * (Integer(2) - ((odd & window) * y))) & window
        self.inverse = y^

    def quotient(self, value: Integer) raises -> Integer:
        var q = (((value >> self.shift) & self.mask) * self.inverse) & self.mask
        return q - self.modulus if q >= self.half else q^


comptime _EXP_QUADRATIC = 2048
"""Orders below which exponentials with word-sized i F_i run the quadratic integer recurrence (measured)."""


def _addmul(mut acc: List[UInt64], h: List[UInt64], m: UInt64):
    """acc += h m in two's complement, carrying through acc's length."""
    var carry = UInt64(0)
    var pa = acc.unsafe_ptr()
    var ph = h.unsafe_ptr()
    for k in range(len(h)):
        var t = UInt128(ph.unsafe_offset(k).unsafe_load()) * UInt128(m) + UInt128(pa.unsafe_offset(k).unsafe_load()) + UInt128(carry)
        pa.unsafe_offset(k).unsafe_store(UInt64(t))
        carry = UInt64(t >> 64)
    var k = len(h)
    while carry and k < len(acc):
        var s = acc[k] + carry
        carry = UInt64(1) if s < carry else UInt64(0)
        acc[k] = s
        k += 1


def _submul(mut acc: List[UInt64], h: List[UInt64], m: UInt64):
    """acc -= h m in two's complement, borrowing through acc's length."""
    var borrow = UInt64(0)
    var pa = acc.unsafe_ptr()
    var ph = h.unsafe_ptr()
    for k in range(len(h)):
        var t = UInt128(ph.unsafe_offset(k).unsafe_load()) * UInt128(m) + UInt128(borrow)
        var low = UInt64(t)
        var x = pa.unsafe_offset(k).unsafe_load()
        pa.unsafe_offset(k).unsafe_store(x - low)
        borrow = UInt64(t >> 64) + (UInt64(1) if x < low else UInt64(0))
    var k = len(h)
    while borrow and k < len(acc):
        var x = acc[k]
        acc[k] = x - borrow
        borrow = UInt64(1) if x < borrow else UInt64(0)
        k += 1


def _exp_words(f: List[ZZ], d: Integer, n: Int, mut budget: Budget) raises PolyError -> Optional[List[ZZ]]:
    """Numerators of exp(f / d) modulo x^n over L = (n - 1)! d^(n - 1), or None when some i F_i leaves a signed word.

    H_j = L exp(f / d)_j is an integer with j d H_j = sum_i i F_i H_(j-i), so
    each H_j is one run of word-by-limbs multiply-adds and an exact division:
    quadratic in n, but without a transform or a rescaling per step.
    """
    var count = min(len(f), n)
    var bits = 0
    while (1 << bits) <= n:
        bits += 1
    _allocation[Int64](count, budget)
    var weights = List[Int64](length=max(count, 1), fill=0)
    var scale = Integer(1)
    # Exact divisions by j d run on words when every j d fits one.
    var word_divisor = d.magnitude_bit_length() + bits <= 63
    var dw = Int64(0)
    try:
        for i in range(1, count):
            if f[i].value.magnitude_bit_length() + bits > 62:
                return None
            weights[i] = f[i].value.to_native_exact[DType.int64]() * Int64(i)
        for j in range(1, n):
            scale = scale * Integer(j) * d
        if word_divisor:
            dw = d.to_native_exact[DType.int64]()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var magnitudes = List[List[UInt64]](capacity=n)
    var negative = List[Bool](capacity=n)
    var first = List[UInt64]()
    _magnitude(scale, first)
    var widest = len(first)
    magnitudes.append(first^)
    negative.append(False)
    var bytes = List[UInt8]()
    for j in range(1, n):
        var terms = min(j, count - 1)
        var width = widest + 2
        _charge_size(_product(max(terms, 1), width), budget)
        var acc = List[UInt64](length=width, fill=0)
        for i in range(1, terms + 1):
            var w = weights[i]
            if w == 0 or len(magnitudes[j - i]) == 0:
                continue
            var m = UInt64(w) if w > 0 else UInt64(0) - UInt64(w)
            if (w < 0) != negative[j - i]:
                _submul(acc, magnitudes[j - i], m)
            else:
                _addmul(acc, magnitudes[j - i], m)
        var below = (acc[width - 1] >> 63) == 1
        if below:
            var carry = UInt64(1)
            for k in range(width):
                var v = ~acc[k] + carry
                carry = UInt64(1) if (carry == 1 and v == 0) else UInt64(0)
                acc[k] = v
        var used = width
        while used and acc[used - 1] == 0:
            used -= 1
        acc.resize(used, 0)
        if used and word_divisor:
            var q = UInt64(j) * UInt64(dw if dw > 0 else -dw)
            var rem = UInt128(0)
            for k in range(used - 1, -1, -1):
                var cur = (rem << 64) | UInt128(acc[k])
                acc[k] = UInt64(cur // UInt128(q))
                rem = cur % UInt128(q)
            if dw < 0:
                below = not below
        elif used:
            if len(bytes) < 8 * used:
                bytes.resize(8 * used, 0)
            var exact: Integer
            try:
                exact = _integer(acc, used, below, bytes) // (Integer(j) * d)
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            below = exact.sign() < 0
            acc.clear()
            _magnitude(exact, acc)
        while len(acc) and acc[len(acc) - 1] == 0:
            _ = acc.pop()
        widest = max(widest, len(acc))
        magnitudes.append(acc^)
        negative.append(below)
    _allocation[ZZ](n, budget)
    var result = List[ZZ](capacity=n)
    for j in range(n):
        ref limbs = magnitudes[j]
        if len(limbs) == 0:
            result.append(ZZ(Integer(0)))
            continue
        if len(bytes) < 8 * len(limbs):
            bytes.resize(8 * len(limbs), 0)
        result.append(ZZ(_integer(limbs, len(limbs), negative[j], bytes)))
    return result^


def _exp_rational(values: _CoefficientView[QQ, _], n: Int, mut budget: Budget) raises PolyError -> Optional[List[QQ]]:
    """exp f modulo x^n for f(0) = 0, or None for other inputs.

    With f = F / d, k! d^k exp(f)_k is an integer, so exp(f) and exp(-f) modulo
    x^m are kept as numerators G and H over L_m = (m - 1)! d^(m - 1). A Newton
    step from k to m = 2k forms g h = 1 + x^k e, takes log g up to x^m from
    g' h (1 - x^k e), where the low part of g' h is f', and sets
    g <- g (1 + x^k t) and h <- h - x^k h (e + t) for t = f - log g. Every new
    numerator is an exact quotient by L_k^2; no gcd runs before the output.
    """
    if n < 2 or len(values) == 0 or not values[0].is_zero():
        return None
    _charge_size(n, budget)
    var count = min(len(values), n)
    var d: Integer
    var f: List[ZZ]
    var k = min(n, _EXP_BASE)
    var lam: Integer
    var g: List[ZZ]
    var h: List[ZZ]
    try:
        var prefix = List[QQ](capacity=count)
        for i in range(count):
            prefix.append(values[i].copy())
        d = _common(_view(prefix))
        if d.magnitude_bit_length() > _EXP_DENOMINATOR_BITS:
            return None
        f = _scaled(_view(prefix), d)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var quadratic = _exp_words(f, d, n, budget) if n < _EXP_QUADRATIC else Optional[List[ZZ]](None)
    try:
        if quadratic:
            # Numerators over L_n already: no Newton steps.
            g = quadratic.take()
            h = List[ZZ]()
            lam = Integer(1)
            k = n
        else:
            lam = Integer(1)
            for j in range(1, k):
                lam = lam * Integer(j) * d
            _charge_size(_product(k, k), budget)
            g = _exp_recurrence(f, d, 1, k, lam)
            h = _exp_recurrence(f, d, -1, k, lam)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    while k < n:
        budget.charge(1)
        var m = min(2 * k, n)
        var w = m - k
        # g h = 1 + x^k E / L^2.
        var em = _middle(g, h, k, m, budget)
        if not em:
            return None
        var e = em.take()
        # g' h over L^2 in degrees k - 1 ... m - 2.
        _allocation[ZZ](k, budget)
        var slope = List[ZZ](capacity=k - 1)
        for j in range(k - 1):
            slope.append(g[j + 1].mul(ZZ(Integer(j + 1)), budget))
        var pm = _middle(slope, h, k - 1, m - 1, budget)
        if not pm:
            return None
        var p = pm.take()
        # f' E in degrees 0 ... w - 2, with f' = F' / d.
        var fe = List[ZZ]()
        if w > 1:
            var fp = List[ZZ](capacity=w - 1)
            for a in range(w - 1):
                fp.append(_at_zz(f, a + 1).mul(ZZ(Integer(a + 1)), budget))
            var low = _low(fp, e, w - 1, budget)
            if not low:
                return None
            fe = low.take()
        _charge_size(_product(w, k), budget)
        var t: List[ZZ]
        var s: List[ZZ]
        var square: Integer
        var widen: Integer
        var pi: Integer
        try:
            square = lam * lam
            pi = Integer(1)
            for j in range(k, m):
                pi = pi * Integer(j)
            widen = d ** (w - 1)
            # t_j = R_j / (d L^2 j), R_j = F_j L^2 j - d P_(j-1) + (F' E)_(j-1-k);
            # over the common d L^2 pi the numerators are R_j pi / j.
            t = List[ZZ](capacity=w)
            s = List[ZZ](capacity=w)
            for i in range(w):
                var j = k + i
                var r = _at(f, j) * square * Integer(j) - d * _at(p, i)
                if i >= 1:
                    r = r + _at(fe, i - 1)
                var scaled = r * (pi // Integer(j))
                s.append(ZZ(_at(e, i) * d * pi + scaled))
                t.append(ZZ(scaled^))
        except err:
            raise PolyError(PolyErrorKind.CAPACITY, String(err))
        # g_(k+i) = (G T)_i d^(w-1) / L^2 and h_(k+i) = -(H S)_i d^(w-1) / L^2 over L_m.
        var gt = _low(g, t, w, budget)
        if not gt:
            return None
        var hs = _low(h, s, w, budget)
        if not hs:
            return None
        var gnew = gt.take()
        var hnew = hs.take()
        _charge_size(_product(_sum(k, w), 2), budget)
        try:
            var grow = pi * (d ** w)
            for j in range(k):
                g[j] = ZZ(g[j].value * grow)
                h[j] = ZZ(h[j].value * grow)
            # Every quotient by L^2 is exact: take it 2-adically.
            var top = 0
            for i in range(w):
                top = max(top, _at(gnew, i).magnitude_bit_length(), _at(hnew, i).magnitude_bit_length())
            var divisor = _ExactDivisor(square, top + widen.magnitude_bit_length() - square.magnitude_bit_length() + 3)
            for i in range(w):
                g.append(ZZ(divisor.quotient(_at(gnew, i) * widen)))
                h.append(ZZ(-divisor.quotient(_at(hnew, i) * widen)))
            lam = lam * grow
        except err:
            raise PolyError(PolyErrorKind.CAPACITY, String(err))
        k = m
    _allocation[QQ](n, budget)
    var result = List[QQ](capacity=n)
    try:
        # exp(f)_j = G_j / L_n with j! d^j exp(f)_j integral: divide by
        # L_n / (j! d^j) first so the gcd runs against j! d^j only.
        var cofactors = List[Integer](length=n, fill=Integer(1))
        for j in range(n - 2, -1, -1):
            cofactors[j] = cofactors[j + 1] * Integer(j + 1) * d
        var denominator = Integer(1)
        for j in range(n):
            if j:
                denominator = denominator * Integer(j) * d
            result.append(QQ(Rational(g[j].value // cofactors[j], denominator)))
    except err:
        raise PolyError(PolyErrorKind.CAPACITY, String(err))
    return result^


comptime _SQRT_DENOMINATOR_BITS = 64
"""Common denominators above this size keep the field routine."""


def _sqrt_rational(values: _CoefficientView[QQ, _], n: Int, mut budget: Budget) raises PolyError -> Optional[List[QQ]]:
    """sqrt(f) modulo x^n for f(0) = 1, or None for other inputs.

    With f = F / d, A(x) = f(4 d x) has integer coefficients F_j (4d)^j / d
    and constant one, and so do sqrt(A) and 1 / sqrt(A), since 4^i binom(1/2, i)
    and 4^i binom(-1/2, i) are integers. G = 1 / sqrt(A) doubles its precision
    by G <- G + x^known G E / 2 with A G^2 = 1 - x^known E, the halving exact,
    and sqrt(f)_j = (A G)_j / (4 d)^j.
    """
    if n < 1 or len(values) == 0 or values[0].value != Rational(1):
        return None
    _charge_size(n, budget)
    var count = min(len(values), n)
    var d: Integer
    var a: List[ZZ]
    try:
        var prefix = List[QQ](capacity=count)
        for i in range(count):
            prefix.append(values[i].copy())
        d = _common(_view(prefix))
        if d.magnitude_bit_length() > _SQRT_DENOMINATOR_BITS:
            return None
        var f = _scaled(_view(prefix), d)
        var base = Integer(4) * d
        var power = Integer(1)
        a = List[ZZ](capacity=count)
        for j in range(count):
            a.append(ZZ(f[j].value * power // d))
            power = power * base
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var g: List[ZZ] = [ZZ(Integer(1))]
    var known = 1
    while known < n:
        budget.charge(1)
        var next = min(2 * known, n)
        var width = next - known
        var sq = _low(g, g, next, budget)
        if not sq:
            return None
        var prefix = List[ZZ](capacity=min(len(a), next))
        for i in range(min(len(a), next)):
            prefix.append(a[i].copy())
        var em = _middle(prefix, sq.take(), known, next, budget)
        if not em:
            return None
        var tm = _low(g, em.take(), width, budget)
        if not tm:
            return None
        var t = tm.take()
        _charge_size(width, budget)
        try:
            for i in range(width):
                g.append(ZZ(-(_at(t, i) // Integer(2))))
        except err:
            raise PolyError(PolyErrorKind.CAPACITY, String(err))
        known = next
    var rm = _low(a, g, n, budget)
    if not rm:
        return None
    var root = rm.take()
    _allocation[QQ](n, budget)
    var result = List[QQ](capacity=n)
    try:
        var base = Integer(4) * d
        var power = Integer(1)
        for j in range(n):
            result.append(QQ(Rational(_at(root, j), power)))
            power = power * base
    except err:
        raise PolyError(PolyErrorKind.CAPACITY, String(err))
    return result^
