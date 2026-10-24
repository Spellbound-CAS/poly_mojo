"""Division of dense polynomials over odd prime fields on word residues.

Polynomials are `List[UInt64]` Montgomery residues as in `_nmod`, lowest
degree first. Short operands divide by dot products; longer ones take the
quotient from the transform division of `_ntt` (a Newton inverse of the
reversed divisor) and the remainder from the coefficients of q * b below the
divisor's degree, by the kernel dense multiplication would use (schoolbook,
Karatsuba, or transforms by `_uses_transform`). Every step charges the shared
budget in proportion to its coefficient work.
"""
from .errors import Budget, PolyError, _allocation, _charge_size, _product as _times, _sum
from ._nmod import _Field, _mulm, _mul, _schoolbook, _trim, _slice, _uses_transform, _subm, _from, _dot, _invm, _reversed
from ._ntt import _modular_product, _quotient_mod, _PRIMES, _PRIME_BITS, _reconstruct_zz, _bit_length
from ._nmod import _field
from .coefficients import CoefficientRing, ZZ
from .coefficient_buffer import _CoefficientView

comptime _DOT = 512
"""Divisor degree and quotient length from which the quotient uses the transform division."""
comptime _DOT_SMALL = 2048
"""The same for primes below 2^28, whose dot products run on SIMD lanes."""
comptime _SMALL_ZZ = 2048
"""Divisor length below which exact ZZ division runs on word primes without a product check."""
comptime _ZZ_PRIMES = 8
"""The most word primes the short ZZ division tries before leaving the case to transforms."""


def _product(a: List[UInt64], b: List[UInt64], limit: Int, f: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    """The coefficients of a * b below `limit`, trimmed, by the kernel dense multiplication would use."""
    if limit <= 0:
        return List[UInt64]()
    var x = _slice(a, 0, limit, budget)
    var y = _slice(b, 0, limit, budget)
    if len(x) == 0 or len(y) == 0:
        return List[UInt64]()
    var shorter = min(len(x), len(y))
    if shorter < 32:
        return _schoolbook(x, y, limit, f, budget)
    if _uses_transform(shorter, f.p):
        # A transform product of Montgomery residues aR and bR is abR^2; one
        # reduction returns it to Montgomery form.
        var c = _modular_product(x, y, min(limit, len(x) + len(y) - 1), f.p, False, budget)
        for i in range(len(c)):
            c[i] = _mulm(c[i], 1, f)
        _trim(c)
        return c^
    return _slice(_mul(x, y, f, budget), 0, limit, budget)


def _dot_divmod(a: List[UInt64], b: List[UInt64], f: _Field, mut budget: Budget) raises PolyError -> Tuple[List[UInt64], List[UInt64]]:
    """Long division by dot products, for len(a) >= len(b) >= 1.

    With B = rev(b), rev(q)_k = (rev(a)_k - sum_(j=1..k) B_j rev(q)_(k-j)) / B_0.
    Writing rev(q)_k at position m - 1 - k fills q in natural order and makes
    each sum one contiguous dot product; the remainder below the divisor's
    degree is a_k - sum_i q_i b_(k-i), again a dot product against B.
    """
    var n = len(b) - 1
    var m = len(a) - n
    _charge_size(_sum(_times(m, n), _times(n, min(m, n))) + 1, budget)
    _allocation[UInt64](_sum(m, n), budget)
    var rb = _reversed(b, n + 1, budget)
    var inverse = _invm(rb[0], f, budget)
    var q = List[UInt64](length=m, fill=0)
    var pb = rb.unsafe_ptr()
    for k in range(m):
        var t = min(k, n)
        var s = _dot(pb.unsafe_offset(1), q.unsafe_ptr().unsafe_offset(m - k), t, f) if t else UInt64(0)
        q[m - 1 - k] = _mulm(_subm(a[len(a) - 1 - k], s, f), inverse, f)
    var r = List[UInt64](capacity=n)
    var pq = q.unsafe_ptr()
    for k in range(n):
        var count = min(k, m - 1) + 1
        r.append(_subm(a[k], _dot(pq, pb.unsafe_offset(n - k), count, f), f))
    _trim(r)
    return (q^, r^)


def _field_divmod(a: List[UInt64], b: List[UInt64], f: _Field, mut budget: Budget) raises PolyError -> Tuple[List[UInt64], List[UInt64]]:
    """(quotient, remainder) of trimmed a by trimmed nonzero b, both trimmed."""
    var n = len(b) - 1
    if len(a) < len(b):
        return (List[UInt64](), a.copy())
    var m = len(a) - n
    var cutoff = _DOT_SMALL if f.p < UInt64(1) << 28 else _DOT
    if n < cutoff or m < cutoff:
        return _dot_divmod(a, b, f, budget)
    _allocation[UInt64](len(a) + len(b), budget)
    var x = List[UInt64](capacity=len(a))
    for v in a:
        x.append(_from(v, f))
    var y = List[UInt64](capacity=len(b))
    for v in b:
        y.append(_from(v, f))
    var plain = _quotient_mod(x, y, f.p, budget)
    var q = List[UInt64](capacity=len(plain))
    for v in plain:
        q.append(_mulm(v, f.r2, f))
    # q * b agrees with a from degree n up; the remainder is a - q * b below n.
    var low = _product(q, b, n, f, budget)
    var r = List[UInt64](capacity=n)
    for i in range(n):
        r.append(_subm(a[i], low[i] if i < len(low) else UInt64(0), f))
    _trim(r)
    return (q^, r^)


def _words_of[C: CoefficientRing](values: _CoefficientView[C, _]) -> Optional[List[Int64]]:
    """The coefficients as words, or None when one exceeds 63 bits."""
    var out = List[Int64](capacity=len(values))
    for i in range(len(values)):
        ref v = rebind[ZZ](values[i]).value
        if v.magnitude_bit_length() > 63:
            return None
        try:
            out.append(v.to_native_exact[DType.int64]())
        except:
            return None
    return out^


def _signed_residues(words: List[Int64], f: _Field) -> List[UInt64]:
    """Montgomery residues of signed words modulo f.p."""
    var out = List[UInt64](capacity=len(words))
    for w in words:
        var magnitude = UInt64(w) if w >= 0 else UInt64(0) - UInt64(w)
        var r = magnitude % f.p
        if w < 0 and r:
            r = f.p - r
        out.append(_mulm(r, f.r2, f))
    return out^


def _bits(words: List[Int64]) -> Int:
    var top = UInt64(0)
    for w in words:
        top = max(top, UInt64(w) if w >= 0 else UInt64(0) - UInt64(w))
    var bits = 0
    while bits < 64 and (UInt64(1) << UInt64(bits)) <= top:
        bits += 1
    return bits


def _zz_exact[C: CoefficientRing](a: _CoefficientView[C, _], b: _CoefficientView[C, _], mut budget: Budget) raises PolyError -> Optional[Optional[List[C]]]:
    """The exact ZZ quotient a / b (inner None when b does not divide a) for short divisors with word coefficients, or None.

    Over word primes p not dividing lc(b), b | a forces a zero remainder
    modulo p. When every remainder vanishes, the symmetric CRT candidate q
    satisfies a - q b = 0 modulo P, the product of the primes, and
    |a - q b| <= |a| + L |q| |b| < P / 2 proves a = q b without a product.
    Requires len(a) >= len(b) >= 2 with nonzero leading coefficients.
    """
    if len(b) >= _SMALL_ZZ:
        return None
    var x = _words_of(a)
    if not x:
        return None
    var y = _words_of(b)
    if not y:
        return None
    var bits_a = _bits(x.value())
    var bits_b = _bits(y.value())
    var m = len(a) - len(b) + 1
    var terms = _bit_length(min(m, len(b)))
    var primes = materialize[_PRIMES]()
    var residues = List[List[UInt64]]()
    var k = 1
    while k <= _ZZ_PRIMES:
        for i in range(len(residues), k):
            budget.charge(1)
            var f = _field(primes[i])
            var u = _signed_residues(x.value(), f)
            var v = _signed_residues(y.value(), f)
            if v[len(v) - 1] == 0:
                return None
            var qr = _field_divmod(u, v, f, budget)
            if len(qr[1]):
                return Optional(Optional[List[C]](None))
            var plain = List[UInt64](length=m, fill=0)
            for j in range(len(qr[0])):
                plain[j] = _from(qr[0][j], f)
            residues.append(plain^)
        var candidate = _reconstruct_zz[C](residues, k, m, budget)
        var bits_q = 0
        for i in range(len(candidate)):
            bits_q = max(bits_q, rebind[ZZ](candidate[i]).value.magnitude_bit_length())
        if max(bits_a, terms + bits_q + bits_b) + 3 <= k * _PRIME_BITS:
            return Optional(Optional(candidate^))
        k *= 2
    return None
