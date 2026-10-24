"""Unit Bezout identities of univariate rational polynomials by multimodular reconstruction.

For coprime f and h over QQ of positive degrees m and n there are unique u
and v with u f + v h = 1, deg u < n and deg v < m: extended Euclid returns
exactly them. Modulo a prime that divides no denominator and neither leading
coefficient, f and h stay coprime unless the prime divides their resultant,
and then the images of u and v are the modular cofactors. Those are combined
by the Chinese remainder theorem and each coefficient is recovered as a
fraction; a candidate is returned only once u f + v h = 1 holds exactly, so
an unlucky prime or a premature reconstruction costs only time. Extended
Euclid over QQ makes every intermediate remainder's coefficients grow, while
here only the answer's do.
"""
from apn_mojo import Integer, Rational, isqrt
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .primes import _residue
from ._rational_groebner import _next_prime, _reconstruct
from ._modular_simd import _Prime32, _subtract32

comptime _BEZOUT_PRIMES = 4096
"""Primes tried before the caller falls back to extended Euclid over QQ."""
comptime _BEZOUT_UNLUCKY = 4
"""Consecutive primes modulo which the inputs share a factor before giving up."""


def _inverse(a: UInt64, p: UInt64) -> UInt64:
    if p < (UInt64(1) << 32):
        return _Prime32.of(p).power[1](a % p, p - 2)
    var result = UInt64(1)
    var base = a % p
    var e = p - 2
    while e:
        if e & 1:
            result = UInt64((UInt128(result) * UInt128(base)) % UInt128(p))
        base = UInt64((UInt128(base) * UInt128(base)) % UInt128(p))
        e >>= 1
    return result


def _trim(mut a: List[UInt64]):
    while len(a) and a[len(a) - 1] == 0:
        _ = a.pop()


def _image(f: DensePoly[QQ], p: UInt64, mut budget: Budget) raises PolyError -> Optional[List[UInt64]]:
    """f modulo p, or None when p divides a denominator."""
    var result = List[UInt64](capacity=len(f.coefficients))
    for i in range(len(f.coefficients)):
        budget.charge(1)
        ref value = f.coefficients[i].value
        var d = _residue(value.denominator(), p, budget)
        if d == 0:
            return None
        result.append(UInt64((UInt128(_residue(value.numerator(), p, budget)) * UInt128(_inverse(d, p))) % UInt128(p)))
    return result^


def _mul_sub(a: List[UInt64], q: List[UInt64], b: List[UInt64], p: UInt64) -> List[UInt64]:
    """a - q b modulo p."""
    if p < (UInt64(1) << 32):
        return _mul_sub_in[True](a, q, b, _Prime32.of(p))
    return _mul_sub_in[False](a, q, b, _Prime32(p, 0))


def _mul_sub_in[small: Bool](a: List[UInt64], q: List[UInt64], b: List[UInt64], field: _Prime32) -> List[UInt64]:
    var p = field.p
    var n = max(len(a), len(q) + len(b) - 1 if len(q) and len(b) else 0)
    var result = List[UInt64](length=n, fill=0)
    for i in range(len(a)):
        result[i] = a[i]
    for i in range(len(q)):
        if q[i] == 0:
            continue
        comptime if small:
            _subtract32(result.unsafe_ptr().unsafe_offset(i).unsafe_origin_cast[MutUntrackedOrigin](), b.unsafe_ptr(), len(b), q[i], field)
        else:
            for j in range(len(b)):
                var t = UInt64((UInt128(q[i]) * UInt128(b[j])) % UInt128(p))
                result[i + j] = (result[i + j] + p - t) % p
    _trim(result)
    return result^


def _modular_bezout(f: List[UInt64], h: List[UInt64], p: UInt64) -> Optional[Tuple[List[UInt64], List[UInt64]]]:
    """Cofactors u, v with u f + v h = 1 modulo p, or None when f and h share a factor there."""
    if p < (UInt64(1) << 32):
        return _modular_bezout_in[True](f, h, _Prime32.of(p))
    return _modular_bezout_in[False](f, h, _Prime32(p, 0))


def _modular_bezout_in[small: Bool](f: List[UInt64], h: List[UInt64], field: _Prime32) -> Optional[Tuple[List[UInt64], List[UInt64]]]:
    var p = field.p
    var r0 = f.copy()
    var r1 = h.copy()
    var u0: List[UInt64] = [UInt64(1)]
    var v0 = List[UInt64]()
    var u1 = List[UInt64]()
    var v1: List[UInt64] = [UInt64(1)]
    while len(r1):
        # q, r with r0 = q r1 + r.
        var r = r0.copy()
        var dq = len(r0) - len(r1)
        var q = List[UInt64](length=max(dq + 1, 0), fill=0)
        var lead = _inverse(r1[len(r1) - 1], p)
        while len(r) >= len(r1) and len(r):
            var shift = len(r) - len(r1)
            var c: UInt64
            comptime if small:
                c = field.mul[1](r[len(r) - 1], lead)
            else:
                c = UInt64((UInt128(r[len(r) - 1]) * UInt128(lead)) % UInt128(p))
            q[shift] = c
            comptime if small:
                _subtract32(r.unsafe_ptr().unsafe_offset(shift).unsafe_origin_cast[MutUntrackedOrigin](), r1.unsafe_ptr(), len(r1), c, field)
            else:
                for j in range(len(r1)):
                    var t = UInt64((UInt128(c) * UInt128(r1[j])) % UInt128(p))
                    r[shift + j] = (r[shift + j] + p - t) % p
            _trim(r)
        _trim(q)
        var u = _mul_sub_in[small](u0, q, u1, field)
        var v = _mul_sub_in[small](v0, q, v1, field)
        r0 = r1^
        r1 = r^
        u0 = u1^
        v0 = v1^
        u1 = u^
        v1 = v^
    if len(r0) != 1:
        return None
    var scale = _inverse(r0[0], p)
    for i in range(len(u0)):
        comptime if small:
            u0[i] = field.mul[1](u0[i], scale)
        else:
            u0[i] = UInt64((UInt128(u0[i]) * UInt128(scale)) % UInt128(p))
    for i in range(len(v0)):
        comptime if small:
            v0[i] = field.mul[1](v0[i], scale)
        else:
            v0[i] = UInt64((UInt128(v0[i]) * UInt128(scale)) % UInt128(p))
    return (u0^, v0^)


def _unit_bezout(f: DensePoly[QQ], h: DensePoly[QQ], mut budget: Budget) raises PolyError -> Optional[Tuple[DensePoly[QQ], DensePoly[QQ]]]:
    """The unique u, v with u f + v h = 1, deg u < deg h, deg v < deg f, or None.

    None means only that this method did not find them; the caller decides by
    extended Euclid. Both inputs have positive degree in one variable.
    """
    var m = f.degree()
    var n = h.degree()
    if m <= 0 or n <= 0:
        return None
    # Residues of u's n and v's m coefficients, lifted by Garner's step.
    var count = n + m
    _allocation[Integer](count, budget)
    var values = List[Integer](length=count, fill=Integer(0))
    var modulus = Integer(1)
    var primes = 0
    var unlucky = 0
    var attempt = 1
    var p = (UInt64(1) << 31) + 1
    for _ in range(_BEZOUT_PRIMES):
        p = _next_prime(p, budget)
        var fp = _image(f, p, budget)
        var hp = _image(h, p, budget)
        if not fp or not hp or fp.value()[m] == 0 or hp.value()[n] == 0:
            continue
        var pair = _modular_bezout(fp.value(), hp.value(), p)
        if not pair:
            unlucky += 1
            if unlucky >= _BEZOUT_UNLUCKY:
                return None
            continue
        unlucky = 0
        ref up = pair.value()[0]
        ref vp = pair.value()[1]
        budget.charge(UInt64(count))
        var inverse = _inverse(_residue(modulus, p, budget), p)
        for k in range(count):
            var r: UInt64
            if k < n:
                r = up[k] if k < len(up) else UInt64(0)
            else:
                r = vp[k - n] if k - n < len(vp) else UInt64(0)
            var x = _residue(values[k], p, budget)
            var t = UInt64((UInt128((r + p - x) % p) * UInt128(inverse)) % UInt128(p))
            if t:
                try:
                    values[k] = values[k] + modulus * Integer(Int(t))
                except e:
                    raise PolyError(PolyErrorKind.CAPACITY, String(e))
        try:
            modulus = modulus * Integer(Int(p))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        primes += 1
        if primes < attempt:
            continue
        attempt = 2 * primes
        var bound: Integer
        try:
            bound = isqrt(modulus // Integer(2))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        var fractions = List[QQ](capacity=count)
        var complete = True
        for k in range(count):
            var fraction = _reconstruct(values[k], modulus, bound, budget)
            if not fraction:
                complete = False
                break
            fractions.append(QQ(fraction.take()))
        if not complete:
            continue
        var u_coefficients = List[QQ](capacity=n)
        for k in range(n):
            u_coefficients.append(fractions[k].copy())
        var v_coefficients = List[QQ](capacity=m)
        for k in range(m):
            v_coefficients.append(fractions[n + k].copy())
        var u = DensePoly[QQ].from_coefficients(f.variable, NoContext(), u_coefficients^, budget)
        var v = DensePoly[QQ].from_coefficients(f.variable, NoContext(), v_coefficients^, budget)
        var combination = u.mul(f, budget).add(v.mul(h, budget), budget)
        if combination.degree() == 0 and combination.coefficients[0] == QQ.one(NoContext()):
            return (u^, v^)
    return None
