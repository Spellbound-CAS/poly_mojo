"""Candidate filters for algebraic-coefficient roots.

A candidate root b of the coefficient annihilator is rejected when interval
evaluation of sum c_i b^i excludes zero on refined witnesses, and accepted
when every coefficient is represented as q_i(b) for rational polynomials q_i
with sum q_i(y) y^i divisible by the minimal polynomial of b. Both proofs use
polynomial identities and interval containments only; no arithmetic
annihilators are formed.
"""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .sparse import SparsePoly
from .algebraic_types import RealAlgebraic, _copy_real, _same_root, _rational, _check_real
from ._algebraic import _refine, _root_cached, _NormalizationCache, _Frame, _qq_gcd
from .descartes import _unique_root
from ._root_ops import _sign_at, _twos, _power_of_two
from .monomial import MonomialOrder
from ._dense_ops import _copy as _copy_poly, _remainder
from ._euclidean import _extended
from ._root_ops import _compare as _compare_q
from .certificates import _dense_equal, ResultantResult
from .division import divide_exact
from .resultant import resultant_value
from ._images import _denominator, _integer_image, _rational_image
from apn_mojo import Integer, Rational, gcd as integer_gcd
from .coefficients import ZZ
from ._coefficient_ops import _gcd_coefficient
from .primes import prime, _residue, _merge_residue
from ._rational_groebner import _reconstruct
from ._nmod import _modulus, _mulmod as _nmod_mulmod, _Field, _field, _mulm, _to, _from, _invm, _scale, _divmod_classic, _trim as _nmod_trim, _sub as _nmod_sub, _mul as _nmod_mul


# ------------------------------------------------------------------ intervals

def _bounds(a: RealAlgebraic, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    budget.charge(1)
    return (a.low.copy(), a.high.copy())


def _interval_mul(a: Tuple[QQ, QQ], b: Tuple[QQ, QQ], mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    var products = [a[0].mul(b[0], budget), a[0].mul(b[1], budget), a[1].mul(b[0], budget), a[1].mul(b[1], budget)]
    var low = products[0].copy()
    var high = products[0].copy()
    for i in range(1, 4):
        budget.charge(1)
        if _compare_q(products[i], low, budget) < 0:
            low = products[i].copy()
        if _compare_q(products[i], high, budget) > 0:
            high = products[i].copy()
    return (low^, high^)


def _interval_add(a: Tuple[QQ, QQ], b: Tuple[QQ, QQ], mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    return (a[0].add(b[0], budget), a[1].add(b[1], budget))


def _interval_value(coefficients: List[RealAlgebraic], x: RealAlgebraic, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    # Horner over closed intervals; every witness interval contains its value.
    var point = _bounds(x, budget)
    var value = (QQ.zero(NoContext()), QQ.zero(NoContext()))
    for i in range(len(coefficients) - 1, -1, -1):
        budget.charge(1)
        value = _interval_add(_interval_mul(value, point, budget), _bounds(coefficients[i], budget), budget)
    return value^


def _dense_interval(p: DensePoly[QQ], x: RealAlgebraic, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    # Horner over the closed interval of x, on integers over one positive
    # scale: the same bounds as rational interval Horner, normalized once.
    if len(p.coefficients) == 0:
        return (QQ.zero(NoContext()), QQ.zero(NoContext()))
    var cleared = _cleared(p, budget)
    ref c = cleared[0]
    var box = List[Tuple[QQ, QQ]]()
    box.append((x.low.copy(), x.high.copy()))
    var scaled = _scaled_box(box, budget)
    ref t = scaled[0]
    var n = len(c) - 1
    var value = (c[n], c[n])
    var scale = Integer(1)
    for i in range(n - 1, -1, -1):
        budget.charge(1)
        value = _scaled_mul(value, (t[0], t[1]), budget)
        scale = _zmul(scale, t[2], budget)
        var term = _zmul(c[i], scale, budget)
        value = (_zadd(value[0], term, budget), _zadd(value[1], term, budget))
    var denominator = _zmul(scale, cleared[1], budget)
    try:
        return (QQ(Rational(value[0], denominator)), QQ(Rational(value[1], denominator)))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _excludes(value: Tuple[QQ, QQ], mut budget: Budget) raises PolyError -> Bool:
    return value[0].value.sign() > 0 or value[1].value.sign() < 0


def _inside(value: Tuple[QQ, QQ], a: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    # Inside an isolating interval, or equal to a rational value. The interval's
    # endpoints are not roots of its definition, so a root of that definition in
    # the closed image interval lies in the open isolating interval.
    if a.rational:
        return _compare_q(value[0], a.low, budget) == 0 and _compare_q(value[1], a.low, budget) == 0
    return _compare_q(a.low, value[0], budget) <= 0 and _compare_q(value[1], a.high, budget) <= 0


# ------------------------------------------------------------------ the field Q(b) = Q[y]/(m)

def _reduce(p: DensePoly[QQ], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    if p.degree() < m.degree():
        return _copy_poly(p, budget)
    if m.degree() < 1:
        return _remainder(p, m, budget)
    # The remainder by integer pseudo-division, normalized once.
    return _field_mul(p, DensePoly[QQ].from_coefficients(p.variable, NoContext(), [QQ.one(NoContext())], budget), m, budget)


def _cleared(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> Tuple[List[Integer], Integer]:
    # p = P / d with integer coefficients P and the positive common denominator d.
    var common = Integer(1)
    for c in p.coefficients:
        budget.charge(1)
        if c.value.denominator() != 1:
            common = _denominator(common, c, budget)
    var scaled = List[Integer](capacity=len(p.coefficients))
    budget.charge(UInt64(len(p.coefficients)))
    try:
        for c in p.coefficients:
            var d = c.value.denominator()
            if d == common:
                scaled.append(c.value.numerator())
            else:
                var factor = common // d
                scaled.append(c.value.numerator() * factor)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return (scaled^, common^)


def _field_mul(a: DensePoly[QQ], b: DensePoly[QQ], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # a b modulo m over Q, computed on integers: clear denominators, multiply,
    # pseudo-divide by the integer multiple M of m, and normalize each
    # remainder coefficient once instead of after every rational operation.
    if a.variable != b.variable or a.is_zero() or b.is_zero() or m.degree() < 1:
        return _reduce(a.mul(b, budget), m, budget)
    var x = _cleared(a, budget)
    var y = _cleared(b, budget)
    var z = _cleared(m, budget)
    ref A = x[0]
    ref B = y[0]
    ref M = z[0]
    var n = len(A) + len(B) - 1
    var dm = len(M) - 1
    budget.charge(UInt64(_product(len(A), len(B))))
    if n - 1 >= dm:
        budget.charge(UInt64(_product(n - dm, len(M))))
    budget.charge(UInt64(n))
    var coefficients = List[QQ]()
    try:
        var product = List[Integer](length=n, fill=Integer(0))
        for i in range(len(A)):
            if not A[i]:
                continue
            for j in range(len(B)):
                product[i + j] = product[i + j] + A[i] * B[j]
        var lead = M[dm]
        var scale = x[1] * y[1]
        if n - 1 >= dm:
            for top in range(n - 1, dm - 1, -1):
                var t = product[top]
                if not t:
                    continue
                if lead != 1:
                    for k in range(top + 1):
                        product[k] = product[k] * lead
                    scale = scale * lead
                for k in range(dm + 1):
                    product[top - dm + k] = product[top - dm + k] - t * M[k]
            n = dm
        for k in range(n):
            coefficients.append(QQ(Rational(product[k], scale)))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return DensePoly[QQ].from_coefficients(a.variable, NoContext(), coefficients^, budget)


def _field_inverse(a: DensePoly[QQ], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    var e = _extended(a, m, budget)
    if e[0].degree() != 0:
        raise PolyError(PolyErrorKind.DOMAIN, "element is not invertible modulo the minimal polynomial")
    # _extended returns a monic gcd, so s*a + t*m == 1.
    return _reduce(e[1], m, budget)


def _compose(f: DensePoly[QQ], q: DensePoly[QQ], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # f(q(y)) modulo m(y), by Horner on integers: the accumulator is P / s
    # with P reduced modulo the integer multiple M of m, put in lowest terms
    # once per step instead of normalizing every rational coefficient.
    if m.degree() < 1 or q.variable != 0 or f.is_zero() or q.is_zero():
        return _compose_rational(f, q, m, budget)
    var x = _cleared(f, budget)
    var y = _cleared(q, budget)
    var z = _cleared(m, budget)
    ref F = x[0]
    ref Q = y[0]
    ref M = z[0]
    var P = List[Integer]()
    var s = Integer(1)
    for i in range(len(F) - 1, -1, -1):
        budget.charge(1)
        budget.charge(UInt64(_product(_sum(len(P), 1), len(Q))))
        budget.charge(UInt64(_product(_sum(len(P), len(Q)), len(M))))
        var reduced = _reduced_product(P, Q, M, budget)
        try:
            # (P / s) q + F_i over the common denominator s c_q L^e.
            var scale = s * y[1] * reduced[1]
            var next = reduced[0].copy()
            if len(next) == 0:
                next.append(Integer(0))
            next[0] = next[0] + F[i] * scale
            while len(next) and not next[len(next) - 1]:
                _ = next.pop()
            var lowest = _lowest_terms(next, scale, budget)
            P = lowest[0].copy()
            s = lowest[1].copy()
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var coefficients = List[QQ]()
    budget.charge(UInt64(len(P)))
    try:
        var denominator = s * x[1]
        for k in range(len(P)):
            coefficients.append(QQ(Rational(P[k], denominator)))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)


def _compose_rational(f: DensePoly[QQ], q: DensePoly[QQ], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    var result = DensePoly[QQ](0, NoContext())
    for i in range(len(f.coefficients) - 1, -1, -1):
        budget.charge(1)
        var constant = DensePoly[QQ].from_coefficients(0, NoContext(), [f.coefficients[i].copy()], budget)
        result = _field_mul(result, q, m, budget).add(constant, budget)
    return result^


def _combination(representations: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # sum q_i(y) y^i modulo m(y).
    var y = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
    var result = DensePoly[QQ](0, NoContext())
    for i in range(len(representations) - 1, -1, -1):
        budget.charge(1)
        result = _field_mul(result, y, m, budget).add(representations[i], budget)
    return _reduce(result, m, budget)


# ------------------------------------------------------------------ checks

def _check_exclusion(coefficients: List[RealAlgebraic], candidate: RealAlgebraic, refined: List[RealAlgebraic], refined_candidate: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    if len(refined) != len(coefficients) or not _same_root(refined_candidate, candidate, budget):
        return False
    for i in range(len(coefficients)):
        budget.charge(1)
        if not _same_root(refined[i], coefficients[i], budget):
            return False
    return _excludes(_interval_value(refined, refined_candidate, budget), budget)


def _check_representation(coefficients: List[RealAlgebraic], candidate: RealAlgebraic, refined: List[RealAlgebraic], refined_candidate: RealAlgebraic, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    var n = len(coefficients)
    if len(refined) != n or len(representations) != n or not _same_root(refined_candidate, candidate, budget):
        return False
    ref m = candidate.polynomial
    for i in range(n):
        budget.charge(1)
        ref q = representations[i]
        q.check_invariants(budget)
        if q.variable != 0 or not _same_root(refined[i], coefficients[i], budget):
            return False
        if coefficients[i].rational:
            # q_i(b) equals the rational coefficient when m divides q_i - c_i.
            var constant = DensePoly[QQ].from_coefficients(0, NoContext(), [coefficients[i].low.copy()], budget)
            if not _reduce(q.sub(constant, budget), m, budget).is_zero():
                return False
        else:
            # q_i(b) is a root of the coefficient's definition inside its
            # isolating interval, hence equal to the coefficient.
            if not _compose(coefficients[i].polynomial, q, m, budget).is_zero() or not _inside(_dense_interval(q, refined_candidate, budget), refined[i], budget):
                return False
    return _combination(representations, m, budget).is_zero()


# ------------------------------------------------------------------ producers

def _try_exclusion(coefficients: List[RealAlgebraic], candidate: RealAlgebraic, rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], RealAlgebraic]]:
    var levels = _PointLevels(List[RealAlgebraic](), budget)
    return _try_exclusion_levels(coefficients, candidate, levels, rounds, budget)


def _try_exclusion_levels(coefficients: List[RealAlgebraic], candidate: RealAlgebraic, mut levels: _PointLevels, rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], RealAlgebraic]]:
    # As _try_exclusion, with levels kept between calls: the coefficients'
    # bisections serve every candidate of one polynomial.
    # With rational coefficients the norm is the input itself, so every
    # candidate is a root and no interval can exclude it.
    var irrational = False
    for c in coefficients:
        budget.charge(1)
        if not c.rational:
            irrational = True
            break
    if not irrational:
        return None
    # The first round r at which Horner's interval value over the coefficients
    # and the candidate, each bisected r times, excludes zero. Horner over
    # intervals is inclusion-monotone and the bisections are nested, so the
    # rounds that exclude zero form a suffix: gallop, then bisect the gap.
    var values = List[RealAlgebraic]()
    for c in coefficients:
        values.append(_copy_real(c, budget))
    values.append(_copy_real(candidate, budget))
    levels.reset(values, budget)
    var failed = -1
    var found = -1
    var step = 1
    var r = 0
    while r < rounds:
        budget.charge(1)
        _ = _has_level(levels, r, budget)
        if _horner_excludes(_level_box(levels, r, budget), budget):
            found = r
            break
        failed = r
        if r == rounds - 1:
            break
        r = min(r + step, rounds - 1)
        step *= 2
    if found < 0:
        return None
    while found - failed > 1:
        budget.charge(1)
        var middle = failed + (found - failed) // 2
        if _horner_excludes(_level_box(levels, middle, budget), budget):
            found = middle
        else:
            failed = middle
    var point = _level_point(levels, found, budget)
    var refined_candidate = point.pop()
    return (point^, refined_candidate^)


def _horner_excludes(box: List[Tuple[QQ, QQ]], mut budget: Budget) raises PolyError -> Bool:
    # Whether Horner's rule over closed intervals, with box[:-1] the
    # coefficient intervals and box[-1] the argument, gives an interval
    # without zero; the bounds are kept as integers over one positive scale.
    var scaled = _scaled_box(box, budget)
    var n = len(scaled) - 1
    ref x = scaled[n]
    var value = (Integer(0), Integer(0))
    var scale = Integer(1)
    for i in range(n - 1, -1, -1):
        budget.charge(1)
        value = _scaled_mul(value, (x[0], x[1]), budget)
        scale = _zmul(scale, x[2], budget)
        ref c = scaled[i]
        value = (_zadd(_zmul(value[0], c[2], budget), _zmul(c[0], scale, budget), budget), _zadd(_zmul(value[1], c[2], budget), _zmul(c[1], scale, budget), budget))
        scale = _zmul(scale, c[2], budget)
    return value[0].sign() > 0 or value[1].sign() < 0


def _try_containment(coefficients: List[RealAlgebraic], candidate: RealAlgebraic, representations: List[DensePoly[QQ]], rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], RealAlgebraic]]:
    var refined = List[RealAlgebraic]()
    for c in coefficients:
        budget.charge(1)
        refined.append(_copy_real(c, budget))
    var point = _copy_real(candidate, budget)
    for _ in range(rounds):
        budget.charge(1)
        var inside = True
        for i in range(len(refined)):
            if not refined[i].rational and not _inside(_dense_interval(representations[i], point, budget), refined[i], budget):
                inside = False
                break
        if inside:
            return (refined^, point^)
        _refine(point, budget)
    return None


def _coefficient_polynomials(f: SparsePoly[QQ], variable: Int64, coordinate: Int64, point: List[Tuple[Int64, RealAlgebraic]], size: Int, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    # The coefficient of variable^i as a dense polynomial in the one
    # irrational coordinate, after substituting every rational coordinate.
    var rows = List[List[QQ]]()
    for _ in range(size):
        rows.append(List[QQ]())
    for term in f.terms:
        budget.charge(1)
        var power = 0
        var degree = 0
        var value = term.coefficient.copy()
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] == variable:
                power = Int(pair[1])
            elif pair[0] == coordinate:
                degree = Int(pair[1])
            else:
                var found = False
                for entry in point:
                    if entry[0] == pair[0]:
                        if not entry[1].rational:
                            raise PolyError(PolyErrorKind.DOMAIN, "a second irrational coordinate")
                        for _ in range(Int(pair[1])):
                            value = value.mul(entry[1].low, budget)
                        found = True
                        break
                if not found:
                    raise PolyError(PolyErrorKind.DOMAIN, "missing coordinate")
        if power >= size:
            raise PolyError(PolyErrorKind.DOMAIN, "coefficient degree mismatch")
        while len(rows[power]) <= degree:
            rows[power].append(QQ.zero(NoContext()))
        rows[power][degree] = rows[power][degree].add(value, budget)
    var result = List[DensePoly[QQ]]()
    for i in range(size):
        budget.charge(1)
        result.append(DensePoly[QQ].from_coefficients(0, NoContext(), rows[i].copy(), budget))
    return result^


def _gcd_over_field(a: List[DensePoly[QQ]], b: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget, *, prime_limit: Int = 512) raises PolyError -> List[DensePoly[QQ]]:
    # Field inverses swell rational coefficients; modular images avoid them
    # and the exact division checks keep the result exact.
    if m.degree() >= 2:
        var fast = _gcd_modular(a, b, m, budget, prime_limit=prime_limit)
        if fast:
            return fast.take()
    return _gcd_euclid(a, b, m, budget)


def _gcd_euclid(a: List[DensePoly[QQ]], b: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    # Monic Euclid in Q(b)[x]; polynomials are coefficient lists, lowest degree first.
    var f = _trim(a, m, budget)
    var g = _trim(b, m, budget)
    while len(g):
        budget.charge(1)
        var inverse = _field_inverse(g[len(g) - 1], m, budget)
        var r = f.copy()
        while len(r) >= len(g):
            budget.charge(1)
            var factor = _field_mul(r[len(r) - 1], inverse, m, budget)
            var shift = len(r) - len(g)
            for i in range(len(g)):
                r[shift + i] = _reduce(r[shift + i].sub(_field_mul(factor, g[i], m, budget), budget), m, budget)
            r = _trim(r, m, budget)
            if len(r) == 0:
                break
        f = g^
        g = r^
    if len(f):
        var inverse = _field_inverse(f[len(f) - 1], m, budget)
        for i in range(len(f)):
            f[i] = _field_mul(f[i], inverse, m, budget)
    return f^


# Modular images of the monic gcd in Q(b)[x]: Euclid modulo word primes in
# (F_p[y]/m)[x], CRT and rational reconstruction, then exact division checks.
comptime _GCD_PRIMES = 512
"""Primes tried before the exact Euclidean fallback."""


def _residue_qq(q: QQ, f: _Field, mut budget: Budget) raises PolyError -> Optional[UInt64]:
    var den = _residue(q.value.denominator(), f.p, budget)
    if den == 0:
        return None
    var num = _residue(q.value.numerator(), f.p, budget)
    return _mulm(_to(num, f), _invm(_to(den, f), f, budget), f)


def _image_words(p: DensePoly[QQ], f: _Field, mut budget: Budget) raises PolyError -> Optional[List[UInt64]]:
    _allocation[UInt64](len(p.coefficients), budget)
    var result = List[UInt64](capacity=len(p.coefficients))
    for c in p.coefficients:
        var r = _residue_qq(c, f, budget)
        if not r:
            return None
        result.append(r.value())
    _nmod_trim(result)
    return result^


def _inverse_words(a: List[UInt64], modulus: List[UInt64], f: _Field, mut budget: Budget) raises PolyError -> Optional[List[UInt64]]:
    # Extended Euclid in F_p[y]; None when a is a zero divisor modulo the modulus.
    var r0 = modulus.copy()
    var r1 = a.copy()
    var t0 = List[UInt64]()
    var t1 = List[UInt64]()
    t1.append(f.one)
    while len(r1):
        budget.charge(1)
        var qr = _divmod_classic(r0, r1, True, f, budget)
        var t2 = _nmod_sub(t0, _nmod_mul(qr[0], t1, f, budget), f, budget)
        r0 = r1^
        r1 = qr[1].copy()
        t0 = t1^
        t1 = t2^
    if len(r0) != 1:
        return None
    return _divmod_classic(_scale(t0, _invm(r0[0], f, budget), f, budget), modulus, False, f, budget)[1].copy()


def _gcd_image(a: List[List[UInt64]], b: List[List[UInt64]], modulus: List[UInt64], f: _Field, mut budget: Budget) raises PolyError -> Optional[List[List[UInt64]]]:
    # Monic Euclid in (F_p[y]/modulus)[x], as _gcd_euclid; None on a zero divisor.
    var x = a.copy()
    var y = b.copy()
    # Products reduce by the modulus's precomputed inverse, not by long division.
    var reduction = _modulus(modulus, f, budget)
    while len(y):
        budget.charge(1)
        var inverse = _inverse_words(y[len(y) - 1], modulus, f, budget)
        if not inverse:
            return None
        var r = x.copy()
        while len(r) >= len(y):
            budget.charge(1)
            var factor = _nmod_mulmod(r[len(r) - 1], inverse.value(), reduction, f, budget)
            var shift = len(r) - len(y)
            for i in range(len(y)):
                r[shift + i] = _nmod_sub(r[shift + i], _nmod_mulmod(factor, y[i], reduction, f, budget), f, budget)
            while len(r) and len(r[len(r) - 1]) == 0:
                _ = r.pop()
        x = y^
        y = r^
    if len(x):
        var inverse = _inverse_words(x[len(x) - 1], modulus, f, budget)
        if not inverse:
            return None
        for i in range(len(x)):
            x[i] = _nmod_mulmod(x[i], inverse.value(), reduction, f, budget)
    return x^


def _half_bound(modulus: Integer) raises PolyError -> Integer:
    # 2 * bound^2 < modulus, so a reconstructed fraction is unique.
    try:
        return Integer(1) << Integer((modulus.magnitude_bit_length() - 2) // 2)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _gcd_modular(a: List[DensePoly[QQ]], b: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget, *, prime_limit: Int = _GCD_PRIMES) raises PolyError -> Optional[List[DensePoly[QQ]]]:
    var x = _trim(a, m, budget)
    var y = _trim(b, m, budget)
    if len(x) == 0 or len(y) == 0:
        return None
    var d = m.degree()
    var degree = -1
    var values = List[Integer]()
    var modulus = Integer(1)
    var probe = 0
    for index in range(prime_limit):
        budget.charge(1)
        var f = _field(prime(index, budget))
        var mm = _image_words(m, f, budget)
        if not mm or len(mm.value()) != _sum(d, 1):
            continue
        var monic = _scale(mm.value(), _invm(mm.value()[d], f, budget), f, budget)
        var xi = List[List[UInt64]]()
        var yi = List[List[UInt64]]()
        var usable = True
        for c in x:
            var w = _image_words(c, f, budget)
            if not w:
                usable = False
                break
            xi.append(w.value().copy())
        for c in y:
            if not usable:
                break
            var w = _image_words(c, f, budget)
            if not w:
                usable = False
                break
            yi.append(w.value().copy())
        if not usable or len(xi[len(xi) - 1]) == 0 or len(yi[len(yi) - 1]) == 0:
            continue
        var image = _gcd_image(xi, yi, monic, f, budget)
        if not image:
            continue
        ref g = image.value()
        var k = len(g) - 1
        if degree >= 0 and k > degree:
            continue
        if degree < 0 or k < degree:
            degree = k
            values = List[Integer](length=_product(k, d), fill=Integer(0))
            modulus = Integer(1)
            probe = 0
        for i in range(k):
            for j in range(d):
                budget.charge(1)
                var r = _from(g[i][j], f) if j < len(g[i]) else UInt64(0)
                values[i * d + j] = _merge_residue(values[i * d + j], modulus, r, f.p, budget)
        modulus = ZZ(modulus^).mul(ZZ(Integer(Int(f.p))), budget).value
        var bound = _half_bound(modulus)
        # The entry that failed last stops further work until it reconstructs:
        # one reconstruction per prime instead of all of them at doubling counts.
        if probe < len(values) and not _reconstruct(values[probe], modulus, bound, budget):
            continue
        var candidate = List[DensePoly[QQ]]()
        var complete = True
        for i in range(k):
            var entries = List[QQ]()
            for j in range(d):
                var q = _reconstruct(values[i * d + j], modulus, bound, budget)
                if not q:
                    probe = i * d + j
                    complete = False
                    break
                entries.append(QQ(q.value().copy()))
            if not complete:
                break
            candidate.append(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget))
        if not complete:
            continue
        candidate.append(DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], budget))
        if _divides(candidate, x, m, budget) and _divides(candidate, y, m, budget):
            return _trim(candidate, m, budget)
    return None


def _trim(p: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    var result = List[DensePoly[QQ]]()
    for c in p:
        result.append(_reduce(c, m, budget))
    while len(result) and result[len(result) - 1].is_zero():
        _ = result.pop()
    return result^


def _represent(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], coefficients: List[RealAlgebraic], candidate: RealAlgebraic, mut budget: Budget) raises PolyError -> Optional[List[DensePoly[QQ]]]:
    """Polynomials q_i with c_i = q_i(b), when b determines the point's one irrational coordinate."""
    ref m = candidate.polynomial
    var irrational = -1
    for i in range(len(point)):
        if not point[i][1].rational:
            if irrational >= 0:
                return None
            irrational = i
    var result = List[DensePoly[QQ]]()
    if irrational < 0:
        for c in coefficients:
            result.append(DensePoly[QQ].from_coefficients(0, NoContext(), [c.low.copy()], budget))
        return result^
    var coordinate = point[irrational][0]
    var rows = _coefficient_polynomials(f, variable, coordinate, point, len(coefficients), budget)
    # F(x, b) as a polynomial in x over Q(b): coefficient of x^k is sum_i row_i[k] y^i.
    var width = 0
    for r in rows:
        width = max(width, len(r.coefficients))
    var specialized = List[DensePoly[QQ]]()
    for k in range(width):
        budget.charge(1)
        var entries = List[QQ]()
        for i in range(len(rows)):
            entries.append(rows[i].coefficients[k].copy() if k < len(rows[i].coefficients) else QQ.zero(NoContext()))
        specialized.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), m, budget))
    var definition = List[DensePoly[QQ]]()
    ref alpha = point[irrational][1]
    for c in alpha.polynomial.coefficients:
        definition.append(DensePoly[QQ].from_coefficients(0, NoContext(), [c.copy()], budget))
    var common = _gcd_over_field(definition, specialized, m, budget)
    if len(common) != 2:
        return None
    # x = -common[0] since common is monic of degree one.
    var a = DensePoly[QQ](0, NoContext()).sub(common[0], budget)
    for i in range(len(rows)):
        budget.charge(1)
        result.append(_compose(rows[i], a, m, budget))
    if not _combination(result, m, budget).is_zero():
        return None
    return result^


def _represent_element(element: RealAlgebraic, representations: List[DensePoly[QQ]], candidate: RealAlgebraic, mut budget: Budget) raises PolyError -> Optional[List[DensePoly[QQ]]]:
    """Polynomials q_i with c_i = q_i(b) from c_i = r_i(element), when b determines the element."""
    if candidate.rational or element.rational:
        return None
    ref me = element.polynomial
    ref mb = candidate.polynomial
    var width = 0
    for r in representations:
        width = max(width, len(r.coefficients))
    var specialized = List[DensePoly[QQ]]()
    for k in range(width):
        budget.charge(1)
        var entries = List[QQ]()
        for i in range(len(representations)):
            entries.append(representations[i].coefficients[k].copy() if k < len(representations[i].coefficients) else QQ.zero(NoContext()))
        specialized.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), mb, budget))
    var common = _gcd_over_field(_constants(me, budget), specialized, mb, budget)
    if len(common) != 2:
        return None
    var g = DensePoly[QQ](0, NoContext()).sub(common[0], budget)
    var result = List[DensePoly[QQ]]()
    for r in representations:
        budget.charge(1)
        result.append(_compose(r, g, mb, budget))
    if not _combination(result, mb, budget).is_zero():
        return None
    return result^


def _irrational_count(point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> Int:
    var count = 0
    for entry in point:
        budget.charge(1)
        if not entry[1].rational:
            count += 1
    return count


# ------------------------------------------------------------------ signs at points

def _scaled_mul(a: Tuple[Integer, Integer], b: Tuple[Integer, Integer], mut budget: Budget) raises PolyError -> Tuple[Integer, Integer]:
    # The interval product min/max by sign cases: the same bounds as the four
    # products, with two multiplications unless both intervals contain zero.
    if b[0].sign() >= 0:
        if a[0].sign() >= 0:
            return (_zmul(a[0], b[0], budget), _zmul(a[1], b[1], budget))
        if a[1].sign() <= 0:
            return (_zmul(a[0], b[1], budget), _zmul(a[1], b[0], budget))
        return (_zmul(a[0], b[1], budget), _zmul(a[1], b[1], budget))
    if b[1].sign() <= 0:
        if a[0].sign() >= 0:
            return (_zmul(a[1], b[0], budget), _zmul(a[0], b[1], budget))
        if a[1].sign() <= 0:
            return (_zmul(a[1], b[1], budget), _zmul(a[0], b[0], budget))
        return (_zmul(a[1], b[0], budget), _zmul(a[0], b[0], budget))
    if a[0].sign() >= 0:
        return (_zmul(a[1], b[0], budget), _zmul(a[1], b[1], budget))
    if a[1].sign() <= 0:
        return (_zmul(a[0], b[1], budget), _zmul(a[0], b[0], budget))
    var low = _zmul(a[0], b[1], budget)
    var other = _zmul(a[1], b[0], budget)
    if other < low:
        low = other^
    var high = _zmul(a[0], b[0], budget)
    var top = _zmul(a[1], b[1], budget)
    if top > high:
        high = top^
    return (low^, high^)


def _zmul(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    return ZZ(a).mul(ZZ(b), budget).value


def _zadd(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    return ZZ(a).add(ZZ(b), budget).value


def _scaled_box(box: List[Tuple[QQ, QQ]], mut budget: Budget) raises PolyError -> List[Tuple[Integer, Integer, Integer]]:
    # Each closed interval [l, h] as (L, H, D) with l = L/D, h = H/D, D > 0.
    var scaled = List[Tuple[Integer, Integer, Integer]]()
    for x in box:
        budget.charge(1)
        var ld = x[0].value.denominator()
        var hd = x[1].value.denominator()
        if ld == hd:
            scaled.append((x[0].value.numerator(), x[1].value.numerator(), ld^))
        else:
            scaled.append((_zmul(x[0].value.numerator(), hd, budget), _zmul(x[1].value.numerator(), ld, budget), _zmul(ld, hd, budget)))
    return scaled^


def _endpoints(point: List[RealAlgebraic], mut budget: Budget) raises PolyError -> List[Tuple[QQ, QQ]]:
    var box = List[Tuple[QQ, QQ]]()
    for x in point:
        budget.charge(1)
        box.append((x.low.copy(), x.high.copy()))
    return box^


def _interval_signs(g: SparsePoly[ZZ], variables: List[Int64], box: List[Tuple[QQ, QQ]], mut budget: Budget) raises PolyError -> Tuple[Int, Int]:
    # The signs of the bounds that term-wise interval arithmetic over the box
    # gives for g: each term's power products in order, then the sum. Every
    # bound is scaled by the same positive integer, so only signs are exact
    # here, and those are all the witnesses use.
    var scaled = _scaled_box(box, budget)
    var top = List[Int](length=len(variables), fill=0)
    var indices = List[List[Int]]()
    for term in g.terms:
        budget.charge(1)
        var row = List[Int]()
        for pair in term.monomial.pairs:
            budget.charge(1)
            var index = -1
            for i in range(len(variables)):
                if variables[i] == pair[0]:
                    index = i
                    break
            if index < 0:
                raise PolyError(PolyErrorKind.DOMAIN, "missing coordinate")
            row.append(index)
            if Int(pair[1]) > top[index]:
                top[index] = Int(pair[1])
        indices.append(row^)
    # powers[i][e] = D_i^e for e up to the largest exponent of coordinate i.
    var powers = List[List[Integer]]()
    for i in range(len(variables)):
        var row = List[Integer]()
        row.append(Integer(1))
        for e in range(top[i]):
            budget.charge(1)
            row.append(_zmul(row[e], scaled[i][2], budget))
        powers.append(row^)
    var low = Integer(0)
    var high = Integer(0)
    for t in range(len(g.terms)):
        budget.charge(1)
        ref term = g.terms[t]
        var value = (term.coefficient.value, term.coefficient.value)
        var used = List[Int](length=len(variables), fill=0)
        for j in range(len(term.monomial.pairs)):
            var index = indices[t][j]
            var e = Int(term.monomial.pairs[j][1])
            used[index] = e
            for _ in range(e):
                budget.charge(1)
                value = _scaled_mul(value, (scaled[index][0], scaled[index][1]), budget)
        var factor = Integer(1)
        for i in range(len(variables)):
            if top[i] > used[i]:
                factor = _zmul(factor, powers[i][top[i] - used[i]], budget)
        low = _zadd(low, _zmul(value[0], factor, budget), budget)
        high = _zadd(high, _zmul(value[1], factor, budget), budget)
    return (low.sign(), high.sign())


def _constant_compose(f: SparsePoly[QQ], variables: List[Int64], representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # _sparse_compose when every representation is a constant: the value of f at a rational point.
    var values = List[QQ]()
    for r in representations:
        values.append(r.coefficients[0].copy() if r.degree() == 0 else QQ.zero(NoContext()))
    var total = QQ.zero(NoContext())
    for term in f.terms:
        budget.charge(1)
        var value = term.coefficient.copy()
        for pair in term.monomial.pairs:
            budget.charge(1)
            var index = -1
            for i in range(len(variables)):
                if variables[i] == pair[0]:
                    index = i
                    break
            if index < 0:
                raise PolyError(PolyErrorKind.DOMAIN, "missing coordinate")
            for _ in range(Int(pair[1])):
                value = value.mul(values[index], budget)
        total = total.add(value, budget)
    return DensePoly[QQ].from_coefficients(0, NoContext(), [total^], budget)


def _sparse_compose(f: SparsePoly[QQ], variables: List[Int64], representations: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # f(r_1(y), ..., r_k(y)) modulo m(y).
    if m.degree() >= 1:
        var constant = True
        for r in representations:
            if r.degree() > 0:
                constant = False
                break
        if constant:
            return _constant_compose(f, variables, representations, budget)
    var total = DensePoly[QQ](0, NoContext())
    for term in f.terms:
        budget.charge(1)
        var value = DensePoly[QQ].from_coefficients(0, NoContext(), [term.coefficient.copy()], budget)
        for pair in term.monomial.pairs:
            budget.charge(1)
            var index = -1
            for i in range(len(variables)):
                if variables[i] == pair[0]:
                    index = i
                    break
            if index < 0:
                raise PolyError(PolyErrorKind.DOMAIN, "missing coordinate")
            for _ in range(Int(pair[1])):
                value = _field_mul(value, representations[index], m, budget)
        total = total.add(value, budget)
    return _reduce(total, m, budget)


def _check_point_interval(f: SparsePoly[QQ], variables: List[Int64], point: List[RealAlgebraic], refined: List[RealAlgebraic], relation: Int, mut budget: Budget) raises PolyError -> Bool:
    if len(refined) != len(point) or relation == 0:
        return False
    for i in range(len(point)):
        budget.charge(1)
        if not _same_root(refined[i], point[i], budget):
            return False
    var value = _interval_signs(_integer_image(f, budget), variables, _endpoints(refined, budget), budget)
    return (relation > 0 and value[0] > 0) or (relation < 0 and value[1] < 0)


def _check_point_zero(f: SparsePoly[QQ], variables: List[Int64], point: List[RealAlgebraic], refined: List[RealAlgebraic], primitive: Int, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    var k = len(point)
    if len(refined) != k or len(representations) != k or primitive < 0 or primitive >= k or point[primitive].rational:
        return False
    ref m = point[primitive].polynomial
    var y = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
    for j in range(k):
        budget.charge(1)
        ref r = representations[j]
        r.check_invariants(budget)
        if r.variable != 0 or not _same_root(refined[j], point[j], budget):
            return False
        if j == primitive:
            if not _dense_equal(r, y, budget):
                return False
        elif point[j].rational:
            var constant = DensePoly[QQ].from_coefficients(0, NoContext(), [point[j].low.copy()], budget)
            if not _reduce(r.sub(constant, budget), m, budget).is_zero():
                return False
        elif not _compose(point[j].polynomial, r, m, budget).is_zero() or not _inside(_dense_interval(r, refined[primitive], budget), refined[j], budget):
            return False
    return _sparse_compose(f, variables, representations, m, budget).is_zero()


struct _PointLevels(Movable):
    """The bisections of one sample point, shared by every input at a cell.

    Each coordinate keeps its own levels, so a point that differs from the
    last one only in some coordinates keeps the others' bisections.
    """
    var coordinates: List[RealAlgebraic]
    """The point itself."""
    var bounds: List[List[Tuple[QQ, QQ]]]
    """bounds[i][r]: coordinate i's interval after r bisections."""
    var ends: List[Tuple[Int, Int]]
    """The defining polynomials' signs at each coordinate's last level's endpoints, 2 when unknown."""
    var exact: Bool
    """Whether every coordinate is rational, so that no level follows the first."""
    var cleared: List[List[Integer]]
    """Each irrational coordinate's defining polynomial with denominators cleared."""

    def __init__(out self, point: List[RealAlgebraic], mut budget: Budget) raises PolyError:
        self.coordinates = List[RealAlgebraic]()
        self.cleared = List[List[Integer]]()
        self.bounds = List[List[Tuple[QQ, QQ]]]()
        self.ends = List[Tuple[Int, Int]]()
        self.exact = True
        self.reset(point, budget)

    def reset(mut self, point: List[RealAlgebraic], mut budget: Budget) raises PolyError:
        """Make this the levels of point, keeping each unchanged coordinate's bisections."""
        while len(self.coordinates) > len(point):
            _ = self.coordinates.pop()
            _ = self.cleared.pop()
            _ = self.bounds.pop()
            _ = self.ends.pop()
        self.exact = True
        for i in range(len(point)):
            budget.charge(1)
            ref c = point[i]
            if not c.rational:
                self.exact = False
            if i < len(self.coordinates):
                ref old = self.coordinates[i]
                if old.rational == c.rational and _compare_q(old.low, c.low, budget) == 0 and _compare_q(old.high, c.high, budget) == 0 and _dense_equal(old.polynomial, c.polynomial, budget):
                    continue
                self.coordinates[i] = _copy_real(c, budget)
                self.cleared[i] = _cleared(c.polynomial, budget)[0].copy() if not c.rational else List[Integer]()
                var first = List[Tuple[QQ, QQ]]()
                first.append((c.low.copy(), c.high.copy()))
                self.bounds[i] = first^
                self.ends[i] = (2, 2)
            else:
                self.coordinates.append(_copy_real(c, budget))
                self.cleared.append(_cleared(c.polynomial, budget)[0].copy() if not c.rational else List[Integer]())
                var first = List[Tuple[QQ, QQ]]()
                first.append((c.low.copy(), c.high.copy()))
                self.bounds.append(first^)
                self.ends.append((2, 2))


def _mid(a: QQ, b: QQ, mut budget: Budget) raises PolyError -> QQ:
    # A point strictly inside (a, b): the midpoint (a + b) / 2 for dyadic
    # endpoints, where only factors of two can cancel, so no gcd is needed.
    budget.charge(1)
    ref x = a.value
    ref y = b.value
    var da = x.denominator()
    var db = y.denominator()
    var ea = _power_of_two(da)
    var eb = _power_of_two(db)
    try:
        if ea >= 0 and eb >= 0:
            var k = max(ea, eb)
            var n = (x.numerator() << (k - ea)) + (y.numerator() << (k - eb))
            if not n:
                return QQ.zero(NoContext())
            var t = min(_twos(n), k + 1)
            return QQ(Rational(_numerator=n >> t, _denominator=Integer(1) << (k + 1 - t)))
        # Otherwise a dyadic point within a quarter width of the midpoint:
        # floor(mid 2^k) / 2^k with 2^-k <= (b - a) / 4, still strictly
        # inside, by one integer division instead of a gcd.
        var n = x.numerator() * db + y.numerator() * da
        var den = da * db
        var k = max((den << 2).magnitude_bit_length() - (y.numerator() * da - x.numerator() * db).magnitude_bit_length() + 1, 1)
        var p = (n << (k - 1)) // den
        if not p:
            return QQ.zero(NoContext())
        var t = min(_twos(p), k)
        return QQ(Rational(_numerator=p >> t, _denominator=Integer(1) << (k - t)))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _sign_cleared(p: List[Integer], x: QQ, mut budget: Budget) raises PolyError -> Int:
    # The sign of p(a/b) for b > 0: that of sum p_i a^i b^(n-i), by homogeneous
    # Horner on integers, with shifts when b is a power of two.
    var n = len(p) - 1
    budget.charge(UInt64(len(p)))
    var a = x.value.numerator()
    var b = x.value.denominator()
    var e = _power_of_two(b)
    try:
        var acc = p[n]
        var power = Integer(1)
        for i in range(n - 1, -1, -1):
            if e >= 0:
                acc = acc * a + (p[i] << (e * (n - i)))
            else:
                power = power * b
                acc = acc * a + p[i] * power
        return acc.sign()
    except err:
        raise PolyError(PolyErrorKind.CAPACITY, String(err))


def _bisect(mut a: RealAlgebraic, cleared: List[Integer], mut ends: Tuple[Int, Int], mut budget: Budget) raises PolyError:
    # One step of _refine with the endpoint signs remembered: the same
    # interval, but one evaluation at a simple root instead of three, each on
    # the definition with its denominators cleared once.
    if a.rational:
        budget.charge(1)
        return
    if ends[0] == 2 or ends[1] == 2:
        var low = _sign_cleared(cleared, a.low, budget)
        var high = _sign_cleared(cleared, a.high, budget)
        ends = (low, high)
    if ends[0] == 0 or ends[1] == 0 or ends[0] == ends[1]:
        _refine(a, budget)
        ends = (2, 2)
        return
    budget.charge(1)
    var middle = _mid(a.low, a.high, budget)
    var at_middle = _sign_cleared(cleared, middle, budget)
    if at_middle == 0:
        a.low = _mid(a.low, middle, budget)
        a.high = _mid(middle, a.high, budget)
        ends = (2, 2)
    elif at_middle == ends[0]:
        a.low = middle^
        ends = (at_middle, ends[1])
    else:
        a.high = middle^
        ends = (ends[0], at_middle)


def _has_level(mut levels: _PointLevels, r: Int, mut budget: Budget) raises PolyError -> Bool:
    # Build every coordinate's bounds up to level r by bisecting it in place.
    if levels.exact:
        return r == 0
    for i in range(len(levels.coordinates)):
        ref c = levels.coordinates[i]
        while len(levels.bounds[i]) <= r:
            budget.charge(1)
            var last = len(levels.bounds[i]) - 1
            var low = c.low.copy()
            var high = c.high.copy()
            c.low = levels.bounds[i][last][0].copy()
            c.high = levels.bounds[i][last][1].copy()
            _bisect(c, levels.cleared[i], levels.ends[i], budget)
            levels.bounds[i].append((c.low.copy(), c.high.copy()))
            c.low = low^
            c.high = high^
    return True


def _level_box(levels: _PointLevels, r: Int, mut budget: Budget) raises PolyError -> List[Tuple[QQ, QQ]]:
    var box = List[Tuple[QQ, QQ]]()
    for i in range(len(levels.coordinates)):
        budget.charge(1)
        box.append((levels.bounds[i][r][0].copy(), levels.bounds[i][r][1].copy()))
    return box^


def _level_point(levels: _PointLevels, r: Int, mut budget: Budget) raises PolyError -> List[RealAlgebraic]:
    var point = List[RealAlgebraic]()
    for i in range(len(levels.coordinates)):
        var c = _copy_real(levels.coordinates[i], budget)
        c.low = levels.bounds[i][r][0].copy()
        c.high = levels.bounds[i][r][1].copy()
        point.append(c^)
    return point^


def _try_point_interval(f: SparsePoly[QQ], variables: List[Int64], point: List[RealAlgebraic], rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], Int]]:
    var levels = _PointLevels(point, budget)
    return _try_point_levels(f, variables, levels, 0, rounds, budget)


def _try_point_levels(f: SparsePoly[QQ], variables: List[Int64], mut levels: _PointLevels, start: Int, rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], Int]]:
    # The first level r in [start, rounds) whose interval value excludes zero,
    # as testing the levels in order would find it. Interval arithmetic is
    # inclusion-monotone and the levels are nested, so a level that excludes
    # zero is followed only by such levels: gallop, then bisect the gap.
    var g = _integer_image(f, budget)
    var failed = start - 1
    var found = -1
    var relation = 0
    var step = 1
    var r = start
    while r < rounds:
        budget.charge(1)
        if not _has_level(levels, r, budget):
            break
        var value = _interval_signs(g, variables, _level_box(levels, r, budget), budget)
        if value[0] > 0 or value[1] < 0:
            found = r
            relation = 1 if value[0] > 0 else -1
            break
        failed = r
        if r == rounds - 1:
            break
        r = min(r + step, rounds - 1)
        step *= 2
    if found < 0:
        return None
    while found - failed > 1:
        budget.charge(1)
        var middle = failed + (found - failed) // 2
        var value = _interval_signs(g, variables, _level_box(levels, middle, budget), budget)
        if value[0] > 0 or value[1] < 0:
            found = middle
            relation = 1 if value[0] > 0 else -1
        else:
            failed = middle
    return (_level_point(levels, found, budget), relation)


def _try_point_zero(f: SparsePoly[QQ], variables: List[Int64], point: List[RealAlgebraic], rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], Int, List[DensePoly[QQ]]]]:
    # The last irrational coordinate is primitive; at most one other may be irrational.
    var primitive = -1
    var other = -1
    for i in range(len(point)):
        if not point[i].rational:
            if primitive >= 0:
                if other >= 0:
                    return None
                other = primitive
            primitive = i
    if primitive < 0:
        return None
    ref m = point[primitive].polynomial
    var representations = List[DensePoly[QQ]]()
    for i in range(len(point)):
        budget.charge(1)
        if i == primitive:
            representations.append(DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget))
        elif point[i].rational:
            representations.append(DensePoly[QQ].from_coefficients(0, NoContext(), [point[i].low.copy()], budget))
        else:
            representations.append(DensePoly[QQ](0, NoContext()))
    if other >= 0:
        # Rows of f as a polynomial in the other coordinate x over Q(theta): gcd with m_x.
        var assignment = List[Tuple[Int64, RealAlgebraic]]()
        for i in range(len(point)):
            if i != other and i != primitive:
                assignment.append((variables[i], _copy_real(point[i], budget)))
        var degree = Int(f.degree_in(variables[primitive], budget))
        var rows = _coefficient_polynomials(f, variables[primitive], variables[other], assignment, degree + 1, budget)
        var width = 0
        for r in rows:
            width = max(width, len(r.coefficients))
        var specialized = List[DensePoly[QQ]]()
        for k in range(width):
            budget.charge(1)
            var entries = List[QQ]()
            for i in range(len(rows)):
                entries.append(rows[i].coefficients[k].copy() if k < len(rows[i].coefficients) else QQ.zero(NoContext()))
            specialized.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), m, budget))
        var definition = List[DensePoly[QQ]]()
        for c in point[other].polynomial.coefficients:
            definition.append(DensePoly[QQ].from_coefficients(0, NoContext(), [c.copy()], budget))
        var common = _gcd_over_field(definition, specialized, m, budget)
        if len(common) != 2:
            return None
        representations[other] = DensePoly[QQ](0, NoContext()).sub(common[0], budget)
    if not _sparse_compose(f, variables, representations, m, budget).is_zero():
        return None
    var refined = List[RealAlgebraic]()
    for c in point:
        refined.append(_copy_real(c, budget))
    for _ in range(rounds):
        budget.charge(1)
        if other < 0 or _inside(_dense_interval(representations[other], refined[primitive], budget), refined[other], budget):
            return (refined^, primitive, representations^)
        _refine(refined[primitive], budget)
    return None


# ------------------------------------------------------------------ zero proofs through a common divisor

def _divides(divisor: List[DensePoly[QQ]], dividend: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> Bool:
    # Division by a monic polynomial in Q(b)[x] leaves remainder zero.
    var r = _trim(dividend, m, budget)
    var n = len(divisor)
    if len(r) < n or n == 0 or m.degree() < 1 or divisor[n - 1].degree() != 0 or divisor[n - 1].coefficients[0].value != Rational(1):
        return _divides_rational(divisor, r, m, budget)
    return _monic_division(divisor, r, m, False, budget)[0]


def _monic_division(divisor: List[DensePoly[QQ]], r: List[DensePoly[QQ]], m: DensePoly[QQ], quotient: Bool, mut budget: Budget) raises PolyError -> Tuple[Bool, List[DensePoly[QQ]]]:
    # Whether a monic divisor divides r in Q(b)[x], with the quotient when
    # asked; r is reduced and trimmed, at least as long as the divisor, and
    # m has positive degree. Each remainder coefficient is kept as P / s with
    # P an integer polynomial reduced modulo the integer multiple M of m: a
    # coefficient changes at most len(divisor) times, so only the one divided
    # by is put in lowest terms, instead of every coefficient of every product.
    var n = len(divisor)
    var z = _cleared(m, budget)
    ref M = z[0]
    var c = Integer(1)
    for g in divisor:
        for q in g.coefficients:
            budget.charge(1)
            if q.value.denominator() != 1:
                c = _denominator(c, q, budget)
    var G = List[List[Integer]]()
    for g in divisor:
        G.append(_scaled_list(g, c, budget))
    var nums = List[List[Integer]]()
    var dens = List[Integer]()
    for q in r:
        var cleared = _cleared(q, budget)
        nums.append(cleared[0].copy())
        dens.append(cleared[1].copy())
    var result = List[DensePoly[QQ]]()
    if quotient:
        for _ in range(len(r) - n + 1):
            result.append(DensePoly[QQ](0, NoContext()))
    while len(nums) >= n:
        budget.charge(1)
        var top = len(nums) - 1
        var shift = top - (n - 1)
        var t = _lowest_terms(nums[top], dens[top], budget)
        ref T = t[0]
        if quotient:
            var coefficients = List[QQ]()
            budget.charge(UInt64(len(T)))
            try:
                for k in range(len(T)):
                    coefficients.append(QQ(Rational(T[k].copy(), t[1].copy())))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            result[shift] = DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)
        for i in range(n - 1):
            budget.charge(UInt64(_product(len(T), len(G[i]))))
            budget.charge(UInt64(_product(_sum(len(T), len(G[i])), len(M))))
            var j = shift + i
            var reduced = _reduced_product(T, G[i], M, budget)
            budget.charge(UInt64(_sum(len(reduced[0]), len(nums[j]))))
            try:
                # P_j / s_j - R / (s_T c L^e) over the common denominator.
                var scale = t[1] * c * reduced[1]
                var size = max(len(nums[j]), len(reduced[0]))
                var updated = List[Integer](capacity=size)
                for k in range(size):
                    var a = nums[j][k] * scale if k < len(nums[j]) else Integer(0)
                    if k < len(reduced[0]):
                        a = a - reduced[0][k] * dens[j]
                    updated.append(a^)
                while len(updated) and not updated[len(updated) - 1]:
                    _ = updated.pop()
                nums[j] = updated^
                dens[j] = dens[j] * scale
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
        _ = nums.pop()
        _ = dens.pop()
        while len(nums) and len(nums[len(nums) - 1]) == 0:
            _ = nums.pop()
            _ = dens.pop()
    return (len(nums) == 0, result^)


def _divides_rational(divisor: List[DensePoly[QQ]], dividend: List[DensePoly[QQ]], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> Bool:
    var r = _trim(dividend, m, budget)
    var n = len(divisor)
    while len(r) >= n:
        budget.charge(1)
        var factor = r[len(r) - 1].copy()
        var shift = len(r) - n
        for i in range(n):
            r[shift + i] = r[shift + i].sub(_field_mul(factor, divisor[i], m, budget), budget)
        while len(r) and r[len(r) - 1].is_zero():
            _ = r.pop()
    return len(r) == 0


def _scaled_list(p: DensePoly[QQ], common: Integer, mut budget: Budget) raises PolyError -> List[Integer]:
    # The integer coefficients of common p, for common a multiple of every denominator.
    budget.charge(UInt64(len(p.coefficients)))
    var result = List[Integer](capacity=len(p.coefficients))
    try:
        for q in p.coefficients:
            var d = q.value.denominator()
            result.append(q.value.numerator() * common if d == 1 else q.value.numerator() * (common // d))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return result^


def _lowest_terms(p: List[Integer], s: Integer, mut budget: Budget) raises PolyError -> Tuple[List[Integer], Integer]:
    # P / s with the common factor of P's coefficients and s removed.
    budget.charge(UInt64(len(p)))
    try:
        var g = s.copy()
        for a in p:
            if g == 1:
                return (p.copy(), s.copy())
            g = integer_gcd(g, a)
        if g == 1:
            return (p.copy(), s.copy())
        var result = List[Integer](capacity=len(p))
        for a in p:
            result.append(a // g)
        return (result^, s // g)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _reduced_product(A: List[Integer], B: List[Integer], M: List[Integer], mut budget: Budget) raises PolyError -> Tuple[List[Integer], Integer]:
    # (R, L^e) with R = L^e (A B mod M) for the leading coefficient L of M, by
    # integer pseudo-division; the charges are the caller's.
    var dm = len(M) - 1
    try:
        var n = len(A) + len(B) - 1
        var product = List[Integer](length=max(n, 0), fill=Integer(0))
        for i in range(len(A)):
            if not A[i]:
                continue
            for j in range(len(B)):
                product[i + j] = product[i + j] + A[i] * B[j]
        var lead = M[dm]
        var scale = Integer(1)
        if n - 1 >= dm:
            for top in range(n - 1, dm - 1, -1):
                var t = product[top]
                if not t:
                    continue
                if lead != 1:
                    for k in range(top + 1):
                        product[k] = product[k] * lead
                    scale = scale * lead
                for k in range(dm + 1):
                    product[top - dm + k] = product[top - dm + k] - t * M[k]
            n = dm
        var result = List[Integer](capacity=max(n, 0))
        for k in range(n):
            result.append(product[k])
        while len(result) and not result[len(result) - 1]:
            _ = result.pop()
        return (result^, scale^)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _divisor_value(divisor: List[DensePoly[QQ]], x: QQ, b: RealAlgebraic, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    # sum [g_k(b)] x^k over intervals.
    var point = (x.copy(), x.copy())
    var value = (QQ.zero(NoContext()), QQ.zero(NoContext()))
    for k in range(len(divisor) - 1, -1, -1):
        budget.charge(1)
        value = _interval_add(_interval_mul(value, point, budget), _dense_interval(divisor[k], b, budget), budget)
    return value^


def _sign_change(divisor: List[DensePoly[QQ]], a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    var left = _divisor_value(divisor, a.low, b, budget)
    var right = _divisor_value(divisor, a.high, b, budget)
    return (left[0].value.sign() > 0 and right[1].value.sign() < 0) or (left[1].value.sign() < 0 and right[0].value.sign() > 0)


def _rows_over_field(f: SparsePoly[QQ], variables: List[Int64], point: List[RealAlgebraic], primitive: Int, other: Int, m: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    # f(x, b) with every rational coordinate substituted, as coefficients in Q(b) by powers of x.
    var assignment = List[Tuple[Int64, RealAlgebraic]]()
    for i in range(len(point)):
        if i != other and i != primitive:
            if not point[i].rational:
                raise PolyError(PolyErrorKind.DOMAIN, "a third irrational coordinate")
            assignment.append((variables[i], _copy_real(point[i], budget)))
    var degree = Int(f.degree_in(variables[primitive], budget))
    var rows = _coefficient_polynomials(f, variables[primitive], variables[other], assignment, degree + 1, budget)
    var width = 0
    for r in rows:
        width = max(width, len(r.coefficients))
    var specialized = List[DensePoly[QQ]]()
    for k in range(width):
        budget.charge(1)
        var entries = List[QQ]()
        for i in range(len(rows)):
            entries.append(rows[i].coefficients[k].copy() if k < len(rows[i].coefficients) else QQ.zero(NoContext()))
        specialized.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), m, budget))
    return specialized^


def _constants(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    var result = List[DensePoly[QQ]]()
    for c in p.coefficients:
        result.append(DensePoly[QQ].from_coefficients(0, NoContext(), [c.copy()], budget))
    return result^


def _check_point_divisor(f: SparsePoly[QQ], variables: List[Int64], point: List[RealAlgebraic], refined: List[RealAlgebraic], primitive: Int, divisor: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    var k = len(point)
    if len(refined) != k or primitive < 0 or primitive >= k or point[primitive].rational or len(divisor) < 2:
        return False
    var other = -1
    for i in range(k):
        budget.charge(1)
        if not _same_root(refined[i], point[i], budget):
            return False
        if i != primitive and not point[i].rational:
            if other >= 0:
                return False
            other = i
    if other < 0:
        return False
    ref m = point[primitive].polynomial
    for g in divisor:
        g.check_invariants(budget)
        if g.variable != 0 or g.degree() >= m.degree():
            return False
    var one = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], budget)
    if not _dense_equal(divisor[len(divisor) - 1], one, budget):
        return False
    # The divisor's roots are roots of the other coordinate's definition and of f(x, b);
    # one of them lies in that coordinate's isolating interval, so it is the coordinate.
    if not _divides(divisor, _constants(point[other].polynomial, budget), m, budget):
        return False
    if not _divides(divisor, _rows_over_field(f, variables, point, primitive, other, m, budget), m, budget):
        return False
    return _sign_change(divisor, refined[other], refined[primitive], budget)


def _try_point_divisor(f: SparsePoly[QQ], variables: List[Int64], point: List[RealAlgebraic], rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], Int, List[DensePoly[QQ]]]]:
    var primitive = -1
    var other = -1
    for i in range(len(point)):
        if not point[i].rational:
            if primitive >= 0:
                if other >= 0:
                    return None
                other = primitive
            primitive = i
    if primitive < 0 or other < 0:
        return None
    ref m = point[primitive].polynomial
    var specialized = _rows_over_field(f, variables, point, primitive, other, m, budget)
    var divisor = _gcd_over_field(_constants(point[other].polynomial, budget), specialized, m, budget)
    if len(divisor) < 2:
        return None
    var refined = List[RealAlgebraic]()
    for c in point:
        refined.append(_copy_real(c, budget))
    for round in range(rounds):
        budget.charge(1)
        if _sign_change(divisor, refined[other], refined[primitive], budget):
            return (refined^, primitive, divisor^)
        # The divisor's coefficients carry the primitive coordinate's interval
        # width; it must shrink faster than the other coordinate's interval.
        _refine(refined[primitive], budget)
        if round % 3 == 0:
            _refine(refined[other], budget)
    return None


# ------------------------------------------------------------------ root filters through a common divisor

def _check_root_divisor(coefficients: List[RealAlgebraic], candidate: RealAlgebraic, refined: List[RealAlgebraic], refined_candidate: RealAlgebraic, element: RealAlgebraic, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    # c_i = r_i(e) for an algebraic e, and a monic G over Q(b) divides m_e(x) and
    # sum_i r_i(x) b^i and changes sign on e's interval: then sum_i c_i b^i = 0.
    var n = len(coefficients)
    if len(refined) != n or len(representations) < n + 2 or element.rational or not _check_real(element, budget) or not _same_root(refined_candidate, candidate, budget) or candidate.rational:
        return False
    ref me = element.polynomial
    ref mb = candidate.polynomial
    for i in range(n):
        budget.charge(1)
        ref r = representations[i]
        r.check_invariants(budget)
        if r.variable != 0 or not _same_root(refined[i], coefficients[i], budget):
            return False
        if coefficients[i].rational:
            var constant = DensePoly[QQ].from_coefficients(0, NoContext(), [coefficients[i].low.copy()], budget)
            if not _reduce(r.sub(constant, budget), me, budget).is_zero():
                return False
        elif not _compose(coefficients[i].polynomial, r, me, budget).is_zero() or not _inside(_dense_interval(r, element, budget), refined[i], budget):
            return False
    var divisor = List[DensePoly[QQ]]()
    for i in range(n, len(representations)):
        representations[i].check_invariants(budget)
        if representations[i].variable != 0:
            return False
        divisor.append(_reduce(representations[i], mb, budget))
    var one = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], budget)
    if not _dense_equal(divisor[len(divisor) - 1], one, budget):
        return False
    # sum_i r_i(x) b^i, as coefficients in Q(b) by powers of x.
    var width = 0
    for i in range(n):
        width = max(width, len(representations[i].coefficients))
    var rows = List[DensePoly[QQ]]()
    for k in range(width):
        budget.charge(1)
        var entries = List[QQ]()
        for i in range(n):
            entries.append(representations[i].coefficients[k].copy() if k < len(representations[i].coefficients) else QQ.zero(NoContext()))
        rows.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), mb, budget))
    if not _divides(divisor, _constants(me, budget), mb, budget) or not _divides(divisor, rows, mb, budget):
        return False
    return _sign_change(divisor, element, refined_candidate, budget)


def _try_root_divisor(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], coefficients: List[RealAlgebraic], candidate: RealAlgebraic, rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], RealAlgebraic, RealAlgebraic, List[DensePoly[QQ]]]]:
    if candidate.rational:
        return None
    var irrational = -1
    for i in range(len(point)):
        if not point[i][1].rational:
            if irrational >= 0:
                return None
            irrational = i
    if irrational < 0:
        return None
    var element = _copy_real(point[irrational][1], budget)
    ref me = point[irrational][1].polynomial
    var rows = _coefficient_polynomials(f, variable, point[irrational][0], point, len(coefficients), budget)
    var representations = List[DensePoly[QQ]]()
    for r in rows:
        representations.append(_reduce(r, me, budget))
    return _element_divisor(element^, representations^, coefficients, candidate, rounds, budget)


def _element_divisor(var element: RealAlgebraic, var representations: List[DensePoly[QQ]], coefficients: List[RealAlgebraic], candidate: RealAlgebraic, rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[List[RealAlgebraic], RealAlgebraic, RealAlgebraic, List[DensePoly[QQ]]]]:
    # A kind-3 witness from c_i = r_i(element): a common divisor of m_e(x) and
    # sum_i r_i(x) b^i over Q(b) that changes sign on the element's interval.
    if candidate.rational or element.rational:
        return None
    ref me = element.polynomial
    ref mb = candidate.polynomial
    var width = 0
    for r in representations:
        width = max(width, len(r.coefficients))
    var specialized = List[DensePoly[QQ]]()
    for k in range(width):
        budget.charge(1)
        var entries = List[QQ]()
        for i in range(len(representations)):
            entries.append(representations[i].coefficients[k].copy() if k < len(representations[i].coefficients) else QQ.zero(NoContext()))
        specialized.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), mb, budget))
    var divisor = _gcd_over_field(_constants(me, budget), specialized, mb, budget)
    if len(divisor) < 2:
        return None
    var refined = List[RealAlgebraic]()
    for c in coefficients:
        refined.append(_copy_real(c, budget))
    var point_b = _copy_real(candidate, budget)
    for round in range(rounds):
        budget.charge(1)
        var ready = _sign_change(divisor, element, point_b, budget)
        if ready:
            for i in range(len(coefficients)):
                if not coefficients[i].rational and not _inside(_dense_interval(representations[i], element, budget), refined[i], budget):
                    ready = False
                    break
        if ready:
            for g in divisor:
                representations.append(g.copy())
            return (refined^, point_b^, element^, representations^)
        _refine(point_b, budget)
        if round % 3 == 0:
            _refine(element, budget)
    return None


# ------------------------------------------------------------------ direct evaluation at one irrational coordinate

def _single_irrational(f: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> Int:
    # The index of the one irrational coordinate that f uses, -1 if none, -2 if several or a variable is unassigned.
    var found = -1
    for v in f.variables(budget):
        budget.charge(1)
        var index = -1
        for i in range(len(point)):
            if point[i][0] == v:
                index = i
                break
        if index < 0:
            return -2
        if not point[index][1].rational:
            if found >= 0 and found != index:
                return -2
            found = index
    return found


def _image(f: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], index: Int, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # f with every rational coordinate substituted, as a dense polynomial in the irrational one.
    var coordinate = point[index][0]
    var row = List[QQ]()
    for term in f.terms:
        budget.charge(1)
        var degree = 0
        var value = term.coefficient.copy()
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] == coordinate:
                degree = Int(pair[1])
                continue
            for entry in point:
                if entry[0] == pair[0]:
                    for _ in range(Int(pair[1])):
                        value = value.mul(entry[1].low, budget)
                    break
        while len(row) <= degree:
            row.append(QQ.zero(NoContext()))
        row[degree] = row[degree].add(value, budget)
    return DensePoly[QQ].from_coefficients(0, NoContext(), row^, budget)


def _check_direct(f: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], value: RealAlgebraic, witness: List[RealAlgebraic], mut budget: Budget) raises PolyError -> Bool:
    var index = _single_irrational(f, point, budget)
    if index < 0 or len(witness) != 1 or not _same_root(witness[0], point[index][1], budget) or not _check_real(value, budget):
        return False
    ref m = point[index][1].polynomial
    var r = _reduce(_image(f, point, index, budget), m, budget)
    # value is a root of its definition at r(alpha); for an irrational value the
    # image interval selects it among that definition's roots.
    if not _vanishes_at(_compose(value.polynomial, r, m, budget), witness[0], budget):
        return False
    return value.rational or _inside(_dense_interval(r, witness[0], budget), value, budget)


def _vanishes_at(c: DensePoly[QQ], alpha: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    """Whether c, reduced modulo alpha's polynomial m, vanishes at alpha.

    c = 0 suffices. Otherwise m may be reducible, so c need not be zero modulo
    m: a divisor g of both m and c, checked by exact division rather than
    trusted from the gcd, that changes sign across alpha's interval has a root
    there, which as a root of m can only be alpha, so c(alpha) = 0.
    """
    if c.is_zero():
        return True
    ref m = alpha.polynomial
    var g = _qq_gcd(m, c, budget)
    if g.degree() < 1 or not _remainder(m, g, budget).is_zero() or not _remainder(c, g, budget).is_zero():
        return False
    var left = _sign_at(g, alpha.low, budget)
    return left != 0 and _sign_at(g, alpha.high, budget) == -left


def _direct_value(f: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> Optional[Tuple[RealAlgebraic, RealAlgebraic]]:
    """The value of f at a point with one irrational coordinate alpha, as r(alpha) for r = f mod m_alpha."""
    var index = _single_irrational(f, point, budget)
    if index < 0:
        return None
    ref alpha = point[index][1]
    ref m = alpha.polynomial
    var r = _reduce(_image(f, point, index, budget), m, budget)
    var witness = _copy_real(alpha, budget)
    if r.degree() <= 0:
        var q = r.coefficients[0].copy() if r.degree() == 0 else QQ.zero(NoContext())
        return (_rational(q, budget), witness^)
    # The minimal polynomial of r(alpha) divides Res_x(m(x), z - r(x)).
    var h = _norm_image(r, m, budget)
    if h.is_zero():
        return None
    for _ in range(4096):
        budget.charge(1)
        var image = _dense_interval(r, witness, budget)
        if _compare_q(image[0], image[1], budget) < 0 and _sign_at(h, image[0], budget) != 0 and _sign_at(h, image[1], budget) != 0 and _unique_root(h, image[0], image[1], budget):
            var value = _root_cached(h, image[0], image[1], cache, budget)
            # The checker needs an irrational value's interval to hold the
            # image, which normalization may have narrowed. The value lies
            # strictly inside it and the image shrinks to the value, so
            # refining the witness alone gets there.
            while not value.rational and not _inside(image, value, budget):
                _refine(witness, budget)
                image = _dense_interval(r, witness, budget)
            return (value^, witness^)
        _refine(witness, budget)
    return None


def _norm_image(r: DensePoly[QQ], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # Res_x(m(x), z - r(x)) = lc(m)^deg(r) * charpoly of r(alpha); deg(r) >= 1 and deg(r) < deg(m).
    var d = m.degree()
    var lead = m.coefficients[d].copy()
    var elementary = _conjugate_symmetric(r, m, budget)
    var scale = QQ.one(NoContext())
    for _ in range(r.degree()):
        budget.charge(1)
        scale = scale.mul(lead, budget)
    var coefficients = List[QQ](capacity=_sum(d, 1))
    for i in range(d + 1):
        budget.charge(1)
        var value = elementary[d - i].mul(scale, budget)
        coefficients.append(value.neg(budget) if (d - i) % 2 == 1 else value^)
    return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)


def _conjugate_symmetric(r: DensePoly[QQ], m: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[QQ]:
    # The elementary symmetric functions e_0..e_d of r(alpha_i) over the roots
    # alpha_i of m, from the traces of r^k by Newton's identities; r reduced.
    var d = m.degree()
    var lead = m.coefficients[d].copy()
    var inverse = lead.inverse(budget)
    _allocation[QQ](_sum(d, 1), budget)
    var monic = List[QQ](capacity=d)
    for j in range(d):
        budget.charge(1)
        monic.append(m.coefficients[j].mul(inverse, budget))
    # Power sums of the roots of m: p_k = -(k c_{d-k} + sum_{i<k} c_{d-i} p_{k-i}).
    var sums = List[QQ](capacity=d)
    sums.append(QQ(Rational(Integer(d))))
    for k in range(1, d):
        budget.charge(1)
        var total = monic[d - k].mul(QQ(Rational(Integer(k))), budget)
        for i in range(1, k):
            budget.charge(1)
            total = total.add(monic[d - i].mul(sums[k - i], budget), budget)
        sums.append(total.neg(budget))
    # Elementary symmetric functions of the conjugates of r(alpha).
    var traces = List[QQ](capacity=_sum(d, 1))
    traces.append(QQ.zero(NoContext()))
    var power = _copy_poly(r, budget)
    for k in range(1, d + 1):
        budget.charge(1)
        if k > 1:
            power = _field_mul(power, r, m, budget)
        var trace = QQ.zero(NoContext())
        for j in range(len(power.coefficients)):
            budget.charge(1)
            trace = trace.add(power.coefficients[j].mul(sums[j], budget), budget)
        traces.append(trace^)
    var elementary = List[QQ](capacity=_sum(d, 1))
    elementary.append(QQ.one(NoContext()))
    for k in range(1, d + 1):
        budget.charge(1)
        var total = QQ.zero(NoContext())
        for i in range(1, k + 1):
            budget.charge(1)
            var term = elementary[k - i].mul(traces[i], budget)
            total = total.add(term, budget) if i % 2 == 1 else total.sub(term, budget)
        elementary.append(total.mul(QQ(Rational(Integer(k))).inverse(budget), budget))
    return elementary^


def _embed_dense(p: DensePoly[QQ], variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    var moved = _copy_poly(p, budget)
    moved.variable = variable
    return moved.to_sparse(MonomialOrder(), budget)


# ------------------------------------------------------------------ tower elimination

def _check_coefficient_images(coefficients: List[RealAlgebraic], element: RealAlgebraic, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    # c_i = r_i(e) for every coefficient.
    if len(representations) != len(coefficients) or element.rational or not _check_real(element, budget):
        return False
    ref me = element.polynomial
    for i in range(len(coefficients)):
        budget.charge(1)
        ref r = representations[i]
        r.check_invariants(budget)
        if r.variable != 0:
            return False
        if coefficients[i].rational:
            var constant = DensePoly[QQ].from_coefficients(0, NoContext(), [coefficients[i].low.copy()], budget)
            if not _reduce(r.sub(constant, budget), me, budget).is_zero():
                return False
        elif not _compose(coefficients[i].polynomial, r, me, budget).is_zero() or not _inside(_dense_interval(r, element, budget), coefficients[i], budget):
            return False
    return True


def _tower_definitions(element: RealAlgebraic, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Tuple[SparsePoly[QQ], SparsePoly[QQ]]:
    # sum_i r_i(x) y^i with y variable zero and x variable one, and m_e(x).
    var order = MonomialOrder()
    var y = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var total = SparsePoly[QQ].zero(order, NoContext(), budget)
    for i in range(len(representations) - 1, -1, -1):
        budget.charge(1)
        total = total.mul(y, budget).add(_embed_dense(representations[i], 1, budget), budget)
    return (total^, _embed_dense(element.polynomial, 1, budget))


def _try_tower(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], coefficients: List[RealAlgebraic], rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[RealAlgebraic, List[DensePoly[QQ]]]]:
    var irrational = -1
    for i in range(len(point)):
        if not point[i][1].rational:
            if irrational >= 0:
                return None
            irrational = i
    if irrational < 0:
        return None
    var element = _copy_real(point[irrational][1], budget)
    ref me = point[irrational][1].polynomial
    var rows = _coefficient_polynomials(f, variable, point[irrational][0], point, len(coefficients), budget)
    var representations = List[DensePoly[QQ]]()
    for r in rows:
        representations.append(_reduce(r, me, budget))
    for _ in range(rounds):
        budget.charge(1)
        if _check_coefficient_images(coefficients, element, representations, budget):
            return (element^, representations^)
        _refine(element, budget)
    return None


# ------------------------------------------------------------------ frames: every coordinate in one generator

def _constant(value: QQ, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    return DensePoly[QQ].from_coefficients(0, NoContext(), [value.copy()], budget)


def _frame_applies(frame: _Frame, point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> Bool:
    # The frame names the point's coordinates in order, possibly with more after them.
    if not frame.valid or not frame.generator or len(point) > len(frame.variables):
        return False
    for j in range(len(point)):
        budget.charge(1)
        if frame.variables[j] != point[j][0]:
            return False
    return True


def _frame_rows(f: SparsePoly[QQ], variable: Int64, variables: List[Int64], representations: List[DensePoly[QQ]], m: DensePoly[QQ], size: Int, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    # The coefficient of variable^k in f(r_1(z), ..., r_n(z), variable) modulo m(z).
    _allocation[DensePoly[QQ]](size, budget)
    var rows = List[DensePoly[QQ]](capacity=size)
    for _ in range(size):
        rows.append(DensePoly[QQ](0, NoContext()))
    for term in f.terms:
        budget.charge(1)
        var power = 0
        var value = _constant(term.coefficient, budget)
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] == variable:
                power = Int(pair[1])
                continue
            var index = -1
            for i in range(len(variables)):
                if variables[i] == pair[0]:
                    index = i
                    break
            if index < 0:
                raise PolyError(PolyErrorKind.DOMAIN, "missing coordinate")
            for _ in range(Int(pair[1])):
                value = _field_mul(value, representations[index], m, budget)
        if power >= size:
            raise PolyError(PolyErrorKind.DOMAIN, "coefficient degree mismatch")
        rows[power] = rows[power].add(value, budget)
    for i in range(size):
        rows[i] = _reduce(rows[i], m, budget)
    return rows^


def _used_irrational(f: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> List[Int]:
    # Indices of the irrational coordinates whose variables occur in f.
    var used = f.variables(budget)
    var result = List[Int]()
    for j in range(len(point)):
        budget.charge(1)
        if point[j][1].rational:
            continue
        for v in used:
            if v == point[j][0]:
                result.append(j)
                break
    return result^


def _substitutions(point: List[Tuple[Int64, RealAlgebraic]], representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    # Rational coordinates by their values, irrational ones by their representations.
    var result = List[DensePoly[QQ]]()
    for j in range(len(point)):
        budget.charge(1)
        if point[j][1].rational:
            result.append(_constant(point[j][1].low, budget))
        else:
            result.append(_copy_poly(representations[j], budget))
    return result^


def _represents(point: List[Tuple[Int64, RealAlgebraic]], used: List[Int], generator: RealAlgebraic, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    # Each used irrational coordinate is the root of its definition that r_j(generator) selects.
    ref m = generator.polynomial
    for j in used:
        budget.charge(1)
        ref r = representations[j]
        r.check_invariants(budget)
        if r.variable != 0 or not _compose(point[j][1].polynomial, r, m, budget).is_zero() or not _inside(_dense_interval(r, generator, budget), point[j][1], budget):
            return False
    return True


def _contained(point: List[Tuple[Int64, RealAlgebraic]], used: List[Int], generator: RealAlgebraic, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    # The interval half of _represents; frames hold exact representations by construction.
    for j in used:
        budget.charge(1)
        if not _inside(_dense_interval(representations[j], generator, budget), point[j][1], budget):
            return False
    return True


def _settle(mut frame: _Frame, point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError:
    # Refine the generator once until every irrational coordinate is contained.
    var used = List[Int]()
    for j in range(len(point)):
        if not point[j][1].rational:
            used.append(j)
    for _ in range(4096):
        budget.charge(1)
        if _contained(point, used, frame.generator.value(), frame.representations, budget):
            return
        _refine(frame.generator.value(), budget)
    frame.valid = False


def _rational_frame(point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> _Frame:
    """A frame for a point whose coordinates are all rational: the generator is
    the root 0 of t on (-1, 1) and every coordinate a constant, so a direct
    evaluation is f at the point over Q."""
    var t = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
    var frame = _Frame()
    frame.generator = RealAlgebraic(t^, QQ.one(NoContext()).neg(budget), QQ.one(NoContext()), False)
    for entry in point:
        budget.charge(1)
        frame.variables.append(entry[0])
        frame.representations.append(_constant(entry[1].low, budget))
    return frame^


def _installed(frame: _Frame, point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> Optional[_Frame]:
    """A settled copy of a usable frame for the point, for a lifting cache."""
    if not frame.valid or not frame.generator:
        return None
    var settled = frame.copy()
    _settle(settled, point, budget)
    if not settled.valid:
        return None
    return settled^


def _check_frame_direct(f: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], value: RealAlgebraic, witness: List[RealAlgebraic], representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    if len(witness) != 1 or len(representations) != len(point) or witness[0].rational or not _check_real(witness[0], budget) or not _check_real(value, budget):
        return False
    var used = _used_irrational(f, point, budget)
    if not _represents(point, used, witness[0], representations, budget):
        return False
    var variables = List[Int64]()
    for entry in point:
        variables.append(entry[0])
    ref m = witness[0].polynomial
    var r = _sparse_compose(f, variables, _substitutions(point, representations, budget), m, budget)
    if not _vanishes_at(_compose(value.polynomial, r, m, budget), witness[0], budget):
        return False
    return value.rational or _inside(_dense_interval(r, witness[0], budget), value, budget)


def _frame_value(f: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], frame: _Frame, mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> Optional[Tuple[RealAlgebraic, RealAlgebraic, List[DensePoly[QQ]]]]:
    """The value of f at a point whose coordinates are polynomials in the frame's generator."""
    if not _frame_applies(frame, point, budget):
        return None
    ref m = frame.generator.value().polynomial
    var representations = List[DensePoly[QQ]]()
    var variables = List[Int64]()
    for j in range(len(point)):
        representations.append(_copy_poly(frame.representations[j], budget))
        variables.append(point[j][0])
    var used = _used_irrational(f, point, budget)
    var r = _sparse_compose(f, variables, _substitutions(point, representations, budget), m, budget)
    var witness = _copy_real(frame.generator.value(), budget)
    var h = DensePoly[QQ](0, NoContext())
    if r.degree() > 0:
        h = _norm_image(r, m, budget)
    for _ in range(4096):
        budget.charge(1)
        if _contained(point, used, witness, representations, budget):
            if r.degree() <= 0:
                var q = r.coefficients[0].copy() if r.degree() == 0 else QQ.zero(NoContext())
                return (_rational(q, budget), witness^, representations^)
            var image = _dense_interval(r, witness, budget)
            if _compare_q(image[0], image[1], budget) < 0 and _sign_at(h, image[0], budget) != 0 and _sign_at(h, image[1], budget) != 0 and _unique_root(h, image[0], image[1], budget):
                var value = _root_cached(h, image[0], image[1], cache, budget)
                while not value.rational and not _inside(image, value, budget):
                    _refine(witness, budget)
                    image = _dense_interval(r, witness, budget)
                return (value^, witness^, representations^)
        _refine(witness, budget)
    return None


def _try_frame_tower(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], coefficients: List[RealAlgebraic], frame: _Frame, rounds: Int, mut budget: Budget) raises PolyError -> Optional[Tuple[RealAlgebraic, List[DensePoly[QQ]]]]:
    # As _try_tower, with every coordinate a polynomial in the frame's generator.
    if not _frame_applies(frame, point, budget):
        return None
    var element = _copy_real(frame.generator.value(), budget)
    var variables = List[Int64]()
    for entry in point:
        variables.append(entry[0])
    var representations = _frame_rows(f, variable, variables, _substitutions(point, frame.representations, budget), element.polynomial, len(coefficients), budget)
    for _ in range(rounds):
        budget.charge(1)
        if _check_coefficient_images(coefficients, element, representations, budget):
            return (element^, representations^)
        _refine(element, budget)
    return None


comptime _TERMINAL_COMBINATION = 16
"""Degree product of generator and sample above which a terminal cell combines them."""


def _child_frame(parent: _Frame, variable: Int64, lifting: Optional[SparsePoly[QQ]], sample: RealAlgebraic, mut budget: Budget, terminal: Bool = False) raises PolyError -> _Frame:
    """The frame of a point extended by one sample, invalid when the sample cannot carry it.

    A terminal cell whose sample does not generate the field combines it with
    the generator only above `_TERMINAL_COMBINATION`: below, its few
    evaluations cost less as arithmetic transcripts than the combination.
    """
    var child = parent.copy()
    if not parent.valid:
        return child^
    child.variables.append(variable)
    if sample.rational:
        child.representations.append(_constant(sample.low, budget))
        return child^
    var z = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
    ref mb = sample.polynomial
    if not parent.generator:
        # Every earlier coordinate is rational: the sample generates the point.
        child.representations.append(z^)
        child.generator = _copy_real(sample, budget)
        return child^
    if not lifting:
        child.valid = False
        return child^
    ref gamma = parent.generator.value()
    ref mg = gamma.polynomial
    # The sample is a root of lifting(r(x), y) at x = gamma; over Q(sample) the
    # gcd of m_gamma(x) and lifting(r(x), sample) is x - c(sample) when the
    # sample determines gamma.
    var size = Int(lifting.value().degree_in(variable, budget)) + 1
    var rows = _frame_rows(lifting.value(), variable, parent.variables, parent.representations, mg, size, budget)
    var width = 0
    for r in rows:
        width = max(width, len(r.coefficients))
    var specialized = List[DensePoly[QQ]]()
    for k in range(width):
        budget.charge(1)
        var entries = List[QQ]()
        for i in range(len(rows)):
            entries.append(rows[i].coefficients[k].copy() if k < len(rows[i].coefficients) else QQ.zero(NoContext()))
        specialized.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), mb, budget))
    var definition = List[DensePoly[QQ]]()
    for c in mg.coefficients:
        definition.append(_constant(c, budget))
    var common = _gcd_over_field(definition, specialized, mb, budget)
    if len(common) != 2:
        # The sample alone does not generate Q(gamma, sample): try sample + t gamma.
        var combine = not terminal or _product(mg.degree(), mb.degree()) > _TERMINAL_COMBINATION
        for t in range(1, 4 if combine else 1):
            var combined = _combined_frame(parent, rows, sample, t, budget)
            if combined:
                var result = combined.take()
                result.variables.append(variable)
                return result^
        child.valid = False
        return child^
    var c = DensePoly[QQ](0, NoContext()).sub(common[0], budget)
    var refined = _copy_real(sample, budget)
    var located = False
    for _ in range(4096):
        budget.charge(1)
        if _inside(_dense_interval(c, refined, budget), gamma, budget):
            located = True
            break
        _refine(refined, budget)
    if not located:
        child.valid = False
        return child^
    for j in range(len(parent.representations)):
        budget.charge(1)
        child.representations[j] = _compose(parent.representations[j], c, mb, budget)
    child.representations.append(z^)
    child.generator = _copy_real(sample, budget)
    return child^


def _interpolate(values: List[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # The polynomial of degree < len(values) taking values[i] at i, by Newton's divided differences.
    var n = len(values)
    var differences = values.copy()
    for level in range(1, n):
        for i in range(n - 1, level - 1, -1):
            budget.charge(1)
            differences[i] = differences[i].sub(differences[i - 1], budget).mul(QQ(Rational(Integer(level))).inverse(budget), budget)
    # Horner in the Newton basis: c_0 + (z - 0)(c_1 + (z - 1)(c_2 + ...)).
    var result = DensePoly[QQ](0, NoContext())
    for i in range(n - 1, -1, -1):
        budget.charge(1)
        var factor = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(Integer(-i))), QQ.one(NoContext())], budget)
        result = result.mul(factor, budget).add(_constant(differences[i], budget), budget)
    return result^


def _combination_norm(rows: List[DensePoly[QQ]], mg: DensePoly[QQ], t: Int, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    """prod_i P(gamma_i, z) over the roots gamma_i of the monic mg, for
    P(x, z) = sum_k rows_k(x) (z - t x)^k: the resultant Res_x(mg, P), taken
    modularly from the integer images. With Pc = cP P and Mc = cM mg,
    Res_x(Pc, Mc) = (-1)^(d e) cM^e cP^d prod_i P(gamma_i, z), d = deg mg and
    e = deg_x P."""
    var order = MonomialOrder()
    var z = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var x = SparsePoly[QQ].variable(1, order, NoContext(), budget)
    var shift = z.sub(x.scale(QQ(Rational(Integer(t))), budget), budget)
    var p = SparsePoly[QQ].zero(order, NoContext(), budget)
    for k in range(len(rows) - 1, -1, -1):
        budget.charge(1)
        p = p.mul(shift, budget).add(_embed_dense(rows[k], 1, budget), budget)
    var m = _embed_dense(mg, 1, budget)
    var cp = Integer(1)
    for term in p.terms:
        cp = _denominator(cp, term.coefficient, budget)
    var cm = Integer(1)
    for term in m.terms:
        cm = _denominator(cm, term.coefficient, budget)
    var d = mg.degree()
    var e = p.degree_in(1, budget)
    var value = resultant_value(_integer_image(p, budget), _integer_image(m, budget), 1, budget)
    var divisor = QQ.one(NoContext())
    for _ in range(d):
        budget.charge(1)
        divisor = divisor.mul(QQ(Rational(cp)), budget)
    for _ in range(e):
        budget.charge(1)
        divisor = divisor.mul(QQ(Rational(cm)), budget)
    if (d * e) % 2 == 1:
        divisor = divisor.neg(budget)
    var norm = _rational_image(value, budget).scale(divisor.inverse(budget), budget)
    return DensePoly[QQ].from_sparse(norm, budget, variable=Int64(0))


def _combined_frame(parent: _Frame, rows: List[DensePoly[QQ]], sample: RealAlgebraic, t: Int, mut budget: Budget) raises PolyError -> Optional[_Frame]:
    """A frame generated by theta = sample + t gamma, where sample is a root of sum_k rows_k(gamma) y^k."""
    ref gamma = parent.generator.value()
    ref mg = gamma.polynomial
    var top = len(rows) - 1
    # N(z) = prod_i P(gamma_i, z) for P(x, z) = sum_k rows_k(x) (z - t x)^k.
    var norm = _combination_norm(rows, mg, t, budget)
    if norm.degree() < 1:
        return None
    # theta lies in sample + t gamma; refine both until one root of N remains there.
    var cache = _NormalizationCache()
    var b = _copy_real(sample, budget)
    var g = _copy_real(gamma, budget)
    var scale = QQ(Rational(Integer(t)))
    var theta = Optional[RealAlgebraic](None)
    for _ in range(4096):
        budget.charge(1)
        var low = b.low.add(g.low.mul(scale, budget), budget)
        var high = b.high.add(g.high.mul(scale, budget), budget)
        if _compare_q(low, high, budget) < 0 and _sign_at(norm, low, budget) != 0 and _sign_at(norm, high, budget) != 0 and _unique_root(norm, low, high, budget):
            var found = _root_cached(norm, low, high, cache, budget)
            while not found.rational and not _inside((low.copy(), high.copy()), found, budget):
                _refine(b, budget)
                _refine(g, budget)
                low = b.low.add(g.low.mul(scale, budget), budget)
                high = b.high.add(g.high.mul(scale, budget), budget)
            theta = found^
            break
        _refine(b, budget)
        _refine(g, budget)
    if not theta:
        return None
    ref mt = theta.value().polynomial
    if mt.degree() < 1:
        return None
    # P(x, theta) over Q(theta), by Horner in (theta - t x).
    var z = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
    var specialized = List[DensePoly[QQ]]()
    for k in range(top, -1, -1):
        budget.charge(1)
        var width = max(_sum(len(specialized), 1), len(rows[k].coefficients))
        var next = List[DensePoly[QQ]]()
        for i in range(width):
            budget.charge(1)
            var entry = DensePoly[QQ](0, NoContext())
            if i < len(specialized):
                entry = _field_mul(specialized[i], z, mt, budget)
            if i >= 1 and i - 1 < len(specialized):
                entry = entry.sub(specialized[i - 1].mul(_constant(scale, budget), budget), budget)
            if i < len(rows[k].coefficients):
                entry = entry.add(_constant(rows[k].coefficients[i], budget), budget)
            next.append(_reduce(entry, mt, budget))
        specialized = next^
    var definition = List[DensePoly[QQ]]()
    for c in mg.coefficients:
        definition.append(_constant(c, budget))
    var common = _gcd_over_field(definition, specialized, mt, budget)
    if len(common) != 2:
        return None
    var c = DensePoly[QQ](0, NoContext()).sub(common[0], budget)
    var beta = _reduce(z.sub(c.mul(_constant(scale, budget), budget), budget), mt, budget)
    var refined = _copy_real(theta.value(), budget)
    var located = False
    for _ in range(4096):
        budget.charge(1)
        if _inside(_dense_interval(c, refined, budget), gamma, budget) and _inside(_dense_interval(beta, refined, budget), sample, budget):
            located = True
            break
        _refine(refined, budget)
    if not located:
        return None
    var frame = parent.copy()
    for j in range(len(parent.representations)):
        budget.charge(1)
        frame.representations[j] = _compose(parent.representations[j], c, mt, budget)
    frame.representations.append(beta^)
    frame.generator = theta.take()
    return frame^


# ------------------------------------------------------------------ tower resultants by norms

def _cleared_product(a: SparsePoly[QQ], b: SparsePoly[QQ], mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    # a b, multiplied over the integers after clearing denominators and
    # rescaled once, instead of normalizing after every rational operation.
    var ca = Integer(1)
    for term in a.terms:
        ca = _denominator(ca, term.coefficient, budget)
    var cb = Integer(1)
    for term in b.terms:
        cb = _denominator(cb, term.coefficient, budget)
    var product = _integer_image(a, budget).mul(_integer_image(b, budget), budget)
    var scale = QQ(Rational(ca)).mul(QQ(Rational(cb)), budget).inverse(budget)
    return _rational_image(product, budget).scale(scale, budget)


def _tower_resultant(element: RealAlgebraic, representations: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Optional[ResultantResult[QQ]]:
    """Res_x(F, m_e) for F = sum_i r_i(x) y^i, with multipliers, as `resultant` of `_tower_definitions`.

    The resultant comes from the integer images of both definitions, whose
    resultant `resultant_value` interpolates modulo primes. It is a multiple of
    the norm of F(e, y), so the multiplier of F is its quotient by F(e, y) over
    Q(e), and that of m_e the exact quotient of the remainder; None when an
    inverse modulo m_e fails or a division is inexact.
    """
    ref m = element.polynomial
    var d = m.degree()
    var top = len(representations) - 1
    while top >= 0 and representations[top].is_zero():
        top -= 1
    if top < 1 or d < 1:
        return None
    var nx = 0
    for r in representations:
        nx = max(nx, r.degree())
    # Res(c_F F, c_m m) = c_F^d c_m^deg_x(F) Res(F, m) over the integers, where
    # resultant_value evaluates and interpolates modulo primes.
    var definitions = _tower_definitions(element, representations, budget)
    var cf = Integer(1)
    for term in definitions[0].terms:
        cf = _denominator(cf, term.coefficient, budget)
    var cm = Integer(1)
    for term in definitions[1].terms:
        cm = _denominator(cm, term.coefficient, budget)
    var integer = resultant_value(_integer_image(definitions[0], budget), _integer_image(definitions[1], budget), 1, budget)
    var divisor = QQ.one(NoContext())
    for _ in range(d):
        budget.charge(1)
        divisor = divisor.mul(QQ(Rational(cf)), budget)
    for _ in range(nx):
        budget.charge(1)
        divisor = divisor.mul(QQ(Rational(cm)), budget)
    var resultant_value_q = _rational_image(integer, budget).scale(divisor.inverse(budget), budget)
    var value = DensePoly[QQ].from_sparse(resultant_value_q, budget, variable=Int64(0))
    # u(e, y) = value(y) / F(e, y) in Q(e)[y], by long division: value / G
    # for the monic G = F / lead, then u = (value / G) / lead.
    var inverse = _field_inverse(representations[top], m, budget)
    var remainder = List[DensePoly[QQ]]()
    for c in value.coefficients:
        remainder.append(_constant(c, budget))
    while len(remainder) and remainder[len(remainder) - 1].is_zero():
        _ = remainder.pop()
    var quotient = List[DensePoly[QQ]]()
    if len(remainder):
        if len(remainder) <= top:
            return None
        var monic = List[DensePoly[QQ]]()
        for i in range(top + 1):
            monic.append(_field_mul(representations[i], inverse, m, budget))
        var division = _monic_division(monic, remainder, m, True, budget)
        if not division[0]:
            return None
        for q in division[1]:
            quotient.append(_field_mul(q, inverse, m, budget))
    var order = MonomialOrder()
    var y = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var u = SparsePoly[QQ].zero(order, NoContext(), budget)
    for j in range(len(quotient) - 1, -1, -1):
        budget.charge(1)
        u = u.mul(y, budget).add(_embed_dense(quotient[j], 1, budget), budget)
    var rest = resultant_value_q.sub(_cleared_product(u, definitions[0], budget), budget)
    var v = _divide_cleared(rest, definitions[1], budget)
    if not v:
        return None
    return ResultantResult[QQ](resultant_value_q^, (u^, v.take()))


def _definition_resultant(first: SparsePoly[QQ], element: RealAlgebraic, mut budget: Budget) raises PolyError -> Optional[ResultantResult[QQ]]:
    """`resultant(first, m_e(x), 1)` for first in y (variable zero) and x (variable one), by `_tower_resultant`.

    The rows of first by powers of y are kept unreduced, so its tower
    definitions are first and m_e themselves; None when first has another
    variable or `_tower_resultant` declines.
    """
    var rows = List[List[QQ]]()
    for term in first.terms:
        budget.charge(1)
        var i = 0
        var j = 0
        for pair in term.monomial.pairs:
            if pair[0] == 0:
                i = Int(pair[1])
            elif pair[0] == 1:
                j = Int(pair[1])
            else:
                return None
        while len(rows) <= i:
            rows.append(List[QQ]())
        while len(rows[i]) <= j:
            rows[i].append(QQ.zero(NoContext()))
        rows[i][j] = term.coefficient.copy()
    var representations = List[DensePoly[QQ]]()
    for r in rows:
        representations.append(DensePoly[QQ].from_coefficients(0, NoContext(), r.copy(), budget))
    return _tower_resultant(element, representations, budget)


def _divide_cleared(a: SparsePoly[QQ], b: SparsePoly[QQ], mut budget: Budget) raises PolyError -> Optional[SparsePoly[QQ]]:
    # divide_exact over Q, on integers: for A = s a and B = c b with content g,
    # a / b = (g A / B) c / (s g), and g A / B is integral exactly when b divides a.
    var s = Integer(1)
    for term in a.terms:
        s = _denominator(s, term.coefficient, budget)
    var c = Integer(1)
    for term in b.terms:
        c = _denominator(c, term.coefficient, budget)
    var B = _integer_image(b, budget)
    var g = ZZ.zero(NoContext())
    for term in B.terms:
        budget.charge(1)
        g = _gcd_coefficient(g, term.coefficient, budget)
    if g.is_zero():
        return divide_exact(a, b, budget)
    var quotient = divide_exact(_integer_image(a, budget).scale(g, budget), B, budget)
    if not quotient:
        return None
    var scale = QQ(Rational(c)).mul(QQ(Rational(s)).mul(QQ(Rational(g.value)), budget).inverse(budget), budget)
    return _rational_image(quotient.value(), budget).scale(scale, budget)
