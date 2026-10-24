"""Exact Sturm counts and dyadic real-root isolation with rational reconstruction."""
from apn_mojo import Rational, Integer, gcd, lcm, div_rem_euclid, div_rem_floor
from .coefficients import QQ, ZZ, NoContext, _work
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .evaluate import _derivative_dense
from .division import divmod
from .certificates import DivisionResult, SturmSequence, RootCountCertificate, _integer_form, _IntegerChain, _integer_variations, _fujiwara_exponent, _SHORT
from ._dense_ops import _copy as _copy_poly
from ._root_ops import _copy, _compare, _midpoint, _variations, _one_sided, _sign_at
from ._images import _integer_dense
from ._integer_ops import _primitive


@fieldwise_init
struct RootInterval(Copyable, Movable, Writable):
    """An open rational interval containing exactly one distinct real root."""
    var low: QQ
    """The lower endpoint, strictly below the root."""
    var high: QQ
    """The upper endpoint, strictly above the root."""


@fieldwise_init
struct RootIsolation(Copyable, Movable, Writable):
    """The distinct real roots as exact rationals and disjoint open intervals."""
    var exact_roots: List[QQ]
    """All rational roots in increasing order, without multiplicities."""
    var intervals: List[RootInterval]
    """Increasing disjoint intervals for the remaining roots, excluding exact roots."""


def sturm_sequence(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> SturmSequence:
    """Construct the unnormalized Sturm sequence and every quotient witness.

    Args:
        f: A nonzero rational polynomial; constants give a one-element sequence.
        budget: The shared work allowance.

    Returns:
        The sequence starting with f and f-prime, including its last nonzero entry.

    Raises:
        PolyError: DOMAIN for zero or malformed input, LIMIT or CAPACITY on work limits.
    """
    var chain = _IntegerChain()
    return _sturm_chain(f, chain, budget, integer=False)


def _sturm_chain(f: DensePoly[QQ], mut chain: _IntegerChain, mut budget: Budget, *, integer: Bool = True) raises PolyError -> SturmSequence:
    # Without integer data (integer=False) the chain stays empty, small
    # rational entries divide as rationals, and integer forms are computed
    # only when the entries grow.
    f.check_invariants(budget)
    if f.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "zero polynomial has no finite root count")
    _allocation[DensePoly[QQ]](1, budget)
    var polynomials: List[DensePoly[QQ]] = [_copy_poly(f, budget)]
    var quotients = List[DensePoly[QQ]]()
    if f.degree() > 0:
        _allocation[DensePoly[QQ]](1, budget)
        polynomials.append(_derivative_dense(f, budget))
        # Each entry is c P with P a primitive integer polynomial and c a
        # rational scalar. rem(c0 P0, c1 P1) = c0 prem(P0, P1) / lc(P1)^(d+1),
        # so the primitive PRS over the integers gives every remainder and
        # quotient of the rational Euclidean sequence exactly, without
        # normalizing rationals inside the division.
        var before = (List[Integer](), Rational(1))
        var current = (List[Integer](), Rational(1))
        if integer:
            before = _primitive_form(polynomials[0], budget)
            current = _primitive_form(polynomials[1], budget)
            chain.first = before[0].copy()
            chain.entries = [before[0].copy(), current[0].copy()]
            chain.bounded = True
            chain.bound = _fujiwara_exponent(before[0])
            chain.second = current[0].copy()
            chain.degrees = [len(before[0]) - 1, len(current[0]) - 1]
            chain.signs = [before[1].sign(), current[1].sign()]
        var forms = integer
        while True:
            budget.charge(1)
            var n = len(polynomials)
            if not integer and _small(polynomials[n - 1], budget):
                var step = _rational_divmod(polynomials[n - 2], polynomials[n - 1], budget)
                var quotient = DensePoly[QQ](f.variable, NoContext())
                swap(quotient, step.quotient)
                _allocation[DensePoly[QQ]](1, budget)
                quotients.append(quotient^)
                if step.remainder.is_zero():
                    break
                _allocation[DensePoly[QQ]](1, budget)
                polynomials.append(step.remainder.neg(budget))
                forms = False
                continue
            if not forms:
                before = _primitive_form(polynomials[n - 2], budget)
                current = _primitive_form(polynomials[n - 1], budget)
                forms = True
            if integer and _rational_cheaper(polynomials[n - 2], polynomials[n - 1], before[0], current[0], budget):
                # Small rational entries divide faster than the integer
                # forms, whose pseudo-remainder carries a power of lc.
                var step = _rational_divmod(polynomials[n - 2], polynomials[n - 1], budget)
                var quotient = DensePoly[QQ](f.variable, NoContext())
                swap(quotient, step.quotient)
                var finished = step.remainder.is_zero()
                if not finished:
                    var following = step.remainder.neg(budget)
                    var form = _primitive_form(following, budget)
                    _integer_step(chain, quotient, before[1], current[1], form[1], budget)
                    chain.degrees.append(len(form[0]) - 1)
                    chain.signs.append(form[1].sign())
                    if len(chain.entries) <= _SHORT:
                        chain.entries.append(form[0].copy())
                    _allocation[DensePoly[QQ]](1, budget)
                    polynomials.append(following^)
                    before = current^
                    current = form^
                _allocation[DensePoly[QQ]](1, budget)
                quotients.append(quotient^)
                if finished:
                    break
                continue
            var division = _pseudo_divide(before[0], current[0], budget)
            var quotient_scale: Rational
            var next_scale: Rational
            var content: Integer
            var primitive: List[Integer]
            budget.charge(UInt64(len(division[1])))
            _allocation[Integer](len(division[1]), budget)
            try:
                var parts = _content_parts(division[1])
                content = parts[0].copy()
                primitive = parts[1].copy()
                var power = Rational(division[2])
                quotient_scale = before[1] / (current[1] * power)
                next_scale = -(before[1] * Rational(content)) / power
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            _allocation[DensePoly[QQ]](1, budget)
            quotients.append(_scaled_poly(division[0], quotient_scale, f.variable, budget))
            if len(division[1]) == 0:
                break
            if integer:
                chain.powers.append(division[2].copy())
                chain.quotients.append(division[0].copy())
                chain.contents.append(content.copy())
                chain.degrees.append(len(primitive) - 1)
                chain.signs.append(next_scale.sign())
                if len(chain.entries) <= _SHORT:
                    chain.entries.append(primitive.copy())
            _allocation[DensePoly[QQ]](1, budget)
            polynomials.append(_scaled_poly(primitive, next_scale, f.variable, budget))
            before = current^
            current = (primitive^, next_scale^)
    if integer and len(chain.entries) != len(chain.degrees):
        chain.entries = List[List[Integer]]()
    return SturmSequence(polynomials^, quotients^)


comptime _LAZY_QUOTIENT = 8
"""Quotient lengths from which a rational division runs on integers over a common denominator."""


def _rational_divmod(f: DensePoly[QQ], g: DensePoly[QQ], mut budget: Budget) raises PolyError -> DivisionResult[QQ]:
    """divmod over QQ; long quotients run on integers over a common denominator.

    With f = cf F and g = cg G for primitive integer F, G and L = lc(G),
    the j-th remainder of F by G times L^j is an integer polynomial. Each
    entry keeps the power of L it was last scaled to, so a step rescales
    and updates only the entries under G's nonzero coefficients; quotient
    and remainder coefficients are normalized once, at the end. The result
    is the unique quotient and remainder, as from divmod.
    """
    var m = f.degree() - g.degree() + 1
    if m < _LAZY_QUOTIENT or g.degree() < 0:
        return divmod(f, g, budget)
    var ff = _primitive_form(f, budget)
    var gg = _primitive_form(g, budget)
    var n = len(gg[0]) - 1
    _charge_size(m * (n + 1), budget)
    _allocation[Integer](len(ff[0]) + 2 * m, budget)
    _allocation[QQ](m + n, budget)
    var parts: Tuple[List[QQ], List[QQ]]
    try:
        parts = _lazy_division(ff[0], ff[1], gg[0], gg[1], m)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return DivisionResult[QQ](DensePoly[QQ].from_coefficients(f.variable, NoContext(), parts[0].copy(), budget), DensePoly[QQ].from_coefficients(f.variable, NoContext(), parts[1].copy(), budget))


def _lazy_division(big_f: List[Integer], cf: Rational, big_g: List[Integer], cg: Rational, m: Int) raises -> Tuple[List[QQ], List[QQ]]:
    """Quotient and remainder coefficients of (cf F) / (cg G) for _rational_divmod."""
    var n = len(big_g) - 1
    var stored = big_f.copy()
    var scaled = List[Int](length=len(stored), fill=0)
    var lead = big_g[n].copy()
    var powers = List[Integer]()
    powers.append(Integer(1))
    var numerators = List[Integer](length=m, fill=Integer(0))
    for j in range(m):
        powers.append(powers[j] * lead)
        var top = len(stored) - 1 - j
        # value = L^j R_j[top]; the quotient coefficient is value / L^(j + 1).
        var value = stored[top] * powers[j - scaled[top]]
        if not value:
            continue
        var shift = top - n
        numerators[shift] = value.copy()
        for i in range(n):
            if not big_g[i]:
                continue
            var k = shift + i
            stored[k] = stored[k] * powers[j + 1 - scaled[k]] - value * big_g[i]
            scaled[k] = j + 1
    # f / g = (cf / cg) F / G and the remainder is cf R.
    var ratio = cf / cg
    var quotient = List[QQ](capacity=m)
    for k in range(m):
        # The coefficient of x^k came from step j = m - 1 - k, over L^(j + 1).
        quotient.append(QQ(ratio * Rational(numerators[k], powers[m - k])) if numerators[k] else QQ.zero(NoContext()))
    var remainder = List[QQ](capacity=n)
    for i in range(n):
        remainder.append(QQ(cf * Rational(stored[i], powers[scaled[i]])) if stored[i] else QQ.zero(NoContext()))
    return (quotient^, remainder^)


def _small(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> Bool:
    # A leading coefficient of at most 512 bytes: rational entries grow
    # with their leading coefficients, and generic inputs outgrow this fast.
    budget.charge(1)
    return p.coefficients[p.degree()].size() <= 512


def _rational_cheaper(a: DensePoly[QQ], b: DensePoly[QQ], pa: List[Integer], pb: List[Integer], mut budget: Budget) raises PolyError -> Bool:
    # Compare the rational entries' sizes with the integer forms' sizes plus
    # the leading-coefficient powers a pseudo-remainder multiplies in.
    budget.charge(UInt64(len(a.coefficients) + len(b.coefficients)))
    var rational = 0
    for c in a.coefficients:
        rational += c.size()
    for c in b.coefficients:
        rational += c.size()
    var integer = 0
    for v in pa:
        integer += v.magnitude_bit_length() // 8 + 1
    for v in pb:
        integer += v.magnitude_bit_length() // 8 + 1
    var lead = pb[len(pb) - 1].magnitude_bit_length() // 8 + 1
    integer += 2 * lead * len(pa)
    return rational < integer


def _integer_step(mut chain: _IntegerChain, quotient: DensePoly[QQ], c0: Rational, c1: Rational, c2: Rational, mut budget: Budget) raises PolyError:
    # From c0 P0 = q c1 P1 - c2 P2: with D clearing the denominators of
    # q c1 / c0 and c2 / c0, D P0 = (D q c1 / c0) P1 + (-D c2 / c0) P2.
    budget.charge(UInt64(len(quotient.coefficients) + 1))
    try:
        var ratio = c1 / c0
        var tail = c2 / c0
        var scaled = List[Rational](capacity=len(quotient.coefficients))
        var common = tail.denominator()
        for c in quotient.coefficients:
            var value = c.value * ratio
            var d = value.denominator()
            if common % d:
                common = lcm(common, d)
            scaled.append(value^)
        var integer = List[Integer](capacity=len(scaled))
        for value in scaled:
            integer.append(value.numerator() * (common // value.denominator()))
        chain.powers.append(common.copy())
        chain.quotients.append(integer^)
        chain.contents.append(-(tail.numerator() * (common // tail.denominator())))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _primitive_form(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> Tuple[List[Integer], Rational]:
    # p = c P with P primitive over the integers and c = content / common denominator.
    var form = _integer_form(p, budget)
    var content: Integer
    try:
        content = _content_of(form[0])
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var primitive = _divide_all(form[0], content, budget)
    try:
        return (primitive^, Rational(content, form[1]))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _content_of(values: List[Integer]) raises -> Integer:
    var result = Integer(0)
    for v in values:
        if v:
            result = gcd(result, v)
            if result == Integer(1):
                break
    return result^


def _content_parts(values: List[Integer]) raises -> Tuple[Integer, List[Integer]]:
    # The positive content and the primitive part in one pass: a gcd of two
    # coefficients is usually the content already, and the exact quotients
    # that confirm it are the primitive part; a nonzero remainder lowers it.
    var result = List[Integer](capacity=len(values))
    if len(values) == 0:
        return (Integer(0), result^)
    var content = gcd(values[len(values) - 1], values[0])
    if not content:
        content = abs(values[len(values) - 1])
    while True:
        result.clear()
        var lowered = False
        for v in values:
            if content == Integer(1):
                result.append(v.copy())
                continue
            var parts = div_rem_floor(v, content)
            if parts[1]:
                content = gcd(content, v)
                lowered = True
                break
            result.append(parts[0].copy())
        if not lowered:
            return (content^, result^)


def _divide_all(values: List[Integer], divisor: Integer, mut budget: Budget) raises PolyError -> List[Integer]:
    budget.charge(UInt64(len(values)))
    _allocation[Integer](len(values), budget)
    var result = List[Integer](capacity=len(values))
    try:
        for v in values:
            result.append(v.copy() if divisor == Integer(1) else v // divisor)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return result^


def _pseudo_divide(a: List[Integer], b: List[Integer], mut budget: Budget) raises PolyError -> Tuple[List[Integer], List[Integer], Integer]:
    # lc(b)^(d+1) a == q b + r with deg r < deg b, d = deg a - deg b; r is trimmed.
    var db = len(b) - 1
    var delta = len(a) - 1 - db
    budget.charge(UInt64((delta + 1) * len(a)))
    _allocation[Integer](len(a) + delta + 1, budget)
    var r = a.copy()
    var q = List[Integer](length=delta + 1, fill=Integer(0))
    var size = 0
    try:
        ref lc = b[db]
        var bits = lc.magnitude_bit_length()
        var shift = bits - 1 if lc.sign() > 0 and (Integer(1) << (bits - 1)) == lc else -1
        for k in range(delta, -1, -1):
            var t = r[db + k].copy()
            for j in range(k + 1, delta + 1):
                q[j] = q[j] * lc if shift < 0 else q[j] << shift
            q[k] = t.copy()
            for j in range(db + k):
                if r[j]:
                    r[j] = r[j] * lc if shift < 0 else r[j] << shift
            r[db + k] = Integer(0)
            if t:
                for j in range(db):
                    r[k + j] = r[k + j] - t * b[j]
            size += r[0].magnitude_bit_length() // 8 + 1
        var length = db
        while length > 0 and not r[length - 1]:
            length -= 1
        r.resize(length, Integer(0))
        var power = Integer(1)
        for _ in range(delta + 1):
            power = power * lc
        _charge_size(size, budget)
        return (q^, r^, power^)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _scaled_poly(values: List[Integer], scale: Rational, variable: Int64, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # The canonical coefficients scale * values[j]: with scale = n / d reduced
    # and g = gcd(values[j], d), n (values[j] / g) and d / g are coprime.
    budget.charge(UInt64(len(values)))
    _allocation[QQ](len(values), budget)
    var result = List[QQ](capacity=len(values))
    var size = 0
    try:
        var n = scale.numerator()
        var d = scale.denominator()
        var unit = d == Integer(1)
        for v in values:
            if not v:
                result.append(QQ.zero(NoContext()))
                continue
            if unit:
                result.append(QQ(Rational(n * v)))
                continue
            var g = gcd(v, d)
            if g == Integer(1):
                result.append(QQ(Rational(_numerator=n * v, _denominator=d.copy())))
            else:
                result.append(QQ(Rational(_numerator=n * (v // g), _denominator=d // g)))
            size += result[len(result) - 1].size()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    _charge_size(size, budget)
    return DensePoly[QQ]._from_canonical(variable, NoContext(), result^)


def sign_changes(sequence: SturmSequence, a: QQ, mut budget: Budget) raises PolyError -> Int:
    """Count changes of sign at an exact point, omitting zero evaluations.

    Args:
        sequence: The polynomial sequence to evaluate.
        a: The exact evaluation point.
        budget: The shared work allowance.

    Returns:
        The number of adjacent opposite nonzero signs.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    return _variations(sequence.polynomials, a, budget)


def count_roots(f: DensePoly[QQ], a: QQ, b: QQ, mut budget: Budget) raises PolyError -> RootCountCertificate:
    """Count distinct real roots in an open interval with nonroot endpoints.

    Args:
        f: A nonzero polynomial, with repeated roots allowed.
        a: The lower endpoint.
        b: The upper endpoint.
        budget: The shared work allowance.

    Returns:
        The count and a complete Sturm certificate.

    Raises:
        PolyError: DOMAIN for zero, endpoint roots, or a not below b; LIMIT or CAPACITY on work limits.
    """
    if _compare(a, b, budget) >= 0 or _sign_at(f, a, budget) == 0 or _sign_at(f, b, budget) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "root count requires ordered nonroot endpoints")
    var integer = _IntegerChain()
    var sequence = _sturm_chain(f, integer, budget)
    var count = _integer_variations(integer, sequence, a, budget) - _integer_variations(integer, sequence, b, budget)
    return RootCountCertificate(sequence^, _copy(a, budget), _copy(b, budget), count)


def _rational_root(f: DensePoly[QQ], low: QQ, high: QQ, bound: ZZ, mut budget: Budget) raises PolyError -> Optional[QQ]:
    # Rational-root theorem: the reduced denominator divides the primitive
    # leading coefficient A. Width < 1/(2*A*A) makes any root of denominator
    # <= A a continued-fraction convergent of the midpoint (Legendre).
    # Exact evaluation admits only actual roots, including negative ones.
    var middle = _midpoint(low, high, budget)
    _charge_size(middle.size(), budget)
    var numerator = ZZ(middle.value.numerator())
    var denominator = ZZ(middle.value.denominator())
    var h0 = ZZ.zero(NoContext())
    var h1 = ZZ.one(NoContext())
    var k0 = ZZ.one(NoContext())
    var k1 = ZZ.zero(NoContext())
    while not denominator.is_zero():
        budget.charge(1)
        _work(numerator.size(), denominator.size(), True, budget)
        var a: ZZ
        var remainder: ZZ
        try:
            var qr = div_rem_euclid(numerator.value, denominator.value)
            a = ZZ(qr[0].copy())
            remainder = ZZ(qr[1].copy())
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        var h = a.mul(h1, budget).add(h0, budget)
        var k = a.mul(k1, budget).add(k0, budget)
        _work(k.size(), bound.size(), False, budget)
        if k.value > bound.value:
            break
        _work(h.size(), k.size(), True, budget)
        var candidate: QQ
        try:
            candidate = QQ(Rational(h.value, k.value))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        if _compare(low, candidate, budget) < 0 and _compare(candidate, high, budget) < 0 and _sign_at(f, candidate, budget) == 0:
            return candidate^
        h0 = h1^
        h1 = h^
        k0 = k1^
        k1 = k^
        numerator = denominator^
        denominator = remainder^
    return None


def isolate_real_roots(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> RootIsolation:
    """Isolate every distinct real root by Descartes' rule of signs.

    Args:
        f: A nonzero rational polynomial; constants have no roots.
        budget: The shared work allowance.

    Returns:
        Every rational root exactly and sorted disjoint intervals of width below 2**-40 for the others.

    Raises:
        PolyError: DOMAIN for zero or malformed input, LIMIT or CAPACITY on work limits.
    """
    from .descartes import _root_isolation
    return _root_isolation(f, budget)


def refine(f: DensePoly[QQ], interval: RootInterval, width: QQ, mut budget: Budget) raises PolyError -> RootInterval:
    """Narrow a validated isolating interval to strictly less than the requested width.

    The interval is validated by Descartes' rule of signs on the squarefree
    part and narrowed by quadratic interval refinement.

    Args:
        f: The polynomial defining the root.
        interval: An open interval with exactly one root and nonroot endpoints.
        width: The positive target width.
        budget: The shared work allowance.

    Returns:
        A narrower open interval containing the same root, including rational roots.

    Raises:
        PolyError: DOMAIN for nonpositive width or an invalid isolating interval, LIMIT or CAPACITY on work limits.
    """
    if width.value.sign() <= 0:
        raise PolyError(PolyErrorKind.DOMAIN, "refinement requires a positive width")
    from .descartes import _refine_interval
    return _refine_interval(f, interval, width, budget)


def sign_on_interval(f: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> Optional[Int]:
    """Determine the strict sign on a closed interval when no root occurs there.

    Args:
        f: The rational polynomial.
        low: The lower endpoint.
        high: The upper endpoint, possibly equal to low.
        budget: The shared work allowance.

    Returns:
        Positive or negative one for a root-free interval, otherwise None.

    Raises:
        PolyError: DOMAIN for reversed endpoints or malformed input, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    var comparison = _compare(low, high, budget)
    if comparison > 0:
        raise PolyError(PolyErrorKind.DOMAIN, "sign interval has reversed endpoints")
    var sign = _sign_at(f, low, budget)
    if sign == 0 or _sign_at(f, high, budget) == 0:
        return None
    if comparison < 0 and count_roots(f, low, high, budget).count != 0:
        return None
    return sign
