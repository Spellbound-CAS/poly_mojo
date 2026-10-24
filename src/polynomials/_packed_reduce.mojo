"""Division by a basis and trace replay on packed exponent vectors.

Both make the choices of their generic counterparts and change only storage,
as `_packed_groebner` does. Division takes the terms of the pending
polynomial in decreasing order and reduces each by the first divisor, in
list order, whose leading monomial divides it, so the quotients and the
remainder are those of `division.divide_by_list`. Over QQ the division is
fraction-free, as in the packed rational engine: the pending polynomial is
an integer polynomial over a common denominator T, a reduction by an element
with integer leading coefficient a scales it by a / gcd(a, c) instead of
dividing by a, and each quotient and remainder term becomes a rational once,
when it is made. Scaling never changes which monomials are present, so the
choices, and the results, are those of the rational division. Trace replay
(`groebner._packed_replay`) sums each step's multiples with `_add_multiple`
in geobuckets.
"""
from apn_mojo import Integer, Rational
from .coefficients import CoefficientField, ModInt, NoContext, QQ, ZZ, _work
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term
from .certificates import ListDivision
from ._sparse_ops import _from_descending
from ._coefficient_ops import _exact, _gcd_coefficient
from ._images import _denominator, _scaled_integer
from ._packed_groebner import _Buckets, _Packed, _Term, _packed, _mul, _div, _divides, _equal, _pack, _monomial, _priorities


def _pack_monomial[W: Int](m: Monomial, order: MonomialOrder, mut budget: Budget) raises PolyError -> _Packed[W]:
    var exponents = SIMD[DType.uint32, W](0)
    for pair in m.pairs:
        budget.charge(1)
        var lane = -1
        for k in range(len(order.variables)):
            if order.variables[k] == pair[0]:
                lane = k
                break
        if lane < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "variable absent from order")
        exponents[lane] = pair[1]
    return _packed(exponents)


def _add_tail[C: CoefficientField, W: Int](mut pending: _Buckets[C, W, False], source: List[_Term[C, W]], shift: _Packed[W], factor: C, mut budget: Budget) raises PolyError:
    """Add every term of `source` but its leading one, times shift * factor."""
    comptime if C == ModInt:
        pending.add_multiple(source, shift, rebind[ModInt](factor), budget)
    else:
        var n = len(source) - 1
        if n <= 0:
            return
        budget.charge(UInt64(n))
        _allocation[_Term[C, W]](n, budget)
        var terms = List[_Term[C, W]](capacity=n)
        for i in range(n):
            ref t = source[i]
            terms.append(_Term[C, W](_mul(t.monomial, shift), t.coefficient.mul(factor, budget)))
        pending.add(terms^, budget)


def _add_multiple[C: CoefficientField, W: Int](mut pending: _Buckets[C, W, False], source: List[_Term[C, W]], shift: _Packed[W], factor: C, mut budget: Budget) raises PolyError:
    """Add all of `source` times shift * factor."""
    if len(source) == 0:
        return
    _add_tail(pending, source, shift, factor, budget)
    ref lead = source[len(source) - 1]
    _allocation[_Term[C, W]](1, budget)
    var top = List[_Term[C, W]](capacity=1)
    top.append(_Term[C, W](_mul(lead.monomial, shift), lead.coefficient.mul(factor, budget)))
    pending.add(top^, budget)


def _drain[C: CoefficientField, W: Int](mut pending: _Buckets[C, W, False], mut budget: Budget) raises PolyError -> List[_Term[C, W]]:
    """The pending polynomial's terms in increasing order, emptying it."""
    var terms = List[_Term[C, W]]()
    while True:
        var next = pending.leading(budget)
        if not next:
            break
        _allocation[_Term[C, W]](1, budget)
        terms.append(next.take())
    terms.reverse()
    return terms^


def _packed_divide[C: CoefficientField, W: Int](f: SparsePoly[C], divisors: List[SparsePoly[C]], mut budget: Budget) raises PolyError -> ListDivision[C]:
    """`divide_by_list(f, divisors)` for nonzero divisors in f's ring and order."""
    ref order = f.order
    var kind = order.kind
    _allocation[List[_Term[C, W]]](len(divisors), budget)
    var packed = List[List[_Term[C, W]]](capacity=len(divisors))
    var leads = List[_Packed[W]](capacity=len(divisors))
    var inverses = List[C](capacity=len(divisors))
    for g in divisors:
        budget.charge(1)
        var p = _pack[C, W](g, order, budget)
        leads.append(p[len(p) - 1].monomial)
        inverses.append(p[len(p) - 1].coefficient.inverse(budget))
        packed.append(p^)
    var pending = _Buckets[C, W, False](kind, f.ring)
    pending.add(_pack[C, W](f, order, budget), budget)
    var by_index = _priorities(order, budget)
    var quotient_terms = List[List[Term[C]]](capacity=len(divisors))
    for _ in range(len(divisors)):
        quotient_terms.append(List[Term[C]]())
    var remainder = List[Term[C]]()
    while True:
        var next = pending.leading(budget)
        if not next:
            break
        var top = next.take()
        var reducer = -1
        for i in range(len(leads)):
            if _divides(leads[i], top.monomial):
                reducer = i
                break
        budget.charge(UInt64(reducer + 1 if reducer >= 0 else len(leads)))
        _allocation[Term[C]](1, budget)
        if reducer < 0:
            remainder.append(Term[C](_monomial(top.monomial, order, by_index, budget), top.coefficient.copy()))
            continue
        var shift = _div(top.monomial, leads[reducer])
        var factor = top.coefficient.mul(inverses[reducer], budget)
        _add_tail(pending, packed[reducer], shift, factor.neg(budget), budget)
        quotient_terms[reducer].append(Term[C](_monomial(shift, order, by_index, budget), factor^))
    _allocation[SparsePoly[C]](len(divisors), budget)
    var quotients = List[SparsePoly[C]](capacity=len(divisors))
    for i in range(len(divisors)):
        budget.charge(1)
        var terms = List[Term[C]]()
        swap(terms, quotient_terms[i])
        quotients.append(_from_descending(terms^, f, budget))
    return ListDivision[C](quotients^, _from_descending(remainder^, f, budget))


def _integral[W: Int](p: SparsePoly[QQ], order: MonomialOrder, mut budget: Budget) raises PolyError -> Tuple[List[_Term[ZZ, W]], ZZ]:
    """p times the least common denominator D of its coefficients, and D."""
    var rational = _pack[QQ, W](p, order, budget)
    var common = Integer(1)
    for t in rational:
        budget.charge(1)
        common = _denominator(common, t.coefficient, budget)
    _allocation[_Term[ZZ, W]](len(rational), budget)
    var terms = List[_Term[ZZ, W]](capacity=len(rational))
    for t in rational:
        budget.charge(1)
        terms.append(_Term[ZZ, W](t.monomial, _scaled_integer(t.coefficient, common, budget)))
    return (terms^, ZZ(common^))


def _ratio(numerator: ZZ, denominator: ZZ, mut budget: Budget) raises PolyError -> QQ:
    _work(numerator.size(), denominator.size(), True, budget)
    try:
        return QQ(Rational(numerator.value.copy(), denominator.value.copy()))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _packed_divide_rational[W: Int](f: SparsePoly[QQ], divisors: List[SparsePoly[QQ]], mut budget: Budget) raises PolyError -> ListDivision[QQ]:
    """`divide_by_list(f, divisors)` over QQ, fraction-free, for nonzero divisors in f's order."""
    ref order = f.order
    var ring = NoContext()
    # Divisor i is integral[i] / scales[i] with integral[i] over ZZ.
    _allocation[List[_Term[ZZ, W]]](len(divisors), budget)
    var integral = List[List[_Term[ZZ, W]]](capacity=len(divisors))
    var denominators = List[ZZ](capacity=len(divisors))
    var leads = List[_Packed[W]](capacity=len(divisors))
    for g in divisors:
        budget.charge(1)
        var cleared = _integral[W](g, order, budget)
        leads.append(cleared[0][len(cleared[0]) - 1].monomial)
        integral.append(cleared[0].copy())
        denominators.append(cleared[1].copy())
    # The pending polynomial is P / T for an integer polynomial P; T > 0.
    var start = _integral[W](f, order, budget)
    var total = start[1].copy()
    var pending = _Buckets[ZZ, W, True](order.kind, ring)
    pending.add(start[0].copy(), budget)
    var one = ZZ.one(ring)
    var by_index = _priorities(order, budget)
    var quotient_terms = List[List[Term[QQ]]](capacity=len(divisors))
    for _ in range(len(divisors)):
        quotient_terms.append(List[Term[QQ]]())
    var remainder = List[Term[QQ]]()
    while True:
        var next = pending.leading(budget)
        if not next:
            break
        var top = next.take()
        var reducer = -1
        for i in range(len(leads)):
            if _divides(leads[i], top.monomial):
                reducer = i
                break
        budget.charge(UInt64(reducer + 1 if reducer >= 0 else len(leads)))
        _allocation[Term[QQ]](1, budget)
        if reducer < 0:
            remainder.append(Term[QQ](_monomial(top.monomial, order, by_index, budget), _ratio(top.coefficient, total, budget)))
            continue
        # P <- s P - e m g with s = a / d and e = c / d for d = gcd(a, c),
        # signed so that s > 0; the term e D / (s T) joins quotient i.
        ref integral_g = integral[reducer]
        ref lead = integral_g[len(integral_g) - 1].coefficient
        var common = _gcd_coefficient(lead, top.coefficient, budget)
        if lead.value.sign() < 0:
            common = common.neg(budget)
        var scale = _exact(lead, common, budget)
        var factor = _exact(top.coefficient, common, budget)
        if scale != one:
            pending.scale(scale, budget)
            total = total.mul(scale, budget)
        var shift = _div(top.monomial, leads[reducer])
        var n = len(integral_g) - 1
        if n > 0:
            budget.charge(UInt64(n))
            _allocation[_Term[ZZ, W]](n, budget)
            var terms = List[_Term[ZZ, W]](capacity=n)
            var negated = factor.neg(budget)
            for k in range(n):
                ref t = integral_g[k]
                terms.append(_Term[ZZ, W](_mul(t.monomial, shift), t.coefficient.mul(negated, budget)))
            pending.add(terms^, budget)
        quotient_terms[reducer].append(Term[QQ](_monomial(shift, order, by_index, budget), _ratio(factor.mul(denominators[reducer], budget), total, budget)))
    _allocation[SparsePoly[QQ]](len(divisors), budget)
    var quotients = List[SparsePoly[QQ]](capacity=len(divisors))
    for i in range(len(divisors)):
        budget.charge(1)
        var terms = List[Term[QQ]]()
        swap(terms, quotient_terms[i])
        quotients.append(_from_descending(terms^, f, budget))
    return ListDivision[QQ](quotients^, _from_descending(remainder^, f, budget))
