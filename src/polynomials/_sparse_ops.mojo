"""Ordered sparse construction and fused monomial-multiple subtraction."""
from .coefficients import CoefficientRing
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _charge_size
from ._coefficient_ops import _copy
from .sparse import SparsePoly, Term
from .monomial import Monomial, MonomialOrder


def _from_descending[C: CoefficientRing](var terms: List[Term[C]], exemplar: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    # Reduction emits distinct nonzero terms in strictly decreasing order.
    for i in range(len(terms) // 2):
        budget.charge(1)
        var temporary = Term[C](Monomial(), C.zero(exemplar.ring))
        swap(temporary, terms[i])
        swap(temporary, terms[len(terms) - 1 - i])
        swap(temporary, terms[i])
    var result = SparsePoly[C].zero(exemplar.order, exemplar.ring, budget)
    result.terms = terms^
    return result^


def _subtract_multiple[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], term: Term[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    f._same(g, budget)
    if term.coefficient.context() != f.ring:
        raise PolyError(PolyErrorKind.DOMAIN, "multiple in another coefficient ring")
    if term.coefficient.is_zero() or g.is_zero():
        return f._copy(budget)
    if f.order.kind == MonomialOrder.KERNEL:
        return f.sub(g.mul(SparsePoly[C].from_terms([term._copy(budget)], f.order, f.ring, budget), budget), budget)
    var result = SparsePoly[C].zero(f.order, f.ring, budget)
    var count = _sum(len(f.terms), len(g.terms))
    _allocation[Term[C]](count, budget)
    result.terms.reserve(count)
    var i = 0
    # A fixed monomial translation preserves every admissible order. Merge
    # that stream directly, retaining only one translated term at a time.
    for source in g.terms:
        budget.charge(1)
        var coefficient = source.coefficient.mul(term.coefficient, budget)
        if coefficient.is_zero():
            continue
        var monomial = source.monomial.mul(term.monomial, budget)
        var comparison = -1
        while i < len(f.terms):
            budget.charge(1)
            comparison = f.order._compare(f.terms[i].monomial, monomial, budget)
            if comparison >= 0:
                break
            result.terms.append(f.terms[i]._copy(budget))
            i += 1
        if i < len(f.terms) and comparison == 0:
            coefficient = f.terms[i].coefficient.sub(coefficient, budget)
            i += 1
        else:
            coefficient = coefficient.neg(budget)
        if not coefficient.is_zero():
            result.terms.append(Term[C](monomial^, coefficient^))
    while i < len(f.terms):
        budget.charge(1)
        result.terms.append(f.terms[i]._copy(budget))
        i += 1
    return result^


def _subtract_multiples[C: CoefficientRing](f: SparsePoly[C], left: Term[C], g: SparsePoly[C], right: Term[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    # Two translated, sorted streams, retaining at most one term from each.
    f._same(g, budget)
    if left.coefficient.context() != f.ring or right.coefficient.context() != f.ring:
        raise PolyError(PolyErrorKind.DOMAIN, "multiples in another coefficient ring")
    if f.order.kind == MonomialOrder.KERNEL:
        var a = SparsePoly[C].from_terms([left._copy(budget)], f.order, f.ring, budget)
        var b = SparsePoly[C].from_terms([right._copy(budget)], f.order, f.ring, budget)
        return f.mul(a, budget).sub(g.mul(b, budget), budget)
    f.order.validate(left.monomial, budget)
    f.order.validate(right.monomial, budget)
    var result = SparsePoly[C].zero(f.order, f.ring, budget)
    var count = _sum(len(f.terms), len(g.terms))
    _allocation[Term[C]](count, budget)
    result.terms.reserve(count)
    var i = len(f.terms) if left.coefficient.is_zero() else 0
    var j = len(g.terms) if right.coefficient.is_zero() else 0
    _charge_size(_sum(left.coefficient.size(), right.coefficient.size()), budget)
    var left_unit = left.coefficient == C.one(f.ring)
    var right_unit = right.coefficient == C.one(f.ring)
    var a = Optional[Term[C]](None)
    var b = Optional[Term[C]](None)
    while i < len(f.terms) or j < len(g.terms) or a or b:
        budget.charge(1)
        if not a and i < len(f.terms):
            var coefficient = _copy(f.terms[i].coefficient, budget) if left_unit else f.terms[i].coefficient.mul(left.coefficient, budget)
            if not coefficient.is_zero():
                a = Term[C](f.terms[i].monomial.mul(left.monomial, budget), coefficient^)
            i += 1
            if not a:
                continue
        if not b and j < len(g.terms):
            var coefficient = _copy(g.terms[j].coefficient, budget) if right_unit else g.terms[j].coefficient.mul(right.coefficient, budget)
            if not coefficient.is_zero():
                b = Term[C](g.terms[j].monomial.mul(right.monomial, budget), coefficient^)
            j += 1
            if not b:
                continue
        var comparison = f.order._compare(a.value().monomial, b.value().monomial, budget) if Bool(a) and Bool(b) else (-1 if Bool(a) else 1)
        if comparison < 0:
            result.terms.append(a.take())
        elif comparison > 0:
            var term = b.take()
            term.coefficient = term.coefficient.neg(budget)
            result.terms.append(term^)
        else:
            var term = a.take()
            var other = b.take()
            term.coefficient = term.coefficient.sub(other.coefficient, budget)
            if not term.coefficient.is_zero():
                result.terms.append(term^)
    return result^
