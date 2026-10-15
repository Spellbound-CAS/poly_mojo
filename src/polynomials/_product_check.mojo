"""Consume exact product groups with a small outer merge and no product buffers."""
from .coefficients import CoefficientRing
from .errors import Budget, PolyError, _allocation, _charge_size, _sum
from .monomial import Monomial
from .sparse import SparsePoly, Term, _product_rows
from ._ordering import _sift_min


from ._product_stream import _ProductGroup, _GroupOrder


@always_inline
def _advance_group[C: CoefficientRing](mut group: _ProductGroup[C], index: Int, left: List[SparsePoly[C]], right: List[SparsePoly[C]], remainder: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    if index == 0:
        if group.position == len(remainder.terms):
            return False
        group.term = remainder.terms[group.position]._copy(budget)
        group.position += 1
        return True
    var i = index - 1
    if len(left[i].terms) <= len(right[i].terms):
        return group.next(left[i], right[i], budget)
    return group.next(right[i], left[i], budget)


@always_inline
def _matches_term[C: CoefficientRing](f: SparsePoly[C], term: Term[C], mut position: Int, mut budget: Budget) raises PolyError -> Bool:
    _charge_size(term.coefficient.size(), budget)
    if term.coefficient.is_zero():
        return True
    if position == len(f.terms):
        return False
    ref expected = f.terms[position]
    _charge_size(_sum(_sum(len(expected.monomial.pairs), len(term.monomial.pairs)), _sum(expected.coefficient.size(), term.coefficient.size())), budget)
    if expected != term:
        return False
    position += 1
    return True


def _matches_one[C: CoefficientRing](f: SparsePoly[C], rows: SparsePoly[C], columns: SparsePoly[C], remainder: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    var group = _ProductGroup[C](f.ring)
    var prepared = _product_rows(rows, columns, budget)
    swap(group.cursors, prepared[0])
    swap(group.heap, prepared[1])
    var available = group.next(rows, columns, budget)
    var remaining = 0
    var position = 0
    while available or remaining < len(remainder.terms):
        budget.charge(1)
        var relation = -1 if available else 1
        if available and remaining < len(remainder.terms):
            relation = f.order._compare(group.term.monomial, remainder.terms[remaining].monomial, budget)
        var term = Term[C](Monomial(), C.zero(f.ring))
        if relation <= 0:
            swap(term, group.term)
            if relation == 0:
                term.coefficient = term.coefficient.add(remainder.terms[remaining].coefficient, budget)
                remaining += 1
        else:
            term = remainder.terms[remaining]._copy(budget)
            remaining += 1
        if not _matches_term(f, term, position, budget):
            return False
        if relation <= 0:
            available = group.next(rows, columns, budget)
    return position == len(f.terms)


def _matches_products[C: CoefficientRing](f: SparsePoly[C], left: List[SparsePoly[C]], right: List[SparsePoly[C]], remainder: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    # Validate before entering: an early mismatch must not hide malformed input.
    # Each product folds its own collisions before the outer merge. This keeps
    # the hot row heap as small as multiplication's, even with many cofactors.
    if len(left) == 1:
        if len(left[0].terms) <= len(right[0].terms):
            return _matches_one(f, left[0], right[0], remainder, budget)
        return _matches_one(f, right[0], left[0], remainder, budget)
    var count = _sum(len(left), 1)
    _allocation[_ProductGroup[C]](count, budget)
    _allocation[Int](count, budget)
    var groups = List[_ProductGroup[C]](capacity=count)
    var heap = List[Int](capacity=count)
    for i in range(count):
        budget.charge(1)
        var group = _ProductGroup[C](f.ring)
        if i:
            var j = i - 1
            var prepared = _product_rows(left[j], right[j], budget) if len(left[j].terms) <= len(right[j].terms) else _product_rows(right[j], left[j], budget)
            swap(group.cursors, prepared[0])
            swap(group.heap, prepared[1])
        if _advance_group(group, i, left, right, remainder, budget):
            heap.append(i)
        groups.append(group^)
    for i in range(len(heap) // 2 - 1, -1, -1):
        budget.charge(1)
        _sift_min[_GroupOrder[C]](heap, groups, i, f.order, budget)
    var position = 0
    while len(heap):
        budget.charge(1)
        var term = Term[C](Monomial(), C.zero(f.ring))
        swap(term, groups[heap[0]].term)
        while True:
            var index = heap[0]
            if not _advance_group(groups[index], index, left, right, remainder, budget):
                var last = heap.pop()
                if len(heap):
                    heap[0] = last
            if not len(heap):
                break
            _sift_min[_GroupOrder[C]](heap, groups, 0, f.order, budget)
            ref following = groups[heap[0]].term
            _charge_size(_sum(len(term.monomial.pairs), len(following.monomial.pairs)), budget)
            if term.monomial != following.monomial:
                break
            term.coefficient = term.coefficient.add(following.coefficient, budget)
        if not _matches_term(f, term, position, budget):
            return False
    return position == len(f.terms)
