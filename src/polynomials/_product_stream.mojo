"""Ordered product groups shared by exact producers and independent consumers."""
from .coefficients import CoefficientRing
from .errors import Budget, PolyError, _allocation, _charge_size, _sum
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term, _ProductCursor, _ProductOrder, _advance_product, _product_rows
from ._ordering import _IndexedOrder, _sift_min


struct _ProductGroup[C: CoefficientRing](Copyable, Movable):
    var cursors: List[_ProductCursor[Self.C]]
    var heap: List[Int]
    var term: Term[Self.C]
    var position: Int

    def __init__(out self, ring: Self.C.Context):
        self.cursors = List[_ProductCursor[Self.C]]()
        self.heap = List[Int]()
        self.term = Term[Self.C](Monomial(), Self.C.zero(ring))
        self.position = 0

    @always_inline
    def next(mut self, rows: SparsePoly[Self.C], columns: SparsePoly[Self.C], mut budget: Budget) raises PolyError -> Bool:
        while len(self.heap):
            budget.charge(1)
            var term = Term[Self.C](Monomial(), Self.C.zero(rows.ring))
            swap(term, self.cursors[self.heap[0]].term)
            while True:
                var index = self.heap[0]
                if not _advance_product(self.cursors[index], rows, columns, budget):
                    var last = self.heap.pop()
                    if len(self.heap):
                        self.heap[0] = last
                if not len(self.heap):
                    break
                _sift_min[_ProductOrder[Self.C]](self.heap, self.cursors, 0, rows.order, budget)
                ref following = self.cursors[self.heap[0]].term
                _charge_size(_sum(len(term.monomial.pairs), len(following.monomial.pairs)), budget)
                if term.monomial != following.monomial:
                    break
                term.coefficient = term.coefficient.add(following.coefficient, budget)
            _charge_size(term.coefficient.size(), budget)
            if not term.coefficient.is_zero():
                self.term = term^
                return True
        return False


struct _GroupOrder[C: CoefficientRing](_IndexedOrder):
    comptime Item = _ProductGroup[Self.C]

    @staticmethod
    @always_inline
    def compare(items: List[Self.Item], a: Int, b: Int, order: MonomialOrder, mut budget: Budget) raises PolyError -> Int:
        return order._compare(items[a].term.monomial, items[b].term.monomial, budget)


def _advance_column[C: CoefficientRing](mut group: _ProductGroup[C], index: Int, f: SparsePoly[C], left: List[SparsePoly[C]], indices: List[Int], rows: List[List[SparsePoly[C]]], column: Int, mut budget: Budget) raises PolyError -> Bool:
    if index == 0:
        if group.position == len(f.terms):
            return False
        group.term = f.terms[group.position]._copy(budget)
        group.position += 1
        return True
    var j = index - 1
    ref right = rows[indices[j]][column]
    var available = group.next(left[j], right, budget) if len(left[j].terms) <= len(right.terms) else group.next(right, left[j], budget)
    if available:
        group.term.coefficient = group.term.coefficient.neg(budget)
    return available


def _subtract_column[C: CoefficientRing](f: SparsePoly[C], left: List[SparsePoly[C]], indices: List[Int], rows: List[List[SparsePoly[C]]], column: Int, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    # All operands are prepared in the same admissible order. Borrow every
    # cofactor row: only product cursors and the final sum own coefficient data.
    var count = _sum(len(left), 1)
    _allocation[_ProductGroup[C]](count, budget)
    _allocation[Int](count, budget)
    var groups = List[_ProductGroup[C]](capacity=count)
    var heap = List[Int](capacity=count)
    for i in range(count):
        var group = _ProductGroup[C](f.ring)
        if i:
            ref right = rows[indices[i - 1]][column]
            ref a = left[i - 1]
            var prepared = _product_rows(a, right, budget) if len(a.terms) <= len(right.terms) else _product_rows(right, a, budget)
            swap(group.cursors, prepared[0])
            swap(group.heap, prepared[1])
        if _advance_column(group, i, f, left, indices, rows, column, budget):
            heap.append(i)
        groups.append(group^)
    if count <= 5:
        # A handful of product groups does not need an outer heap. Coalesce
        # every equal head before advancing, avoiding two heap repairs per
        # cancellation. Product rows still use their collision-folding heaps.
        _allocation[Bool](count, budget)
        var active = List[Bool](length=count, fill=False)
        for i in heap:
            active[i] = True
        var result = SparsePoly[C].zero(f.order, f.ring, budget)
        while True:
            var selected = -1
            for i in range(count):
                budget.charge(1)
                if active[i] and (selected < 0 or f.order._compare(groups[i].term.monomial, groups[selected].term.monomial, budget) < 0):
                    selected = i
            if selected < 0:
                return result^
            var term = Term[C](Monomial(), C.zero(f.ring))
            swap(term, groups[selected].term)
            for i in range(count):
                budget.charge(1)
                if not active[i]:
                    continue
                if i != selected:
                    _charge_size(_sum(len(term.monomial.pairs), len(groups[i].term.monomial.pairs)), budget)
                    if term.monomial != groups[i].term.monomial:
                        continue
                    term.coefficient = term.coefficient.add(groups[i].term.coefficient, budget)
                active[i] = _advance_column(groups[i], i, f, left, indices, rows, column, budget)
            _charge_size(term.coefficient.size(), budget)
            if not term.coefficient.is_zero():
                _allocation[Term[C]](1, budget)
                result.terms.append(term^)
    for i in range(len(heap) // 2 - 1, -1, -1):
        _sift_min[_GroupOrder[C]](heap, groups, i, f.order, budget)
    var result = SparsePoly[C].zero(f.order, f.ring, budget)
    while len(heap):
        var term = Term[C](Monomial(), C.zero(f.ring))
        swap(term, groups[heap[0]].term)
        while True:
            budget.charge(1)
            var i = heap[0]
            if not _advance_column(groups[i], i, f, left, indices, rows, column, budget):
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
        _charge_size(term.coefficient.size(), budget)
        if not term.coefficient.is_zero():
            _allocation[Term[C]](1, budget)
            result.terms.append(term^)
    return result^
