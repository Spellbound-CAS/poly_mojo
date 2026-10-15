"""Bounded descending streams for exact sparse reduction tails."""
from std.sys import size_of
from .coefficients import CoefficientField
from .errors import Budget, PolyError, _allocation, _sum, _product, _charge_size
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term
from ._ordering import _IndexedOrder, _sift_min
from ._sparse_ops import _from_descending


@fieldwise_init
struct _Tail[C: CoefficientField](Copyable, Movable):
    var divisor: Int
    var position: Int
    var factor: Term[Self.C]
    var term: Term[Self.C]


struct _TailOrder[C: CoefficientField](_IndexedOrder):
    comptime Item = _Tail[Self.C]

    @staticmethod
    @always_inline
    def compare(items: List[Self.Item], a: Int, b: Int, order: MonomialOrder, mut budget: Budget) raises PolyError -> Int:
        return -order._compare(items[a].term.monomial, items[b].term.monomial, budget)


struct _ReductionStream[C: CoefficientField](Movable):
    var base: SparsePoly[Self.C]
    var cursors: List[_Tail[Self.C]]
    var heap: List[Int]
    var free: List[Int]

    def __init__(out self, var f: SparsePoly[Self.C], mut budget: Budget) raises PolyError:
        self.base = f^
        self.cursors = List[_Tail[Self.C]]()
        self.heap = List[Int]()
        self.free = List[Int]()
        self._base(budget)

    def _base(mut self, mut budget: Budget) raises PolyError:
        if len(self.base.terms):
            _allocation[_Tail[Self.C]](1, budget)
            _allocation[Int](1, budget)
            self.cursors.append(_Tail[Self.C](-1, 0, Term[Self.C](Monomial(), Self.C.zero(self.base.ring)), self.base.terms.pop()))
            self.heap.append(0)

    def _advance(mut self, i: Int, divisors: List[SparsePoly[Self.C]], mut budget: Budget) raises PolyError -> Bool:
        if self.cursors[i].divisor < 0:
            if not len(self.base.terms):
                return False
            self.cursors[i].term = self.base.terms.pop()
            return True
        self.cursors[i].position -= 1
        var position = self.cursors[i].position
        if position < 0:
            return False
        var divisor = self.cursors[i].divisor
        ref original = divisors[divisor].terms[position]
        ref factor = self.cursors[i].factor
        self.cursors[i].term = Term[Self.C](original.monomial.mul(factor.monomial, budget), original.coefficient.mul(factor.coefficient, budget))
        return True

    def next(mut self, divisors: List[SparsePoly[Self.C]], mut budget: Budget) raises PolyError -> Optional[Term[Self.C]]:
        while len(self.heap):
            budget.charge(1)
            var term = Term[Self.C](Monomial(), Self.C.zero(self.base.ring))
            swap(term, self.cursors[self.heap[0]].term)
            while True:
                var index = self.heap[0]
                if not self._advance(index, divisors, budget):
                    var last = self.heap.pop()
                    if len(self.heap):
                        self.heap[0] = last
                    _allocation[Int](1, budget)
                    self.free.append(index)
                    self.cursors[index].factor = Term[Self.C](Monomial(), Self.C.zero(self.base.ring))
                    self.cursors[index].term = Term[Self.C](Monomial(), Self.C.zero(self.base.ring))
                if not len(self.heap):
                    break
                _sift_min[_TailOrder[Self.C]](self.heap, self.cursors, 0, self.base.order, budget)
                ref following = self.cursors[self.heap[0]].term
                _charge_size(_sum(len(term.monomial.pairs), len(following.monomial.pairs)), budget)
                if term.monomial != following.monomial:
                    break
                term.coefficient = term.coefficient.add(following.coefficient, budget)
            _charge_size(term.coefficient.size(), budget)
            if not term.coefficient.is_zero():
                return Optional(term^)
        return None

    def append(mut self, divisor: Int, factor: Term[Self.C], divisors: List[SparsePoly[Self.C]], mut budget: Budget) raises PolyError:
        var position = len(divisors[divisor].terms) - 2
        if position < 0:
            return
        # Bound the live row frontier independently of the number of quotient
        # terms. Merge outstanding rows and recycle slots before admitting more.
        if len(self.heap) == 64:
            var terms = List[Term[Self.C]]()
            while True:
                var next = self.next(divisors, budget)
                if not next:
                    break
                _ = _product(_sum(len(terms), 1), size_of[Term[Self.C]]())
                _allocation[Term[Self.C]](1, budget)
                terms.append(next.take())
            var base = _from_descending(terms^, self.base, budget)
            self.cursors = List[_Tail[Self.C]]()
            self.heap = List[Int]()
            self.free = List[Int]()
            self.base = base^
            self._base(budget)
        var negative = Term[Self.C](Monomial(), factor.coefficient.neg(budget))
        _charge_size(len(factor.monomial.pairs), budget)
        negative.monomial = factor.monomial.copy()
        ref original = divisors[divisor].terms[position]
        var current = Term[Self.C](original.monomial.mul(negative.monomial, budget), original.coefficient.mul(negative.coefficient, budget))
        var cursor = _Tail[Self.C](divisor, position, negative^, current^)
        var index: Int
        if len(self.free):
            index = self.free.pop()
            self.cursors[index] = cursor^
        else:
            index = len(self.cursors)
            _allocation[_Tail[Self.C]](_sum(index, 1), budget)
            self.cursors.append(cursor^)
        _allocation[Int](_sum(len(self.heap), 1), budget)
        self.heap.append(index)
        var child = len(self.heap) - 1
        while child:
            budget.charge(1)
            var parent = (child - 1) // 2
            if _TailOrder[Self.C].compare(self.cursors, index, self.heap[parent], self.base.order, budget) >= 0:
                break
            self.heap[child] = self.heap[parent]
            child = parent
        self.heap[child] = index
