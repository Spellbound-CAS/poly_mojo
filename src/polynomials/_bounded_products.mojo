"""Balanced owned sums and degree-filtered sparse multiplication."""
from .coefficients import CoefficientRing
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term
from .errors import Budget, PolyError, _allocation, _sum, _product, _charge_size


struct _PolynomialSum[C: CoefficientRing](Movable):
    var bins: List[SparsePoly[Self.C]]
    var occupied: List[Bool]

    def __init__(out self):
        self.bins = List[SparsePoly[Self.C]]()
        self.occupied = List[Bool]()

    def add(mut self, var value: SparsePoly[Self.C], mut budget: Budget) raises PolyError:
        var level = 0
        while level < len(self.bins) and self.occupied[level]:
            budget.charge(1)
            value = self.bins[level].add(value, budget)
            self.bins[level].terms = List[Term[Self.C]]()
            self.occupied[level] = False
            level += 1
        if level == len(self.bins):
            _allocation[SparsePoly[Self.C]](1, budget)
            _allocation[Bool](1, budget)
            self.bins.append(value^)
            self.occupied.append(True)
        else:
            self.bins[level] = value^
            self.occupied[level] = True

    def finish(mut self, exemplar: SparsePoly[Self.C], mut budget: Budget) raises PolyError -> SparsePoly[Self.C]:
        var result = SparsePoly[Self.C].zero(exemplar.order, exemplar.ring, budget)
        for i in range(len(self.bins)):
            budget.charge(1)
            if self.occupied[i]:
                if result.is_zero():
                    swap(result, self.bins[i])
                else:
                    result = self.bins[i].add(result, budget)
                    self.bins[i].terms = List[Term[Self.C]]()
                self.occupied[i] = False
        return result^


def _parameter_degrees[C: CoefficientRing](f: SparsePoly[C], main: Int64, mut budget: Budget) raises PolyError -> List[Int]:
    _allocation[Int](len(f.terms), budget)
    var result = List[Int](capacity=len(f.terms))
    for term in f.terms:
        var degree = 0
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] != main:
                degree = _sum(degree, Int(UInt64(pair[1])))
        result.append(degree)
    return result^


def _truncated_product[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C], main: Int64, degree: Int, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    a._same(b, budget)
    var left = _parameter_degrees(a, main, budget)
    var right = _parameter_degrees(b, main, budget)
    var maximum_left = 0
    var maximum_right = 0
    for d in left:
        budget.charge(1)
        maximum_left = max(maximum_left, d)
    for d in right:
        budget.charge(1)
        maximum_right = max(maximum_right, d)
    if maximum_left <= degree and maximum_right <= degree - maximum_left:
        return a.mul(b, budget)
    var count = _product(len(left), len(right))
    if a.order.kind != MonomialOrder.KERNEL and 64 <= count <= 256:
        # For tiny, almost unfiltered products the existing row heap beats
        # sorting a raw batch. Its complete result has at most 256 terms.
        # Compact that owned result without copying retained coefficients.
        var retained = 0
        for d in left:
            for e in right:
                budget.charge(1)
                if d <= degree and e <= degree - d:
                    retained += 1
        if retained * 8 >= count * 7:
            var result = a.mul(b, budget)
            var degrees = _parameter_degrees(result, main, budget)
            var write = 0
            for read in range(len(result.terms)):
                budget.charge(1)
                if degrees[read] <= degree:
                    if write != read:
                        var kept = Term[C](Monomial(), C.zero(a.ring))
                        swap(kept, result.terms[read])
                        result.terms[write] = kept^
                    write += 1
            while len(result.terms) > write:
                budget.charge(1)
                _ = result.terms.pop()
            return result^
    var sum = _PolynomialSum[C]()
    var batch = List[Term[C]]()
    for i in range(len(left)):
        budget.charge(1)
        if left[i] > degree:
            continue
        for j in range(len(right)):
            budget.charge(1)
            if right[j] > degree - left[i]:
                continue
            var coefficient = a.terms[i].coefficient.mul(b.terms[j].coefficient, budget)
            _charge_size(coefficient.size(), budget)
            if coefficient.is_zero():
                continue
            _allocation[Term[C]](1, budget)
            batch.append(Term[C](a.terms[i].monomial.mul(b.terms[j].monomial, budget), coefficient^))
            if len(batch) == 256:
                sum.add(SparsePoly[C].from_terms(batch^, a.order, a.ring, budget), budget)
                batch = List[Term[C]]()
    if len(batch):
        sum.add(SparsePoly[C].from_terms(batch^, a.order, a.ring, budget), budget)
    return sum.finish(a, budget)
