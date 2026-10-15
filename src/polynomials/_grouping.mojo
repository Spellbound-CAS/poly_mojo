"""Stable occupied coefficient groups and one-variable monomial partitions."""
from .errors import Budget, PolyError, _allocation, _charge_size
from .coefficients import CoefficientRing
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term, _sorted_indices
from ._ordering import _IndexedOrder, _stable_indices


struct _ExponentOrder(_IndexedOrder):
    comptime Item = UInt32

    @staticmethod
    def compare(items: List[UInt32], a: Int, b: Int, order: MonomialOrder, mut budget: Budget) raises PolyError -> Int:
        budget.charge(1)
        return Int(items[a] > items[b]) - Int(items[a] < items[b])


struct _CoefficientGroups(Movable):
    var variable: Int64
    var exponents: List[UInt32]
    var indices: List[Int]
    var starts: List[Int]

    def __init__[C: CoefficientRing](out self, f: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError:
        self.variable = variable
        _allocation[UInt32](len(f.terms), budget)
        self.exponents = List[UInt32](capacity=len(f.terms))
        for term in f.terms:
            self.exponents.append(term.monomial.exponent(variable, budget))
        self.indices = _stable_indices[_ExponentOrder](self.exponents, f.order, budget)
        self.starts = List[Int]()
        for i in range(len(self.indices)):
            budget.charge(1)
            if i == 0 or self.exponents[self.indices[i]] != self.exponents[self.indices[i - 1]]:
                _allocation[Int](1, budget)
                self.starts.append(i)
        _allocation[Int](1, budget)
        self.starts.append(len(self.indices))

    def exponent(self, group: Int) -> UInt32:
        return self.exponents[self.indices[self.starts[group]]]

    def coefficient[C: CoefficientRing](self, f: SparsePoly[C], group: Int, mut budget: Budget) raises PolyError -> SparsePoly[C]:
        var begin = self.starts[group]
        var end = self.starts[group + 1]
        _allocation[Term[C]](end - begin, budget)
        var terms = List[Term[C]](capacity=end - begin)
        for i in range(begin, end):
            budget.charge(1)
            ref term = f.terms[self.indices[i]]
            _charge_size(term.coefficient.size(), budget)
            terms.append(Term[C](_without(term.monomial, self.variable, budget), term.coefficient.copy()))
        return SparsePoly[C].from_terms(terms^, f.order, f.ring, budget)


def _separate(m: Monomial, variable: Int64, mut budget: Budget) raises PolyError -> Tuple[UInt32, Monomial]:
    _allocation[Tuple[Int64, UInt32]](len(m.pairs), budget)
    var rest = Monomial()
    rest.pairs.reserve(len(m.pairs))
    var exponent = UInt32(0)
    for pair in m.pairs:
        budget.charge(1)
        if pair[0] == variable:
            exponent = pair[1]
        else:
            rest.pairs.append(pair)
    return (exponent, rest^)


def _without(m: Monomial, variable: Int64, mut budget: Budget) raises PolyError -> Monomial:
    var parts = _separate(m, variable, budget)
    var result = Monomial()
    swap(result, parts[1])
    return result^


def _index_groups[C: CoefficientRing](keys: List[Term[C]], mut budget: Budget, first: Bool = True) raises PolyError -> List[List[Int]]:
    var kernel = MonomialOrder()
    # Coefficients are payload only; their contexts do not participate in sorting.
    if len(keys) == 0:
        return List[List[Int]]()
    var indices = _sorted_indices(keys, kernel, keys[0].coefficient.context(), budget)
    _allocation[Int](len(keys) if first else 0, budget)
    var mapping = List[Int](length=len(keys) if first else 0, fill=-1)
    var groups = List[List[Int]]()
    var previous = -1
    for index in indices:
        budget.charge(1)
        if previous < 0 or kernel._compare(keys[previous].monomial, keys[index].monomial, budget) != 0:
            _allocation[List[Int]](1, budget)
            groups.append(List[Int]())
        var group = len(groups) - 1
        _allocation[Int](1, budget)
        groups[group].append(index)
        if first:
            mapping[index] = group
        previous = index
    if not first:
        return groups^
    _allocation[Bool](len(groups), budget)
    var seen = List[Bool](length=len(groups), fill=False)
    _allocation[List[Int]](len(groups), budget)
    var ordered = List[List[Int]](capacity=len(groups))
    for group in mapping:
        budget.charge(1)
        if not seen[group]:
            seen[group] = True
            var indices = List[Int]()
            swap(indices, groups[group])
            ordered.append(indices^)
    return ordered^


def _partition[C: CoefficientRing](var keys: List[Term[C]], var values: List[Monomial], exemplar: SparsePoly[C], mut budget: Budget, first: Bool = True) raises PolyError -> Tuple[List[Monomial], List[SparsePoly[C]]]:
    var groups = _index_groups(keys, budget, first)
    _allocation[Monomial](len(groups), budget)
    _allocation[SparsePoly[C]](len(groups), budget)
    var monomials = List[Monomial](capacity=len(groups))
    var coefficients = List[SparsePoly[C]](capacity=len(groups))
    for group in groups:
        budget.charge(1)
        var key = Monomial()
        swap(key, keys[group[0]].monomial)
        monomials.append(key^)
        _allocation[Term[C]](len(group), budget)
        var terms = List[Term[C]](capacity=len(group))
        for index in group:
            budget.charge(1)
            var monomial = Monomial()
            var coefficient = C.zero(exemplar.ring)
            swap(monomial, values[index])
            swap(coefficient, keys[index].coefficient)
            terms.append(Term[C](monomial^, coefficient^))
        coefficients.append(SparsePoly[C].from_terms(terms^, exemplar.order, exemplar.ring, budget))
    return (monomials^, coefficients^)


def _coefficient_groups[C: CoefficientRing](f: SparsePoly[C], variable: Int64, univariate: Bool, mut budget: Budget, first: Bool = True) raises PolyError -> Tuple[List[Monomial], List[SparsePoly[C]]]:
    _allocation[Term[C]](len(f.terms), budget)
    _allocation[Monomial](len(f.terms), budget)
    var keys = List[Term[C]](capacity=len(f.terms))
    var values = List[Monomial](capacity=len(f.terms))
    for term in f.terms:
        budget.charge(1)
        var parts = _separate(term.monomial, variable, budget)
        var single = Monomial([(variable, parts[0])], budget) if parts[0] else Monomial()
        var rest = Monomial()
        swap(rest, parts[1])
        _charge_size(term.coefficient.size(), budget)
        if univariate:
            keys.append(Term[C](rest^, term.coefficient.copy()))
            values.append(single^)
        else:
            keys.append(Term[C](single^, term.coefficient.copy()))
            values.append(rest^)
    return _partition(keys^, values^, f, budget, first)
