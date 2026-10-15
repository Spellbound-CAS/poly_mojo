"""Internal views of sparse polynomials as univariates over smaller rings."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _charge_size
from .coefficients import CoefficientRing, CoefficientField, ZZ
from apn_mojo import Integer
from .monomial import Monomial
from .sparse import SparsePoly, Term
from ._grouping import _coefficient_groups, _without


def _constant[C: CoefficientRing](p: SparsePoly[C], value: C, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    return SparsePoly[C].constant(value, p.order, budget)


def _one[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    return _constant(p, C.one(p.ring), budget)


def _zero[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    return SparsePoly[C].zero(p.order, p.ring, budget)


def _split[C: CoefficientRing](p: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> List[SparsePoly[C]]:
    if p.is_zero():
        return List[SparsePoly[C]]()
    var count = _sum(p.degree_in(variable, budget), 1)
    var grouped = _coefficient_groups(p, variable, False, budget)
    _allocation[SparsePoly[C]](count, budget)
    var result = List[SparsePoly[C]](capacity=count)
    for _ in range(count):
        budget.charge(1)
        result.append(_zero(p, budget))
    for i in range(len(grouped[0])):
        budget.charge(1)
        var exponent = Int(UInt64(grouped[0][i].exponent(variable, budget)))
        swap(result[exponent], grouped[1][i])
    return result^


def _join[C: CoefficientRing](values: List[SparsePoly[C]], variable: Int64, exemplar: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var count = 0
    for value in values:
        budget.charge(1)
        exemplar._same(value, budget)
        count = _sum(count, len(value.terms))
    _allocation[Term[C]](count, budget)
    var terms = List[Term[C]](capacity=count)
    for i in range(len(values)):
        budget.charge(1)
        if UInt64(i) > UInt64(UInt32.MAX):
            raise PolyError(PolyErrorKind.CAPACITY, "univariate degree exceeds UInt32")
        if values[i].is_zero():
            continue
        var shift = Monomial([(variable, UInt32(i))], budget) if i else Monomial()
        for term in values[i].terms:
            budget.charge(1)
            _charge_size(term.coefficient.size(), budget)
            terms.append(Term[C](term.monomial.mul(shift, budget), term.coefficient.copy()))
    return SparsePoly[C].from_terms(terms^, exemplar.order, exemplar.ring, budget)


def _shift[C: CoefficientRing](p: SparsePoly[C], variable: Int64, degree: UInt32, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    if not degree or p.is_zero():
        return p._copy(budget)
    var monomial = Monomial([(variable, degree)], budget)
    _allocation[Term[C]](len(p.terms), budget)
    var terms = List[Term[C]](capacity=len(p.terms))
    for term in p.terms:
        budget.charge(1)
        _charge_size(term.coefficient.size(), budget)
        terms.append(Term[C](term.monomial.mul(monomial, budget), term.coefficient.copy()))
    return SparsePoly[C].from_terms(terms^, p.order, p.ring, budget)


def _lc[C: CoefficientRing](p: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var degree = p.degree_in(variable, budget)
    if degree < 0:
        return _zero(p, budget)
    return _coefficient(p, variable, UInt32(degree), budget)


def _coefficient[C: CoefficientRing](p: SparsePoly[C], variable: Int64, exponent: UInt32, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var terms = List[Term[C]]()
    for term in p.terms:
        budget.charge(1)
        if term.monomial.exponent(variable, budget) != exponent:
            continue
        _allocation[Term[C]](1, budget)
        _charge_size(term.coefficient.size(), budget)
        terms.append(Term[C](_without(term.monomial, variable, budget), term.coefficient.copy()))
    return SparsePoly[C].from_terms(terms^, p.order, p.ring, budget)


def _variables[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], mut budget: Budget) raises PolyError -> List[Int64]:
    var variables = f.variables(budget)
    for value in g.variables(budget):
        budget.charge(1)
        var position = 0
        while position < len(variables) and variables[position] < value:
            budget.charge(1)
            position += 1
        if position == len(variables) or variables[position] != value:
            _ = _sum(len(variables), 1)
            _allocation[Int64](1, budget)
            variables.append(value)
            var j = len(variables) - 1
            while j > position:
                budget.charge(1)
                variables[j] = variables[j - 1]
                j -= 1
            variables[position] = value
    return variables^


def _is_unit[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    budget.charge(1)
    if len(p.terms) != 1 or len(p.terms[0].monomial.pairs):
        return False
    comptime if C == ZZ:
        _charge_size(p.terms[0].coefficient.size(), budget)
        ref value = rebind[ZZ](p.terms[0].coefficient).value
        return value == Integer(1) or value == Integer(-1)
    elif conforms_to(C, CoefficientField):
        return not p.terms[0].coefficient.is_zero()
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "unit test requires ZZ or a coefficient field")
