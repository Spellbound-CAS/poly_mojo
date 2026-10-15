"""Lossless numerical copies with explicit variable maps and canonical array axes."""
from polynomials import Budget, PolyError, PolyErrorKind, NumericDomain, NumericDensePoly, NumericSparsePoly, NumericTerm, Monomial, MonomialOrder, NumericBox, numeric_bounds
from ._common import _allocation
from ._interop import _mapped_order, _remap_pairs


def _zero(box: NumericBox) -> Bool:
    return box.real_low.is_zero() and box.real_high.is_zero() and box.imag_low.is_zero() and box.imag_high.is_zero()


def _dense[D: NumericDomain](source: NumericDensePoly[D], variable: Int64, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    _allocation[D.Scalar](len(source.coefficients), budget)
    var values = List[D.Scalar](capacity=len(source.coefficients))
    for i in range(len(source.coefficients)):
        var bounds = numeric_bounds[D](source.coefficients[i], budget)
        if i == len(source.coefficients) - 1 and _zero(bounds):
            raise PolyError(PolyErrorKind.DOMAIN, "noncanonical trailing zero coefficient")
        values.append(source.coefficients[i].copy())
    return NumericDensePoly[D].from_coefficients(variable, values^, budget)


def _sift[D: NumericDomain](mut terms: List[NumericTerm[D]], order: MonomialOrder, start: Int, end: Int, mut budget: Budget) raises PolyError:
    var parent = start
    while parent < end // 2:
        var child = parent * 2 + 1
        if child + 1 < end and order.compare(terms[child].monomial, terms[child + 1].monomial, budget) < 0:
            child += 1
        if order.compare(terms[parent].monomial, terms[child].monomial, budget) >= 0: break
        terms.swap_elements(parent, child)
        parent = child


def _sparse[D: NumericDomain](source: NumericSparsePoly[D], mapping: List[Int64], into_core: Bool, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    var order = _mapped_order(source.order, mapping, into_core, budget)
    _allocation[NumericTerm[D]](len(source.terms), budget)
    var terms = List[NumericTerm[D]](capacity=len(source.terms))
    for i in range(len(source.terms)):
        ref term = source.terms[i]
        source.order.validate(term.monomial, budget)
        if i and source.order.compare(source.terms[i - 1].monomial, term.monomial, budget) >= 0:
            raise PolyError(PolyErrorKind.DOMAIN, "noncanonical numerical term order")
        if _zero(numeric_bounds[D](term.coefficient, budget)):
            raise PolyError(PolyErrorKind.DOMAIN, "zero numerical sparse coefficient")
        terms.append(NumericTerm[D](Monomial(_remap_pairs(term.monomial, mapping, into_core, budget), budget), term.coefficient.copy()))
    for i in range(len(terms) // 2 - 1, -1, -1):
        _sift(terms, order, i, len(terms), budget)
    for end in range(len(terms) - 1, 0, -1):
        terms.swap_elements(0, end)
        _sift(terms, order, 0, end, budget)
    var result = NumericSparsePoly[D](order^)
    result.terms = terms^
    return result^
