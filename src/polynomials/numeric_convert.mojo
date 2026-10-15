"""Explicit bridges between exact and numerical polynomial representations."""
from apn_mojo import Rational
from .coefficients import CoefficientRing, ZZ, QQ
from .exact_complex import QQi
from .dense import DensePoly
from .sparse import SparsePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .numeric_domains import NumericDomain, _context_cost
from .numeric_dense import NumericDensePoly
from .numeric_sparse import NumericSparsePoly, NumericTerm
from .numeric_bounds import NumericBox, numeric_bounds, numeric_value


def _require_exact_embedding[C: CoefficientRing]() raises PolyError:
    comptime if C != QQ and C != ZZ and C != QQi:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "only ZZ, QQ, and QQi have an exact numerical embedding")


def _exact_box[C: CoefficientRing](value: C, mut budget: Budget) raises PolyError -> NumericBox:
    _charge_size(value.size(), budget)
    var q: QQ
    comptime if C == QQ:
        q = rebind[QQ](value).copy()
    elif C == ZZ:
        q = QQ(Rational(rebind[ZZ](value).value))
    elif C == QQi:
        ref z = rebind[QQi](value).value
        return NumericBox.point(QQ(z.real()), QQ(z.imag()), budget)
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "only ZZ, QQ, and QQi have an exact numerical embedding")
    return NumericBox.point(q, QQ(Rational(0)), budget)


def convert_numeric_dense[E: NumericDomain, D: NumericDomain](poly: NumericDensePoly[D], context: E.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[E]:
    """Convert numerical coefficients with enclosing or certainly rounded semantics.

    Parameters:
        E: The destination domain.
        D: The source domain.

    Args:
        poly: The source polynomial.
        context: Explicit destination settings.
        budget: The shared work allowance.

    Returns:
        The converted polynomial in the same variable.

    Raises:
        PolyError: DOMAIN for uncertain point conversion, LIMIT or CAPACITY on work limits.
    """
    _ = _context_cost[E](context)
    _allocation[E.Scalar](len(poly.coefficients), budget)
    var values = List[E.Scalar](capacity=len(poly.coefficients))
    for value in poly.coefficients:
        values.append(numeric_value[E](numeric_bounds[D](value, budget), context, budget))
    return NumericDensePoly[E].from_coefficients(poly.variable, values^, budget)


def convert_numeric_sparse[E: NumericDomain, D: NumericDomain](poly: NumericSparsePoly[D], context: E.Context, mut budget: Budget) raises PolyError -> NumericSparsePoly[E]:
    """Convert sparse numerical coefficients explicitly between domains.

    Parameters:
        E: The destination domain.
        D: The source domain.

    Args:
        poly: The source polynomial.
        context: Explicit destination settings.
        budget: The shared work allowance.

    Returns:
        The converted polynomial in the same monomial order.

    Raises:
        PolyError: DOMAIN for uncertain point conversion, LIMIT or CAPACITY on work limits.
    """
    _allocation[NumericTerm[E]](len(poly.terms), budget)
    var terms = List[NumericTerm[E]](capacity=len(poly.terms))
    for term in poly.terms:
        _charge_size(len(term.monomial.pairs), budget)
        terms.append(NumericTerm[E](term.monomial.copy(), numeric_value[E](numeric_bounds[D](term.coefficient, budget), context, budget)))
    return NumericSparsePoly[E].from_terms(terms^, poly.order, context, budget)


def numeric_from_dense[D: NumericDomain, C: CoefficientRing](poly: DensePoly[C], context: D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    """Convert an exact integer, rational, or Gaussian rational polynomial directly into a numerical domain.

    Parameters:
        D: The destination numerical domain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The exact source polynomial.
        context: Explicit destination settings.
        budget: The shared work allowance.

    Returns:
        A rounded point polynomial or rigorous coefficient enclosure.

    Raises:
        PolyError: UNSUPPORTED for other exact domains, DOMAIN for invalid conversion, LIMIT or CAPACITY on work limits.
    """
    _require_exact_embedding[C]()
    poly.check_invariants(budget)
    _ = _context_cost[D](context)
    _allocation[D.Scalar](len(poly.coefficients), budget)
    var values = List[D.Scalar](capacity=len(poly.coefficients))
    for index in range(len(poly.coefficients)):
        ref value = poly.coefficients[index]
        values.append(numeric_value[D](_exact_box(value, budget), context, budget))
    return NumericDensePoly[D].from_coefficients(poly.variable, values^, budget)


def numeric_from_sparse[D: NumericDomain, C: CoefficientRing](poly: SparsePoly[C], context: D.Context, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    """Convert exact sparse integer, rational, or Gaussian rational coefficients to numerical values.

    Parameters:
        D: The destination numerical domain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The exact source polynomial.
        context: Explicit destination settings.
        budget: The shared work allowance.

    Returns:
        A numerical polynomial preserving the monomial order.

    Raises:
        PolyError: UNSUPPORTED for other exact domains, DOMAIN for invalid conversion, LIMIT or CAPACITY on work limits.
    """
    _require_exact_embedding[C]()
    poly.check_invariants(budget)
    _allocation[NumericTerm[D]](len(poly.terms), budget)
    var terms = List[NumericTerm[D]](capacity=len(poly.terms))
    for term in poly.terms:
        _charge_size(len(term.monomial.pairs), budget)
        terms.append(NumericTerm[D](term.monomial.copy(), numeric_value[D](_exact_box(term.coefficient, budget), context, budget)))
    return NumericSparsePoly[D].from_terms(terms^, poly.order, context, budget)
