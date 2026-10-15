"""Polynomial calculus and substitutions with explicit numerical contexts."""
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum, _product
from .numeric_domains import NumericDomain, _numeric_cost, _context_cost
from .numeric_dense import NumericDensePoly
from .numeric_sparse import NumericSparsePoly, NumericTerm
from .monomial import Monomial
from ._assignment import _AssignmentPlan
from ._toolkit_numeric import _NumericOps, _validate
from ._toolkit_ops import _integral_values, _compose_values, _jet_values


def numeric_integrate_dense[D: NumericDomain](poly: NumericDensePoly[D], context: D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    """Integrate coefficient-wise using an explicit numerical context.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The polynomial.
        context: The rounding or enclosing arithmetic settings.
        budget: The shared work allowance.

    Returns:
        A rounded or enclosing antiderivative with zero constant.

    Raises:
        PolyError: DOMAIN for invalid values or contexts, LIMIT or CAPACITY on work limits.
    """
    _validate[D](poly.coefficients, context, budget)
    var values = _integral_values[_NumericOps[D]](_view(poly.coefficients), context, budget)
    return NumericDensePoly[D].from_coefficients(poly.variable, values^, budget)


def numeric_definite_integral[D: NumericDomain](poly: NumericDensePoly[D], low: D.Scalar, high: D.Scalar, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    """Evaluate a numerical antiderivative difference at explicit endpoints.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The polynomial.
        low: The lower endpoint.
        high: The upper endpoint.
        context: The explicit arithmetic settings.
        budget: The shared work allowance.

    Returns:
        A rounded value or uniform enclosure of the integral difference.

    Raises:
        PolyError: DOMAIN for invalid values or contexts, LIMIT or CAPACITY on work limits.
    """
    var integral = numeric_integrate_dense(poly, context, budget)
    return _NumericOps[D].apply(integral.evaluate(high, context, budget), integral.evaluate(low, context, budget), 1, context, budget)


def numeric_compose_dense[D: NumericDomain](outer: NumericDensePoly[D], inner: NumericDensePoly[D], context: D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    """Compose numerical polynomials by a fixed Horner schedule.

    Parameters:
        D: The numerical domain.

    Args:
        outer: The outer polynomial.
        inner: The replacement in the same variable.
        context: The explicit arithmetic settings.
        budget: The shared work allowance.

    Returns:
        The rounded composition or coefficient enclosures of every composition.

    Raises:
        PolyError: DOMAIN for incompatible values, LIMIT or CAPACITY on work limits.
    """
    outer._same(inner)
    _validate[D](outer.coefficients, context, budget)
    _validate[D](inner.coefficients, context, budget)
    var count = _sum(_product(max(0, outer.degree_bound()), max(0, inner.degree_bound())), 1)
    var values = _compose_values[_NumericOps[D]](_view(outer.coefficients), _view(inner.coefficients), count, context, budget)
    return NumericDensePoly[D].from_coefficients(outer.variable, values^, budget)


def numeric_affine_transform[D: NumericDomain](poly: NumericDensePoly[D], scale: D.Scalar, shift: D.Scalar, context: D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    """Compose with scale*x + shift using the requested precision.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The polynomial.
        scale: The finite scalar or enclosure multiplying x.
        shift: The finite translation or enclosure.
        context: The explicit arithmetic settings.
        budget: The shared work allowance.

    Returns:
        The rounded or enclosing transformed polynomial.

    Raises:
        PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
    """
    var inner = NumericDensePoly[D].from_coefficients(poly.variable, [_NumericOps[D].copy(shift, budget), _NumericOps[D].copy(scale, budget)], budget)
    return numeric_compose_dense(poly, inner, context, budget)


def numeric_jet_evaluate[D: NumericDomain](poly: NumericDensePoly[D], point: D.Scalar, order: Int, context: D.Context, mut budget: Budget) raises PolyError -> List[D.Scalar]:
    """Compute numerical Taylor coefficients by generalized Horner evaluation.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The polynomial.
        point: The expansion point or enclosure.
        order: The nonnegative number of Taylor coefficients.
        context: The explicit arithmetic settings.
        budget: The shared work allowance.

    Returns:
        Rounded coefficients or uniform enclosures, including requested zero padding.

    Raises:
        PolyError: DOMAIN for invalid values or negative order, LIMIT or CAPACITY on work limits.
    """
    _validate[D](poly.coefficients, context, budget)
    _charge_size(_numeric_cost[D](point), budget)
    return _jet_values[_NumericOps[D]](_view(poly.coefficients), point, order, context, budget)


def numeric_taylor_shift[D: NumericDomain](poly: NumericDensePoly[D], point: D.Scalar, context: D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    """Translate the variable using numerical Taylor coefficients.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The polynomial.
        point: The shift or its enclosure.
        context: The explicit arithmetic settings.
        budget: The shared work allowance.

    Returns:
        A rounded or enclosing polynomial for p(x + point).

    Raises:
        PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
    """
    var values = numeric_jet_evaluate(poly, point, len(poly.coefficients), context, budget)
    return NumericDensePoly[D].from_coefficients(poly.variable, values^, budget)


def numeric_integrate_sparse[D: NumericDomain](poly: NumericSparsePoly[D], variable: Int64, context: D.Context, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    """Integrate a sparse numerical polynomial in one variable.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The input polynomial.
        variable: The integration variable allowed by its order.
        context: The explicit arithmetic settings.
        budget: The shared work allowance.

    Returns:
        A rounded or enclosing antiderivative with zero integration constant.

    Raises:
        PolyError: DOMAIN for invalid values or variables, LIMIT or CAPACITY on work limits.
    """
    var unit = Monomial([(variable, UInt32(1))], budget)
    poly.order.validate(unit, budget)
    _allocation[NumericTerm[D]](len(poly.terms), budget)
    var terms = List[NumericTerm[D]](capacity=len(poly.terms))
    for term in poly.terms:
        poly.order.validate(term.monomial, budget)
        var exponent = term.monomial.exponent(variable, budget)
        if exponent == UInt32.MAX:
            raise PolyError(PolyErrorKind.CAPACITY, "integration exponent exceeds UInt32")
        var coefficient = _NumericOps[D].apply(term.coefficient, _NumericOps[D].integer(Int(exponent) + 1, context, budget), 3, context, budget)
        terms.append(NumericTerm[D](term.monomial.mul(unit, budget), coefficient^))
    return NumericSparsePoly[D].from_terms(terms^, poly.order, context, budget)


def numeric_substitute_many[D: NumericDomain](poly: NumericSparsePoly[D], assignments: List[Tuple[Int64, NumericSparsePoly[D]]], context: D.Context, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    """Replace numerical variables simultaneously in a fixed term order.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The input polynomial.
        assignments: Distinct variable IDs and replacement polynomials.
        context: The explicit arithmetic settings.
        budget: The shared work allowance.

    Returns:
        A rounded composition or an enclosure of every represented composition.

    Raises:
        PolyError: DOMAIN for duplicates or incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    _ = _context_cost[D](context)
    var plan = _AssignmentPlan()
    for i in range(len(assignments)):
        poly._same(assignments[i][1], budget)
        poly.order.validate(Monomial([(assignments[i][0], UInt32(1))], budget), budget)
        for term in assignments[i][1].terms:
            poly.order.validate(term.monomial, budget)
            _charge_size(_numeric_cost[D](term.coefficient), budget)
        plan.assign(assignments[i][0], i, budget)
    _charge_size(len(poly.order.variables), budget)
    var result = NumericSparsePoly[D](poly.order.copy())
    var powers = List[NumericSparsePoly[D]]()
    for term in poly.terms:
        poly.order.validate(term.monomial, budget)
        var product = NumericSparsePoly[D].from_terms([NumericTerm[D](Monomial(), _NumericOps[D].copy(term.coefficient, budget))], poly.order, context, budget)
        for pair in term.monomial.pairs:
            var slot = plan.power(pair, budget)
            if slot == len(powers):
                var base = NumericSparsePoly[D].from_terms([NumericTerm[D](Monomial([(pair[0], UInt32(1))], budget), _NumericOps[D].integer(1, context, budget))], poly.order, context, budget)
                var assigned = plan.assignment(pair[0], budget)
                if assigned >= 0:
                    ref replacement = assignments[assigned][1]
                    _allocation[NumericTerm[D]](len(replacement.terms), budget)
                    for t in replacement.terms:
                        _charge_size(_sum(len(t.monomial.pairs), _numeric_cost[D](t.coefficient)), budget)
                    base = replacement.copy()
                var power = base.pow(UInt64(pair[1]), context, budget)
                _allocation[NumericSparsePoly[D]](1, budget)
                powers.append(power^)
            product = product.mul(powers[slot], context, budget)
        result = result.add(product, context, budget)
    return result^
