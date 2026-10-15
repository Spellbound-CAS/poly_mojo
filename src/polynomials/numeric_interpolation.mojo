"""Rounded and enclosing interpolation with explicit numerical contexts."""
from .coefficient_buffer import _view
from .errors import Budget, PolyError, _allocation
from .numeric_domains import NumericDomain
from .numeric_dense import NumericDensePoly
from ._toolkit_numeric import _NumericOps, _validate
from ._toolkit_ops import _values, _divided_differences, _newton_values, _hermite_values


def numeric_multipoint_evaluate[D: NumericDomain](poly: NumericDensePoly[D], points: List[D.Scalar], context: D.Context, mut budget: Budget) raises PolyError -> List[D.Scalar]:
    """Evaluate numerical samples using a fixed Horner schedule per point.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The polynomial or coefficient family.
        points: The evaluation points or enclosures.
        context: Explicit rounding or enclosure settings.
        budget: The shared work allowance.

    Returns:
        Rounded values or uniform enclosures in point order.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    _validate[D](poly.coefficients, context, budget)
    _validate[D](points, context, budget)
    _allocation[D.Scalar](len(points), budget)
    var values = List[D.Scalar](capacity=len(points))
    for point in points:
        values.append(poly.evaluate(point, context, budget))
    return values^


def numeric_interpolate[D: NumericDomain](variable: Int64, points: List[D.Scalar], values: List[D.Scalar], context: D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    """Interpolate numerical data by divided differences.

    Parameters:
        D: The numerical domain.

    Args:
        variable: The output polynomial variable.
        points: Distinct nodes, or pairwise separated node enclosures.
        values: Corresponding values or enclosures.
        context: Explicit rounding or enclosure settings.
        budget: The shared work allowance.

    Returns:
        A rounded interpolant or coefficient enclosures for every data selection.

    Raises:
        PolyError: DOMAIN for invalid data or node differences containing zero, LIMIT or CAPACITY on work limits.
    """
    _validate[D](points, context, budget)
    _validate[D](values, context, budget)
    var differences = _divided_differences[_NumericOps[D]](_view(points), _view(values), context, budget)
    var nodes = _values[_NumericOps[D]](_view(points), max(0, len(points) - 1), context, budget)
    var coefficients = _newton_values[_NumericOps[D]](_view(differences), _view(nodes), context, budget)
    return NumericDensePoly[D].from_coefficients(variable, coefficients^, budget)


def numeric_hermite_interpolate[D: NumericDomain](variable: Int64, points: List[D.Scalar], jets: List[List[D.Scalar]], context: D.Context, mut budget: Budget, *, derivatives: Bool = True) raises PolyError -> NumericDensePoly[D]:
    """Interpolate numerical values and consecutive derivative constraints.

    Parameters:
        D: The numerical domain.

    Args:
        variable: The output variable.
        points: Distinct nodes or separated node enclosures.
        jets: Nonempty consecutive constraints beginning with the value.
        context: Explicit rounding or enclosure settings.
        budget: The shared work allowance.
        derivatives: Select raw derivatives; false selects Taylor coefficients.

    Returns:
        A rounded interpolant or coefficient enclosures for every data selection.

    Raises:
        PolyError: DOMAIN for invalid data or inseparable nodes, LIMIT or CAPACITY on work limits.
    """
    _validate[D](points, context, budget)
    for jet in jets:
        _validate[D](jet, context, budget)
    var coefficients = _hermite_values[_NumericOps[D]](_view(points), jets, derivatives, context, budget)
    return NumericDensePoly[D].from_coefficients(variable, coefficients^, budget)
