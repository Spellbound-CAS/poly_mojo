"""Exact residual bounds, coefficient sensitivity and certified Newton candidates."""
from apn_mojo import Rational, Complex, ComplexContext
from .coefficients import CoefficientRing, QQ
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from .numeric_domains import NumericDomain, ComplexDomain, _numeric_cost, _numeric_apply
from .numeric_dense import NumericDensePoly
from .numeric_bounds import NumericBox, numeric_bounds, numeric_value, _box_point, _box_valid
from .numeric_evaluate import _exact_dense_boxes, _numeric_dense_boxes
from .numeric_options import NumericOptions, NumericOutcome, _norm_lower, _norm_upper, _accurate
from ._numeric_roots import _exact_box_evaluate, _derivative_boxes
from .numeric_complex_roots import ComplexRootRegion, _representative, _contour, _check_region
from ._coefficient_ops import _copy
from ._root_ops import _compare


@fieldwise_init
struct RootDiagnostics(Copyable, Movable):
    """Rational residual and local sensitivity bounds at a specified point."""
    var residual: NumericBox
    """An enclosure of p(z), uniformly over any coefficient family."""
    var derivative: NumericBox
    """An enclosure of p'(z)."""
    var backward_error_low: QQ
    """Lower bound on the absolute complex coefficient infinity-norm backward error."""
    var backward_error_high: QQ
    """Upper bound on that backward error, |p(z)| divided by sum |z|^k."""
    var newton_correction_upper: Optional[QQ]
    """Upper bound for |p(z)/p'(z)| when the derivative enclosure excludes zero."""
    var coefficient_sensitivity_upper: Optional[QQ]
    """Upper bound for sum |z|^k / |p'(z)|; a root condition bound only when z is a root."""


def _limits(residual: NumericBox, derivative: NumericBox, point: NumericBox, count: Int, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ, Optional[QQ], Optional[QQ]]:
    var low = QQ(Rational(0))
    var high = QQ(Rational(0))
    var plow = QQ(Rational(1))
    var phigh = QQ(Rational(1))
    var zlow = _norm_lower(point, budget)
    var zhigh = _norm_upper(point, budget)
    for _ in range(max(1, count)):
        low = low.add(plow, budget)
        high = high.add(phigh, budget)
        plow = plow.mul(zlow, budget)
        phigh = phigh.mul(zhigh, budget)
    var residual_low = _norm_lower(residual, budget)
    var residual_high = _norm_upper(residual, budget)
    var derivative_low = _norm_lower(derivative, budget)
    var correction = Optional[QQ]()
    var sensitivity = Optional[QQ]()
    if not derivative_low.is_zero():
        var inverse = derivative_low.inverse(budget)
        correction = residual_high.mul(inverse, budget)
        sensitivity = high.mul(inverse, budget)
    return (residual_low.mul(high.inverse(budget), budget), residual_high.mul(low.inverse(budget), budget), correction^, sensitivity^)


def _diagnostics(coefficients: List[NumericBox], point: NumericBox, mut budget: Budget) raises PolyError -> RootDiagnostics:
    if not _box_valid(point, budget) or not _box_point(point, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "root diagnostics require a finite point rectangle")
    var residual = _exact_box_evaluate(coefficients, point, budget)
    var derivative = _exact_box_evaluate(_derivative_boxes(coefficients, budget), point, budget)
    var limits = _limits(residual, derivative, point, len(coefficients), budget)
    # Tuple element access borrows in Mojo 1.1; charge the returned bound copies.
    var lower = _copy(limits[0], budget)
    var upper = _copy(limits[1], budget)
    var correction = Optional[QQ]()
    var sensitivity = Optional[QQ]()
    if limits[2]:
        correction = _copy(limits[2].value(), budget)
        sensitivity = _copy(limits[3].value(), budget)
    return RootDiagnostics(residual^, derivative^, lower^, upper^, correction^, sensitivity^)


def _check_diagnostics(coefficients: List[NumericBox], point: NumericBox, result: RootDiagnostics, mut budget: Budget) raises PolyError -> Bool:
    if not _box_valid(point, budget) or not _box_point(point, budget):
        return False
    if not result.residual.contains(_exact_box_evaluate(coefficients, point, budget), budget) or not result.derivative.contains(_exact_box_evaluate(_derivative_boxes(coefficients, budget), point, budget), budget):
        return False
    var limits = _limits(result.residual, result.derivative, point, len(coefficients), budget)
    if result.backward_error_low.value.sign() < 0 or _compare(result.backward_error_low, limits[0], budget) > 0 or _compare(result.backward_error_high, limits[1], budget) < 0:
        return False
    if Bool(result.newton_correction_upper) != Bool(limits[2]) or Bool(result.coefficient_sensitivity_upper) != Bool(limits[3]):
        return False
    if result.newton_correction_upper:
        if _compare(result.newton_correction_upper.value(), limits[2].value(), budget) < 0 or _compare(result.coefficient_sensitivity_upper.value(), limits[3].value(), budget) < 0:
            return False
    return True


@fieldwise_init
struct PolishedRootResult(Copyable, Movable):
    """A numerical Newton candidate and an optional independently checkable root witness."""
    var candidate: Complex
    """The rounded advisory candidate; its closeness alone is not a root proof."""
    var iterations: Int
    """The number of completed Newton updates."""
    var region: Optional[ComplexRootRegion]
    """A rectangle with a contour proof of exactly one root, when available."""
    var outcome: NumericOutcome
    """SUCCESS only when one root is validated and the requested width is met."""


def _polish(coefficients: List[NumericBox], seed: Complex, radius: QQ, iterations: Int, context: ComplexContext, options: NumericOptions, mut budget: Budget) raises PolyError -> PolishedRootResult:
    options._validate()
    _charge_size(radius.size(), budget)
    if radius.value.sign() <= 0 or iterations < 0:
        raise PolyError(PolyErrorKind.DOMAIN, "polishing requires a positive radius and nonnegative iteration limit")
    var representative = _representative(coefficients, budget)
    _allocation[Complex](len(coefficients), budget)
    var values = List[Complex](capacity=len(coefficients))
    for box in representative:
        values.append(numeric_value[ComplexDomain](box, context, budget))
    var poly = NumericDensePoly[ComplexDomain].from_coefficients(0, values^, budget)
    _charge_size(_numeric_cost[ComplexDomain](seed), budget)
    var candidate = seed.copy()
    var completed = 0
    for _ in range(iterations):
        var evaluated = poly.evaluate_with_derivative(candidate, context, budget)
        ref value = evaluated[0]
        ref slope = evaluated[1]
        _charge_size(_sum(_numeric_cost[ComplexDomain](value), _numeric_cost[ComplexDomain](slope)), budget)
        if value.is_zero() or slope.is_zero():
            break
        var correction = _numeric_apply[ComplexDomain](value, slope, 3, context, budget)
        candidate = _numeric_apply[ComplexDomain](candidate, correction, 1, context, budget)
        completed += 1
    var center = numeric_bounds[ComplexDomain](candidate, budget)
    var rectangle = NumericBox(center.real_low.sub(radius, budget), center.real_high.add(radius, budget), center.imag_low.sub(radius, budget), center.imag_high.add(radius, budget))
    var subdivisions = 0
    var region = _contour(coefficients, representative, rectangle, options, options.max_subdivisions, subdivisions, budget)
    if not region or region.value().count != 1:
        return PolishedRootResult(candidate^, completed, None, NumericOutcome.UNRESOLVED_BOUNDARY)
    var outcome = NumericOutcome.SUCCESS if _accurate(rectangle, options, budget) else NumericOutcome.SUBDIVISION_LIMIT
    return PolishedRootResult(candidate^, completed, region^, outcome)


def _check_polished(coefficients: List[NumericBox], result: PolishedRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    options._validate()
    if result.iterations < 0 or not result.region:
        return False
    ref region = result.region.value()
    if region.count != 1 or not region.rectangle.contains(numeric_bounds[ComplexDomain](result.candidate, budget), budget):
        return False
    if not _check_region(coefficients, _representative(coefficients, budget), region, budget):
        return False
    var accurate = _accurate(region.rectangle, options, budget)
    return result.outcome == (NumericOutcome.SUCCESS if accurate else NumericOutcome.SUBDIVISION_LIMIT)


def root_diagnostics_exact[C: CoefficientRing](poly: DensePoly[C], point: NumericBox, mut budget: Budget) raises PolyError -> RootDiagnostics:
    """Bound residuals, coefficient backward error and local derivative sensitivity.

    Parameters:
        C: ZZ, QQ, or QQi.

    Args:
        poly: The input polynomial.
        point: The exact point rectangle where diagnostics are requested.
        budget: The shared work allowance.

    Returns:
        Uniform rational bounds; local sensitivity becomes root conditioning only at a root.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _diagnostics(_exact_dense_boxes(poly, budget), point, budget)


def check_root_diagnostics_exact[C: CoefficientRing](poly: DensePoly[C], point: NumericBox, result: RootDiagnostics, mut budget: Budget) raises PolyError -> Bool:
    """Check exact residual and local sensitivity bounds.

    Parameters:
        C: ZZ, QQ, or QQi.

    Args:
        poly: The input polynomial.
        point: The exact point rectangle where diagnostics are requested.
        result: The claimed residual and local bounds.
        budget: The shared work allowance.

    Returns:
        Whether the supplied bounds or root witness are valid.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _check_diagnostics(_exact_dense_boxes(poly, budget), point, result, budget)


def polish_root_exact[C: CoefficientRing](poly: DensePoly[C], seed: Complex, radius: QQ, iterations: Int, context: ComplexContext, options: NumericOptions, mut budget: Budget) raises PolyError -> PolishedRootResult:
    """Polish a Newton candidate and attempt an independent one-root contour proof.

    Parameters:
        C: ZZ, QQ, or QQi.

    Args:
        poly: The input polynomial.
        seed: The advisory complex starting point.
        radius: A positive exact radius for the final validation rectangle.
        iterations: The nonnegative Newton update limit.
        context: Explicit candidate arithmetic settings.
        options: Contour precision, subdivision and accuracy limits.
        budget: The shared work allowance.

    Returns:
        An advisory candidate and a root certificate only if validation succeeds.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _polish(_exact_dense_boxes(poly, budget), seed, radius, iterations, context, options, budget)


def check_polished_root_exact[C: CoefficientRing](poly: DensePoly[C], result: PolishedRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    """Check the one-root contour and candidate containment independently.

    Parameters:
        C: ZZ, QQ, or QQi.

    Args:
        poly: The input polynomial.
        result: The proposed polished root and local witness.
        options: The requested contour precision and accuracy.
        budget: The shared work allowance.

    Returns:
        Whether the supplied bounds or root witness are valid.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _check_polished(_exact_dense_boxes(poly, budget), result, options, budget)


def root_diagnostics_numeric[D: NumericDomain](poly: NumericDensePoly[D], point: NumericBox, mut budget: Budget) raises PolyError -> RootDiagnostics:
    """Bound residuals, coefficient backward error and local derivative sensitivity.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        poly: The input polynomial or coefficient family.
        point: The exact point rectangle where diagnostics are requested.
        budget: The shared work allowance.

    Returns:
        Uniform rational bounds; local sensitivity becomes root conditioning only at a root.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _diagnostics(_numeric_dense_boxes(poly, budget), point, budget)


def check_root_diagnostics_numeric[D: NumericDomain](poly: NumericDensePoly[D], point: NumericBox, result: RootDiagnostics, mut budget: Budget) raises PolyError -> Bool:
    """Check exact residual and local sensitivity bounds.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        poly: The input polynomial or coefficient family.
        point: The exact point rectangle where diagnostics are requested.
        result: The claimed residual and local bounds.
        budget: The shared work allowance.

    Returns:
        Whether the supplied bounds or root witness are valid.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _check_diagnostics(_numeric_dense_boxes(poly, budget), point, result, budget)


def polish_root_numeric[D: NumericDomain](poly: NumericDensePoly[D], seed: Complex, radius: QQ, iterations: Int, context: ComplexContext, options: NumericOptions, mut budget: Budget) raises PolyError -> PolishedRootResult:
    """Polish a Newton candidate and attempt an independent one-root contour proof.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        poly: The input polynomial or coefficient family.
        seed: The advisory complex starting point.
        radius: A positive exact radius for the final validation rectangle.
        iterations: The nonnegative Newton update limit.
        context: Explicit candidate arithmetic settings.
        options: Contour precision, subdivision and accuracy limits.
        budget: The shared work allowance.

    Returns:
        An advisory candidate and a root certificate only if validation succeeds.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _polish(_numeric_dense_boxes(poly, budget), seed, radius, iterations, context, options, budget)


def check_polished_root_numeric[D: NumericDomain](poly: NumericDensePoly[D], result: PolishedRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    """Check the one-root contour and candidate containment independently.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        poly: The input polynomial or coefficient family.
        result: The proposed polished root and local witness.
        options: The requested contour precision and accuracy.
        budget: The shared work allowance.

    Returns:
        Whether the supplied bounds or root witness are valid.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _check_polished(_numeric_dense_boxes(poly, budget), result, options, budget)
