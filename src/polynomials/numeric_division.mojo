"""Rounded long division and uniformly enclosing division with local certificates."""
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .numeric_domains import NumericDomain, BallDomain, ComplexBallDomain, _numeric_zero, _context_cost
from .numeric_dense import NumericDensePoly
from .numeric_bounds import NumericBox, numeric_bounds, numeric_value, _box_zero, _box_charge, _box_sub, _box_mul, _point
from .numeric_evaluate import _numeric_dense_boxes
from ._toolkit_numeric import _NumericOps
from ._toolkit_ops import _values
from ._numeric_division_ops import _copy_boxes, _box_divide, _residual_boxes


@fieldwise_init
struct NumericDivisionStep(Copyable, Movable):
    """One enclosing leading-term elimination in descending quotient order."""
    var quotient: NumericBox
    """Enclosure of the exact quotient coefficient at this step."""
    var tail: List[NumericBox]
    """Updated lower coefficients, ascending within the shifted divisor support."""


@fieldwise_init
struct NumericDivisionResult[D: NumericDomain](Copyable, Movable):
    """Numerical quotient and remainder with explicit residual or family evidence.

    Parameters:
        D: The numerical coefficient domain.
    """
    var quotient: NumericDensePoly[Self.D]
    """The rounded quotient or coefficient enclosures of all exact quotients."""
    var remainder: NumericDensePoly[Self.D]
    """The rounded remainder or enclosure, with degree below the divisor."""
    var residual: List[NumericBox]
    """Exact rectangular coefficient bounds for f-q*g-r, including high-degree rounding residuals."""
    var steps: List[NumericDivisionStep]
    """Local interval eliminations for ball domains; empty for rounded point domains."""
    var uniform: Bool
    """True only when q and r enclose exact Euclidean division for every input pair."""


def _boxes[D: NumericDomain](poly: NumericDensePoly[D], mut budget: Budget) raises PolyError -> List[NumericBox]:
    var result = _numeric_dense_boxes(poly, budget)
    if len(poly.coefficients) and _numeric_zero[D](poly.coefficients[len(poly.coefficients) - 1]):
        raise PolyError(PolyErrorKind.DOMAIN, "noncanonical numerical division input")
    return result^


def numeric_divmod[D: NumericDomain](f: NumericDensePoly[D], g: NumericDensePoly[D], context: D.Context, mut budget: Budget) raises PolyError -> NumericDivisionResult[D]:
    """Divide numerical polynomials with a uniformly nonzero leading divisor.

    Parameters:
        D: Float, Complex, Ball or ComplexBall domain.

    Args:
        f: The dividend polynomial or coefficient family.
        g: The divisor in the same variable, with leading enclosure excluding zero.
        context: Explicit rounding or enclosing arithmetic settings.
        budget: The shared work allowance.

    Returns:
        Rounded point q/r plus a rigorous residual, or enclosures of exact family q/r with elimination witnesses.

    Raises:
        PolyError: DOMAIN for zero or uncertain-degree divisors, LIMIT or CAPACITY on work limits.
    """
    f._same(g)
    _charge_size(_context_cost[D](context), budget)
    var a = _boxes(f, budget)
    var b = _boxes(g, budget)
    if not len(b) or _box_zero(b[len(b) - 1], budget):
        raise PolyError(PolyErrorKind.DOMAIN, "numerical divisor must have a uniformly nonzero leading coefficient")
    var degree = len(b) - 1
    var count = max(0, len(a) - degree)
    var q = _values[_NumericOps[D]](_view(List[D.Scalar]()), count, context, budget)
    var working = List[D.Scalar]()
    var current = List[NumericBox]()
    var uniform = False
    comptime if D == BallDomain or D == ComplexBallDomain:
        current = _copy_boxes(a, budget)
        uniform = True
    else:
        working = _values[_NumericOps[D]](_view(f.coefficients), len(f.coefficients), context, budget)
    _allocation[NumericDivisionStep](count if uniform else 0, budget)
    var steps = List[NumericDivisionStep]()
    for shift in range(count - 1, -1, -1):
        comptime if D == BallDomain or D == ComplexBallDomain:
            var exact = _box_divide(current[shift + degree], b[degree], budget)
            q[shift] = numeric_value[D](exact, context, budget)
            var quotient = numeric_bounds[D](q[shift], budget)
            _allocation[NumericBox](degree, budget)
            var tail = List[NumericBox](capacity=degree)
            for j in range(degree):
                var next = _box_sub(current[shift + j], _box_mul(quotient, b[j], budget), budget)
                var rounded = numeric_value[D](next, context, budget)
                var enclosure = numeric_bounds[D](rounded, budget)
                _box_charge(enclosure, budget)
                current[shift + j] = enclosure.copy()
                tail.append(enclosure^)
            current[shift + degree] = _point(0)
            steps.append(NumericDivisionStep(quotient^, tail^))
        else:
            q[shift] = _NumericOps[D].apply(working[shift + degree], g.coefficients[degree], 3, context, budget)
            for j in range(degree):
                var product = _NumericOps[D].apply(q[shift], g.coefficients[j], 2, context, budget)
                working[shift + j] = _NumericOps[D].apply(working[shift + j], product, 1, context, budget)
            # Its rounding error is retained separately in the full residual.
            working[shift + degree] = _NumericOps[D].integer(0, context, budget)
    var length = min(len(a), degree)
    _allocation[D.Scalar](length, budget)
    var r = List[D.Scalar](capacity=length)
    for j in range(length):
        comptime if D == BallDomain or D == ComplexBallDomain:
            r.append(numeric_value[D](current[j], context, budget))
        else:
            r.append(_NumericOps[D].copy(working[j], budget))
    var quotient = NumericDensePoly[D].from_coefficients(f.variable, q^, budget)
    var remainder = NumericDensePoly[D].from_coefficients(f.variable, r^, budget)
    var residual = _residual_boxes(a, b, _boxes(quotient, budget), _boxes(remainder, budget), budget)
    return NumericDivisionResult[D](quotient^, remainder^, residual^, steps^, uniform)


def check_numeric_division[D: NumericDomain](f: NumericDensePoly[D], g: NumericDensePoly[D], result: NumericDivisionResult[D], mut budget: Budget) raises PolyError -> Bool:
    """Check division bounds by exact residual identities and local interval steps.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        f: The original dividend.
        g: The original divisor.
        result: The proposed quotient, remainder and evidence.
        budget: The shared work allowance.

    Returns:
        Whether degree bounds, residual enclosures and all claimed family guarantees hold.

    Raises:
        PolyError: DOMAIN for malformed values or zero divisors, LIMIT or CAPACITY on work limits.
    """
    f._same(g)
    if result.quotient.variable != f.variable or result.remainder.variable != f.variable:
        return False
    var a = _boxes(f, budget)
    var b = _boxes(g, budget)
    var q = _boxes(result.quotient, budget)
    var r = _boxes(result.remainder, budget)
    if not len(b) or _box_zero(b[len(b) - 1], budget):
        return False
    var degree = len(b) - 1
    var count = max(0, len(a) - degree)
    if len(q) > count or len(r) > degree:
        return False
    var enclosing = False
    comptime if D == BallDomain or D == ComplexBallDomain:
        enclosing = True
    if result.uniform != enclosing or len(result.steps) != (count if enclosing else 0):
        return False
    if enclosing:
        var current = _copy_boxes(a, budget)
        for k in range(count):
            var shift = count - 1 - k
            ref step = result.steps[k]
            if len(step.tail) != degree or not step.quotient.contains(_box_divide(current[shift + degree], b[degree], budget), budget):
                return False
            var stored = _point(0)
            if shift < len(q):
                _box_charge(q[shift], budget)
                stored = q[shift].copy()
            if not stored.contains(step.quotient, budget):
                return False
            for j in range(degree):
                var next = _box_sub(current[shift + j], _box_mul(step.quotient, b[j], budget), budget)
                if not step.tail[j].contains(next, budget):
                    return False
                _box_charge(step.tail[j], budget)
                current[shift + j] = step.tail[j].copy()
            current[shift + degree] = _point(0)
        for j in range(min(len(a), degree)):
            var stored = _point(0)
            if j < len(r):
                _box_charge(r[j], budget)
                stored = r[j].copy()
            if not stored.contains(current[j], budget):
                return False
    var residual = _residual_boxes(a, b, q, r, budget)
    if len(result.residual) != len(residual):
        return False
    for i in range(len(residual)):
        if not result.residual[i].contains(residual[i], budget):
            return False
    return True
