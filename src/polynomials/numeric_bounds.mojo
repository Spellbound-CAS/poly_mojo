"""Exact rectangular bounds and public numerical scalar conversions."""
from apn_mojo import Rational, Float, Complex, Ball, ComplexBall, FloatFormat, ArithmeticContext, ComplexContext, BallContext, add, subtract, to_float_if_certain
from .errors import Budget, PolyError, PolyErrorKind, _charge_size, _sum, _product
from .coefficients import QQ
from ._root_ops import _compare, _midpoint
from .numeric_domains import NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, _numeric_cost, _context_cost


@fieldwise_init
struct NumericBox(Copyable, Movable, Equatable, Writable):
    """A closed complex rectangle with exact rational endpoints."""
    var real_low: QQ
    """The real lower endpoint."""
    var real_high: QQ
    """The real upper endpoint."""
    var imag_low: QQ
    """The imaginary lower endpoint."""
    var imag_high: QQ
    """The imaginary upper endpoint."""

    @staticmethod
    def real_interval(low: QQ, high: QQ, mut budget: Budget) raises PolyError -> Self:
        """Construct a real interval with zero imaginary part.

        Args:
            low: The real lower endpoint.
            high: The real upper endpoint.
            budget: The shared work allowance.

        Returns:
            The validated exact rectangle.

        Raises:
            PolyError: DOMAIN for reversed bounds, LIMIT on exhaustion.
        """
        _charge_size(_sum(low.size(), high.size()), budget)
        if _compare(low, high, budget) > 0:
            raise PolyError(PolyErrorKind.DOMAIN, "reversed numerical interval")
        return Self(low.copy(), high.copy(), QQ(Rational(0)), QQ(Rational(0)))

    @staticmethod
    def point(real: QQ, imag: QQ, mut budget: Budget) raises PolyError -> Self:
        """Construct a singleton rectangle.

        Args:
            real: The exact real coordinate.
            imag: The exact imaginary coordinate.
            budget: The shared work allowance.

        Returns:
            The exact point enclosure.

        Raises:
            PolyError: LIMIT on exhaustion.
        """
        _charge_size(_sum(real.size(), imag.size()), budget)
        return Self(real.copy(), real.copy(), imag.copy(), imag.copy())

    def __eq__(self, other: Self) -> Bool:
        """Compare endpoint representations under the protocol budget exception.

        Args:
            other: The other rectangle.

        Returns:
            Whether all four endpoints agree.
        """
        return self.real_low == other.real_low and self.real_high == other.real_high and self.imag_low == other.imag_low and self.imag_high == other.imag_high

    def contains(self, other: Self, mut budget: Budget) raises PolyError -> Bool:
        """Check exact set containment, rejecting malformed bounds.

        Args:
            other: The proposed contained rectangle.
            budget: The shared work allowance.

        Returns:
            Whether both rectangles are valid and containment holds.

        Raises:
            PolyError: LIMIT on exhaustion.
        """
        return _box_contains(self, other, budget)


def _box_charge(box: NumericBox, mut budget: Budget) raises PolyError:
    _charge_size(_sum(_sum(box.real_low.size(), box.real_high.size()), _sum(box.imag_low.size(), box.imag_high.size())), budget)


def _box_valid(box: NumericBox, mut budget: Budget) raises PolyError -> Bool:
    _box_charge(box, budget)
    return _compare(box.real_low, box.real_high, budget) <= 0 and _compare(box.imag_low, box.imag_high, budget) <= 0


def _box_contains(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> Bool:
    if not _box_valid(a, budget) or not _box_valid(b, budget):
        return False
    return _compare(a.real_low, b.real_low, budget) <= 0 and _compare(b.real_high, a.real_high, budget) <= 0 and _compare(a.imag_low, b.imag_low, budget) <= 0 and _compare(b.imag_high, a.imag_high, budget) <= 0


def _box_zero(box: NumericBox, mut budget: Budget) raises PolyError -> Bool:
    _box_charge(box, budget)
    return box.real_low.value.sign() <= 0 and box.real_high.value.sign() >= 0 and box.imag_low.value.sign() <= 0 and box.imag_high.value.sign() >= 0


def _box_point(box: NumericBox, mut budget: Budget) raises PolyError -> Bool:
    _box_charge(box, budget)
    return box.real_low == box.real_high and box.imag_low == box.imag_high


def _box_real(box: NumericBox, mut budget: Budget) raises PolyError -> Bool:
    _box_charge(box, budget)
    return box.imag_low.is_zero() and box.imag_high.is_zero()


def _point(value: Int) -> NumericBox:
    return NumericBox(QQ(Rational(value)), QQ(Rational(value)), QQ(Rational(0)), QQ(Rational(0)))


def _interval_product(al: QQ, ah: QQ, bl: QQ, bh: QQ, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    var p = al.mul(bl, budget)
    var low = p.copy()
    var high = p.copy()
    for candidate in [al.mul(bh, budget), ah.mul(bl, budget), ah.mul(bh, budget)]:
        _charge_size(_sum(candidate.size(), _sum(low.size(), high.size())), budget)
        if _compare(candidate, low, budget) < 0:
            low = candidate.copy()
        if _compare(candidate, high, budget) > 0:
            high = candidate.copy()
    return (low^, high^)


def _box_add(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    return NumericBox(a.real_low.add(b.real_low, budget), a.real_high.add(b.real_high, budget), a.imag_low.add(b.imag_low, budget), a.imag_high.add(b.imag_high, budget))


def _box_sub(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    return NumericBox(a.real_low.sub(b.real_high, budget), a.real_high.sub(b.real_low, budget), a.imag_low.sub(b.imag_high, budget), a.imag_high.sub(b.imag_low, budget))


def _box_mul(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    var rr = _interval_product(a.real_low, a.real_high, b.real_low, b.real_high, budget)
    if _box_real(a, budget) and _box_real(b, budget):
        return NumericBox(rr[0].copy(), rr[1].copy(), QQ(Rational(0)), QQ(Rational(0)))
    var ii = _interval_product(a.imag_low, a.imag_high, b.imag_low, b.imag_high, budget)
    var ri = _interval_product(a.real_low, a.real_high, b.imag_low, b.imag_high, budget)
    var ir = _interval_product(a.imag_low, a.imag_high, b.real_low, b.real_high, budget)
    return NumericBox(rr[0].sub(ii[1], budget), rr[1].sub(ii[0], budget), ri[0].add(ir[0], budget), ri[1].add(ir[1], budget))


def _box_div_real(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    if not _box_real(a, budget) or not _box_real(b, budget) or _box_zero(b, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "interval division requires a nonzero real divisor")
    var inv = NumericBox.real_interval(b.real_high.inverse(budget), b.real_low.inverse(budget), budget)
    return _box_mul(a, inv, budget)


def _box_center(a: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    return NumericBox.point(_midpoint(a.real_low, a.real_high, budget), _midpoint(a.imag_low, a.imag_high, budget), budget)


def _ball_endpoints(value: Ball, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    # A finite Float is an integer times a power of two. Align both operands
    # down to the smaller least-significant exponent and allow one carry bit.
    # At that precision both sums are exact; conversion then removes powers
    # of two directly, avoiding general rational denominator normalization.
    var midpoint = value.midpoint()
    var radius: Float
    try:
        if value.is_exact():
            var q = QQ(midpoint.to_rational_exact())
            return (q.copy(), q^)
        radius = value.radius()
        if midpoint.is_zero():
            var q = radius.to_rational_exact()
            return (QQ(-q), QQ(q^))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var me: Int
    var re: Int
    try:
        me = midpoint.exponent()
        re = radius.exponent()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    # _evaluation_bounds precharges _numeric_cost, establishing that each
    # precision + abs(exponent) fits Int before these signed subtractions.
    var low = min(me - midpoint.precision(), re - radius.precision())
    var high = max(me, re)
    var distance = _sum(high, -low) if high > 0 and low < 0 else high - low
    var bits = _sum(distance, 1)
    _charge_size(_sum(bits, bits), budget)
    try:
        var context = ArithmeticContext(format=FloatFormat(bits, emin=low + 1, emax=high + 1), trap_inexact=True, trap_underflow=True, trap_overflow=True, trap_invalid=True)
        return (QQ(subtract(midpoint, radius, context=context).to_rational_exact()), QQ(add(midpoint, radius, context=context).to_rational_exact()))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _evaluation_bounds[D: NumericDomain](value: D.Scalar, mut budget: Budget) raises PolyError -> NumericBox:
    # Producer-only acceleration. Public bounds and independent checkers use
    # rational midpoint +/- radius, without invoking binary arithmetic.
    _charge_size(_numeric_cost[D](value), budget)
    comptime if D == BallDomain:
        var ends = _ball_endpoints(rebind[Ball](value), budget)
        return NumericBox(ends[0].copy(), ends[1].copy(), QQ(Rational(0)), QQ(Rational(0)))
    elif D == ComplexBallDomain:
        var z = rebind[ComplexBall](value)
        var real = _ball_endpoints(z.real(), budget)
        var imag = _ball_endpoints(z.imag(), budget)
        return NumericBox(real[0].copy(), real[1].copy(), imag[0].copy(), imag[1].copy())
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "evaluation bounds require a ball domain")


def _rational_ball_endpoints(value: Ball) raises -> Tuple[QQ, QQ]:
    var midpoint = value.midpoint_rational()
    if value.is_exact():
        return (QQ(midpoint.copy()), QQ(midpoint^))
    var radius = value.radius().to_rational_exact()
    return (QQ(midpoint - radius), QQ(midpoint + radius))


def numeric_bounds[D: NumericDomain](value: D.Scalar, mut budget: Budget) raises PolyError -> NumericBox:
    """Extract exact rational endpoints from a finite apn_mojo numerical value.

    Parameters:
        D: The numerical domain adapter.

    Args:
        value: The point value or bounded enclosure.
        budget: The shared work allowance.

    Returns:
        The exact stored point or represented rectangle.

    Raises:
        PolyError: DOMAIN for nonfinite values, LIMIT or CAPACITY on work limits.
    """
    _charge_size(_numeric_cost[D](value), budget)
    try:
        comptime if D == FloatDomain:
            var x = QQ(rebind[Float](value).to_rational_exact())
            return NumericBox(x.copy(), x.copy(), QQ(Rational(0)), QQ(Rational(0)))
        elif D == ComplexDomain:
            var z = rebind[Complex](value)
            var x = QQ(z.real().to_rational_exact())
            var y = QQ(z.imag().to_rational_exact())
            return NumericBox(x.copy(), x.copy(), y.copy(), y.copy())
        elif D == BallDomain:
            var ends = _rational_ball_endpoints(rebind[Ball](value))
            return NumericBox(ends[0].copy(), ends[1].copy(), QQ(Rational(0)), QQ(Rational(0)))
        elif D == ComplexBallDomain:
            var z = rebind[ComplexBall](value)
            var real = _rational_ball_endpoints(z.real())
            var imag = _rational_ball_endpoints(z.imag())
            return NumericBox(real[0].copy(), real[1].copy(), imag[0].copy(), imag[1].copy())
        else:
            raise Error("unknown numerical domain")
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _real_ball(low: QQ, high: QQ, precision: Int, mut budget: Budget) raises PolyError -> Ball:
    var midpoint = _midpoint(low, high, budget)
    var radius = _midpoint(QQ(Rational(0)), high.sub(low, budget), budget)
    _charge_size(_sum(precision, _sum(midpoint.size(), radius.size())), budget)
    try:
        return Ball(midpoint.value, radius.value, precision=precision)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def numeric_value[D: NumericDomain](bounds: NumericBox, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    """Convert exact bounds to an enclosure or a certainly rounded point value.

    Parameters:
        D: The destination numerical domain adapter.

    Args:
        bounds: The exact point or interval to convert.
        context: Explicit output precision and rounding settings.
        budget: The shared work allowance.

    Returns:
        An enclosing ball or one rounded point shared by every enclosed value.

    Raises:
        PolyError: DOMAIN for uncertain rounding or nonreal real-domain input, LIMIT or CAPACITY on work limits.
    """
    _charge_size(_context_cost[D](context), budget)
    if not _box_valid(bounds, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical rectangle")
    comptime if D == BallDomain or D == FloatDomain:
        if not _box_real(bounds, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "nonreal numerical value in a real domain")
    comptime if D == BallDomain:
        return rebind_var[D.Scalar](_real_ball(bounds.real_low, bounds.real_high, rebind[BallContext](context).precision().value(), budget))
    elif D == ComplexBallDomain:
        var bits = rebind[BallContext](context).precision().value()
        var real = _real_ball(bounds.real_low, bounds.real_high, bits, budget)
        var imag = _real_ball(bounds.imag_low, bounds.imag_high, bits, budget)
        try:
            return rebind_var[D.Scalar](ComplexBall(real, imag))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    elif D == FloatDomain:
        var c = rebind[ArithmeticContext](context)
        var low: Float
        if bounds.real_low == bounds.real_high:
            try:
                low = Float(bounds.real_low.value, context=c)
            except e:
                raise PolyError(PolyErrorKind.DOMAIN, String(e))
        else:
            # Preserve dyadic endpoints exactly when forming the working ball;
            # extra outward error for other rationals can only reject certainty.
            var bits = max(_sum(c.format().precision(), 32), _product(_sum(bounds.real_low.size(), bounds.real_high.size()), 8))
            var enclosure = _real_ball(bounds.real_low, bounds.real_high, bits, budget)
            _charge_size(_sum(bits, c.format().precision()), budget)
            var certain: Optional[Float]
            try:
                certain = to_float_if_certain(enclosure, context=c)
            except e:
                raise PolyError(PolyErrorKind.DOMAIN, String(e))
            if not certain:
                raise PolyError(PolyErrorKind.DOMAIN, "interval does not have a certainly rounded Float")
            low = certain.take()
        if not low.is_finite():
            raise PolyError(PolyErrorKind.DOMAIN, "interval does not round to a unique finite Float")
        return rebind_var[D.Scalar](low^)
    elif D == ComplexDomain:
        var c = rebind[ComplexContext](context)
        var real_bounds = NumericBox.real_interval(bounds.real_low, bounds.real_high, budget)
        var imag_bounds = NumericBox.real_interval(bounds.imag_low, bounds.imag_high, budget)
        var real = numeric_value[FloatDomain](real_bounds, c.real(), budget)
        var imag = numeric_value[FloatDomain](imag_bounds, c.imag(), budget)
        try:
            return rebind_var[D.Scalar](Complex(real, imag, context=c))
        except e:
            raise PolyError(PolyErrorKind.DOMAIN, String(e))
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "unknown numerical domain")
