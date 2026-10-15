"""Rounded residuals and uniform family quotient/remainder certificates."""
from apn_mojo import Rational, ArithmeticContext, FloatFormat, ComplexContext, BallContext, Ball, ComplexBall, Float, Complex
from polynomials import Budget, PolyErrorKind, QQ, NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, NumericDensePoly, NumericBox, numeric_value, numeric_bounds, numeric_divmod, check_numeric_division
from support.check import Checker


def point[D: NumericDomain](value: Int, context: D.Context, mut work: Budget) raises -> D.Scalar:
    return numeric_value[D](NumericBox.point(QQ(Rational(value)), QQ(Rational(0)), work), context, work)


def domains[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    var f = NumericDensePoly[D].from_coefficients(0, [point[D](1, context, work), point[D](2, context, work), point[D](1, context, work)], work)
    var g = NumericDensePoly[D].from_coefficients(0, [point[D](1, context, work), point[D](1, context, work)], work)
    var result = numeric_divmod(f, g, context, work)
    c.check(check_numeric_division(f, g, result, work), "checked numerical division")
    c.check(result.quotient.degree_bound() == 1 and result.remainder.degree_bound() == -1, "exact dyadic division shape")
    for value in result.quotient.coefficients:
        c.check(numeric_bounds[D](value, work).contains(NumericBox.point(QQ(Rational(1)), QQ(Rational(0)), work), work), "hand-derived quotient coefficient")
    var wrong = result.copy()
    wrong.quotient.coefficients[0] = point[D](99, context, work)
    c.check(not check_numeric_division(f, g, wrong, work), "corrupt quotient rejected")
    var tiny = Budget(0)
    try:
        _ = numeric_divmod(f, g, context, tiny)
        c.check(False, "division ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "division budget preserves input")


def main() raises:
    var c = Checker("numeric_division")
    var work = Budget()
    domains[FloatDomain](ArithmeticContext(), c, work)
    domains[ComplexDomain](ComplexContext(ArithmeticContext()), c, work)
    domains[BallDomain](BallContext(128), c, work)
    domains[ComplexBallDomain](BallContext(128), c, work)
    var context = BallContext(128)
    var interval = Ball(Rational(3, 2), Rational(1, 2))
    var f = NumericDensePoly[BallDomain].from_coefficients(0, [interval.copy(), Ball(0), Ball(1)], work)
    var g = NumericDensePoly[BallDomain].from_coefficients(0, [interval.copy(), Ball(1)], work)
    var family = numeric_divmod(f, g, context, work)
    c.check(family.uniform and check_numeric_division(f, g, family, work), "uniform family division certificate")
    c.check(numeric_bounds[BallDomain](family.quotient.coefficients[0], work).contains(NumericBox.real_interval(QQ(Rational(-2)), QQ(Rational(-1)), work), work), "family quotient enclosure")
    c.check(numeric_bounds[BallDomain](family.remainder.coefficients[0], work).contains(NumericBox.real_interval(QQ(Rational(2)), QQ(Rational(6)), work), work), "family remainder enclosure")
    var bad = family.copy()
    bad.steps[0].quotient = NumericBox.point(QQ(Rational(99)), QQ(Rational(0)), work)
    c.check(not check_numeric_division(f, g, bad, work), "forged uniform division step rejected")
    var uncertain = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(1), Ball(0, Rational(1))], work)
    try:
        _ = numeric_divmod(f, uncertain, context, work)
        c.check(False, "uncertain divisor degree accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "leading enclosure containing zero rejected")
    var low = ArithmeticContext(format=FloatFormat(4))
    var a = NumericDensePoly[FloatDomain].from_coefficients(0, [Float(1), Float(1)], work)
    var b = NumericDensePoly[FloatDomain].from_coefficients(0, [Float(3)], work)
    var rounded = numeric_divmod(a, b, low, work)
    c.check(not rounded.uniform and check_numeric_division(a, b, rounded, work), "rounded division checked through residual")
    c.check(rounded.residual[0] != NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), work), "rounding residual is retained")
    var ca = NumericDensePoly[ComplexBallDomain].from_coefficients(0, [ComplexBall(1)], work)
    var cb = NumericDensePoly[ComplexBallDomain].from_coefficients(0, [ComplexBall(Complex(1, 1))], work)
    var complex_result = numeric_divmod(ca, cb, context, work)
    c.check(check_numeric_division(ca, cb, complex_result, work) and numeric_bounds[ComplexBallDomain](complex_result.quotient.coefficients[0], work).contains(NumericBox.point(QQ(Rational(1, 2)), QQ(Rational(-1, 2)), work), work), "complex interval division certificate")
    c.finish()
