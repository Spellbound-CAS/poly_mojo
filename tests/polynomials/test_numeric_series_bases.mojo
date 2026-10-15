"""Numerical basis and series operations across all four scalar domains."""
from apn_mojo import Rational, ArithmeticContext, ComplexContext, BallContext
from polynomials import Budget, PolyErrorKind, QQ, NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, NumericBox, NumericDensePoly, NumericSeries, NumericBasisPolynomial, numeric_value, numeric_bounds
from support.check import Checker


def point[D: NumericDomain](value: Int, context: D.Context, mut work: Budget) raises -> D.Scalar:
    return numeric_value[D](NumericBox.point(QQ(Rational(value)), QQ(Rational(0)), work), context, work)


def contains[D: NumericDomain](value: D.Scalar, real: QQ, mut work: Budget) raises -> Bool:
    return numeric_bounds[D](value, work).contains(NumericBox.point(real, QQ(Rational(0)), work), work)


def domains[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    var f = NumericDensePoly[D].from_coefficients(0, [point[D](1, context, work), point[D](2, context, work), point[D](3, context, work)], work)
    var nodes: List[D.Scalar] = [point[D](1, context, work), point[D](2, context, work)]
    var newton = NumericBasisPolynomial[D].newton(f, nodes, context, work)
    var chebyshev = NumericBasisPolynomial[D].chebyshev(f, context, work)
    var bernstein = NumericBasisPolynomial[D].bernstein(f, 2, context, work)
    c.check(contains[D](newton.coefficients[1], QQ(Rational(11)), work), "numeric Newton coefficient")
    c.check(contains[D](chebyshev.coefficients[0], QQ(Rational(5, 2)), work), "numeric Chebyshev coefficient")
    c.check(contains[D](bernstein.coefficients[1], QQ(Rational(2)), work), "numeric Bernstein coefficient")
    c.check(contains[D](newton.evaluate(point[D](2, context, work), context, work), QQ(Rational(17)), work), "numeric Newton evaluation")
    c.check(contains[D](chebyshev.evaluate(point[D](2, context, work), context, work), QQ(Rational(17)), work), "numeric Clenshaw evaluation")
    c.check(contains[D](bernstein.evaluate(point[D](2, context, work), context, work), QQ(Rational(17)), work), "numeric de Casteljau evaluation")
    var mono = chebyshev.to_polynomial(context, work)
    for i in range(3):
        c.check(contains[D](mono.coefficients[i], QQ(Rational(i + 1)), work), "numeric basis expansion")
    var half = numeric_value[D](NumericBox.point(QQ(Rational(1, 2)), QQ(Rational(0)), work), context, work)
    var split = bernstein.subdivide(half, context, work)
    c.check(contains[D](split[0].evaluate(point[D](1, context, work), context, work), QQ(Rational(11, 4)), work), "numerical subdivision endpoint")
    c.check(contains[D](split[1].evaluate(point[D](0, context, work), context, work), QQ(Rational(11, 4)), work), "numerical shared subdivision endpoint")
    var linear = NumericDensePoly[D].from_coefficients(0, [point[D](1, context, work), point[D](-1, context, work)], work)
    var geometric = NumericSeries[D].from_polynomial(linear, 5, context, work).inverse(context, work)
    for coefficient in geometric.polynomial.coefficients:
        c.check(contains[D](coefficient, QQ(Rational(1)), work), "numerical geometric series")
    var x = NumericDensePoly[D].from_coefficients(0, [point[D](0, context, work), point[D](1, context, work)], work)
    var series = NumericSeries[D].from_polynomial(x, 3, context, work)
    var exponential = series.exp(context, work)
    c.check(contains[D](exponential.polynomial.coefficients[2], QQ(Rational(1, 2)), work), "numerical exponential coefficients")
    c.check(contains[D](exponential.log(context, work).polynomial.coefficients[1], QQ(Rational(1)), work), "numerical formal logarithm")
    var quadratic = NumericDensePoly[D].from_coefficients(0, [point[D](0, context, work), point[D](1, context, work), point[D](1, context, work)], work)
    var reverted = NumericSeries[D].from_polynomial(quadratic, 5, context, work).revert(context, work)
    c.check(contains[D](reverted.polynomial.coefficients[4], QQ(Rational(-5)), work), "numerical reversion coefficient")
    var product = geometric.mul(NumericSeries[D].from_polynomial(linear, 5, context, work), context, work)
    c.check(contains[D](product.polynomial.coefficients[0], QQ(Rational(1)), work), "numerical truncated multiplication")
    var quotient = product.divide(NumericSeries[D].from_polynomial(linear, 5, context, work), context, work)
    c.check(contains[D](quotient.polynomial.coefficients[4], QQ(Rational(1)), work), "numerical truncated division")
    var difference = quotient.add(geometric, context, work).sub(geometric, context, work)
    c.check(contains[D](difference.polynomial.coefficients[2], QQ(Rational(1)), work), "numerical series sum and difference")
    var composed = NumericSeries[D].from_polynomial(quadratic, 5, context, work).compose(reverted, context, work)
    c.check(contains[D](composed.polynomial.coefficients[1], QQ(Rational(1)), work), "numerical series composition")
    c.check(series.integrate(context, work).derivative(context, work).order == 3 and series.truncate(1, context, work).order == 1, "numerical series precision propagation")
    var tiny = Budget(0)
    try:
        _ = bernstein.evaluate(half, context, tiny)
        c.check(False, "numerical basis ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "numerical basis budget")


def main() raises:
    var c = Checker("numeric_series_bases")
    var work = Budget()
    domains[FloatDomain](ArithmeticContext(), c, work)
    domains[ComplexDomain](ComplexContext(ArithmeticContext()), c, work)
    domains[BallDomain](BallContext(128), c, work)
    domains[ComplexBallDomain](BallContext(128), c, work)
    c.finish()
