"""Interpolation in all numerical domains and separated interval families."""
from apn_mojo import Rational, ArithmeticContext, ComplexContext, BallContext
from polynomials import Budget, PolyErrorKind, QQ, NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, NumericBox, NumericDensePoly, numeric_value, numeric_bounds, numeric_multipoint_evaluate, numeric_interpolate, numeric_hermite_interpolate
from support.check import Checker


def point[D: NumericDomain](value: Int, context: D.Context, mut work: Budget) raises -> D.Scalar:
    return numeric_value[D](NumericBox.point(QQ(Rational(value)), QQ(Rational(0)), work), context, work)


def domains[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    var points: List[D.Scalar] = [point[D](0, context, work), point[D](1, context, work), point[D](2, context, work)]
    var values: List[D.Scalar] = [point[D](1, context, work), point[D](4, context, work), point[D](9, context, work)]
    var poly = numeric_interpolate[D](0, points, values, context, work)
    var expected: List[Int] = [1, 2, 1]
    for i in range(3):
        c.check(numeric_bounds[D](poly.coefficients[i], work).contains(NumericBox.point(QQ(Rational(expected[i])), QQ(Rational(0)), work), work), "hand-derived interpolation coefficients")
    var evaluated = numeric_multipoint_evaluate(poly, points, context, work)
    for i in range(3):
        c.check(numeric_bounds[D](evaluated[i], work).contains(numeric_bounds[D](values[i], work), work), "multipoint interpolation samples")
    _ = points.pop()
    var jets: List[List[D.Scalar]] = [[point[D](1, context, work), point[D](2, context, work)], [point[D](4, context, work), point[D](5, context, work)]]
    var hermite = numeric_hermite_interpolate[D](0, points, jets, context, work)
    var wanted: List[Int] = [1, 2, 0, 1]
    for i in range(4):
        c.check(numeric_bounds[D](hermite.coefficients[i], work).contains(NumericBox.point(QQ(Rational(wanted[i])), QQ(Rational(0)), work), work), "numerical Hermite coefficients")
    points[1] = points[0].copy()
    _ = values.pop()
    try:
        _ = numeric_interpolate[D](0, points, values, context, work)
        c.check(False, "repeated numeric nodes accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "repeated numeric nodes rejected")
    var tiny = Budget(0)
    try:
        _ = numeric_hermite_interpolate[D](0, points, jets, context, tiny)
        c.check(False, "numeric interpolation ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "numeric interpolation budget")


def main() raises:
    var c = Checker("numeric_interpolation")
    var work = Budget()
    domains[FloatDomain](ArithmeticContext(), c, work)
    domains[ComplexDomain](ComplexContext(ArithmeticContext()), c, work)
    domains[BallDomain](BallContext(128), c, work)
    domains[ComplexBallDomain](BallContext(128), c, work)
    # Every line through (0,[1,2]) and (1,[3,4]) has intercept [1,2] and slope [1,3].
    var context = BallContext(128)
    var zero = QQ(Rational(0))
    var points: List[BallDomain.Scalar] = [point[BallDomain](0, context, work), point[BallDomain](1, context, work)]
    var values: List[BallDomain.Scalar] = [numeric_value[BallDomain](NumericBox(QQ(Rational(1)), QQ(Rational(2)), zero.copy(), zero.copy()), context, work), numeric_value[BallDomain](NumericBox(QQ(Rational(3)), QQ(Rational(4)), zero.copy(), zero.copy()), context, work)]
    var family = numeric_interpolate[BallDomain](0, points, values, context, work)
    c.check(numeric_bounds[BallDomain](family.coefficients[0], work).contains(NumericBox(QQ(Rational(1)), QQ(Rational(2)), zero.copy(), zero.copy()), work), "family intercept enclosure")
    c.check(numeric_bounds[BallDomain](family.coefficients[1], work).contains(NumericBox(QQ(Rational(1)), QQ(Rational(3)), zero.copy(), zero.copy()), work), "family slope enclosure")
    c.finish()
