"""Shared calculus schedules across rounded and enclosing numerical domains."""
from apn_mojo import Rational, ArithmeticContext, ComplexContext, BallContext, Complex
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, NumericDensePoly, NumericSparsePoly, MonomialOrder, NumericBox, numeric_value, numeric_bounds, numeric_from_dense, numeric_integrate_dense, numeric_definite_integral, numeric_compose_dense, numeric_affine_transform, numeric_jet_evaluate, numeric_taylor_shift, numeric_integrate_sparse, numeric_substitute_many
from support.check import Checker


def point[D: NumericDomain](value: Int, context: D.Context, mut work: Budget) raises -> D.Scalar:
    return numeric_value[D](NumericBox.point(QQ(Rational(value)), QQ(Rational(0)), work), context, work)


def domains[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    var exact = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1)), QQ(Rational(2)), QQ(Rational(3))], work)
    var poly = numeric_from_dense[D, QQ](exact, context, work)
    var integral = numeric_integrate_dense(poly, context, work)
    var back = integral.derivative(context, work)
    for k in range(3):
        c.check(numeric_bounds[D](back.coefficients[k], work).contains(NumericBox.point(exact.coefficients[k], QQ(Rational(0)), work), work), "antiderivative encloses original coefficients")
    c.check(numeric_bounds[D](numeric_definite_integral(poly, point[D](0, context, work), point[D](2, context, work), context, work), work).contains(NumericBox.point(QQ(Rational(14)), QQ(Rational(0)), work), work), "hand-derived integral")
    var transformed = numeric_affine_transform(poly, point[D](2, context, work), point[D](1, context, work), context, work)
    var expected: List[Int] = [6, 16, 12]
    for k in range(3):
        c.check(numeric_bounds[D](transformed.coefficients[k], work).contains(NumericBox.point(QQ(Rational(expected[k])), QQ(Rational(0)), work), work), "hand-derived affine composition")
    var shift = numeric_taylor_shift(poly, point[D](1, context, work), context, work)
    var jet = numeric_jet_evaluate(poly, point[D](1, context, work), 5, context, work)
    c.check(len(jet) == 5 and numeric_bounds[D](jet[4], work) == NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), work), "zero padded numerical jet")
    for k in range(3):
        c.check(numeric_bounds[D](shift.coefficients[k], work) == numeric_bounds[D](jet[k], work), "shift and jet schedules agree")
    var order = MonomialOrder(MonomialOrder.LEX, [0, 1], work)
    var sparse = NumericSparsePoly[D].from_dense(poly, order, context, work)
    var sparse_back = numeric_integrate_sparse(sparse, 1, context, work).derivative(1, context, work)
    var assignments: List[Tuple[Int64, D.Scalar]] = [(Int64(0), point[D](2, context, work)), (Int64(1), point[D](3, context, work))]
    c.check(numeric_bounds[D](sparse_back.evaluate(assignments, context, work), work).contains(NumericBox.point(QQ(Rational(17)), QQ(Rational(0)), work), work), "sparse numerical antiderivative")
    var replacement_dense = NumericDensePoly[D].from_coefficients(0, [point[D](1, context, work), point[D](2, context, work)], work)
    var replacement = NumericSparsePoly[D].from_dense(replacement_dense, order, context, work)
    var substitutions: List[Tuple[Int64, NumericSparsePoly[D]]] = [(Int64(0), replacement.copy())]
    var composed = numeric_substitute_many(sparse, substitutions, context, work).to_dense(0, context, work)
    for k in range(3):
        c.check(numeric_bounds[D](composed.coefficients[k], work).contains(NumericBox.point(QQ(Rational(expected[k])), QQ(Rational(0)), work), work), "sparse simultaneous composition")
    var tiny = Budget(0)
    try:
        _ = numeric_integrate_dense(poly, context, tiny)
        c.check(False, "numerical calculus ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and poly.degree_bound() == 2, "numerical budget preserves polynomial")


def main() raises:
    var c = Checker("numeric_calculus")
    var work = Budget()
    domains[FloatDomain](ArithmeticContext(), c, work)
    domains[ComplexDomain](ComplexContext(ArithmeticContext()), c, work)
    domains[BallDomain](BallContext(128), c, work)
    domains[ComplexBallDomain](BallContext(128), c, work)
    var context = ComplexContext(ArithmeticContext())
    var complex = NumericDensePoly[ComplexDomain].from_coefficients(0, [Complex(0, 1), Complex(1)], work)
    var composed = numeric_compose_dense(complex, complex, context, work)
    c.check(composed.coefficients[0] == Complex(0, 2), "complex imaginary terms retained")
    c.finish()
