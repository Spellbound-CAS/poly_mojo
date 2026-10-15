"""Hand-derived arithmetic checks across all four apn_mojo numerical domains."""
from apn_mojo import Float, Complex, Ball, ComplexBall, Rational, FloatFormat, ArithmeticContext, ComplexContext, BallContext, RoundingMode, contains
from polynomials import Budget, QQ, ZZ, ModInt, Modulus, NoContext, Monomial, MonomialOrder, DensePoly, PolyErrorKind
from polynomials.numeric_domains import NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain
from polynomials.numeric_dense import NumericDensePoly
from polynomials.numeric_sparse import NumericSparsePoly, NumericTerm
from polynomials.numeric_bounds import NumericBox, numeric_bounds, numeric_value
from polynomials.numeric_convert import convert_numeric_dense, numeric_from_dense
from support.check import Checker


def sparse[D: NumericDomain](poly: NumericDensePoly[D], point: D.Scalar, context: D.Context, mut c: Checker, mut budget: Budget) raises:
    var order = MonomialOrder()
    var p = NumericSparsePoly[D].from_dense(poly, order, context, budget)
    var assignments: List[Tuple[Int64, D.Scalar]] = [(poly.variable, point.copy())]
    var dense_value = numeric_bounds[D](poly.evaluate(point, context, budget), budget)
    c.check(numeric_bounds[D](p.evaluate(assignments, context, budget), budget).contains(dense_value, budget), "sparse-dense exact fixture evaluation")
    var twice = p.add(p, context, budget)
    var back = twice.sub(p, context, budget).to_dense(poly.variable, context, budget)
    c.check(numeric_bounds[D](back.evaluate(point, context, budget), budget).contains(dense_value, budget), "sparse addition/subtraction")
    var product = p.mul(p, context, budget).to_dense(poly.variable, context, budget)
    c.check(numeric_bounds[D](product.evaluate(point, context, budget), budget).contains(numeric_bounds[D](poly.mul(poly, context, budget).evaluate(point, context, budget), budget), budget), "sparse multiplication")
    var slope = p.derivative(poly.variable, context, budget).to_dense(poly.variable, context, budget)
    c.check(numeric_bounds[D](slope.evaluate(point, context, budget), budget).contains(numeric_bounds[D](poly.derivative(context, budget).evaluate(point, context, budget), budget), budget), "sparse differentiation")
    var minus = p.neg(context, budget).neg(context, budget)
    c.check(numeric_bounds[D](minus.evaluate(assignments, context, budget), budget).contains(dense_value, budget), "sparse double negation")
    var two = numeric_value[D](NumericBox.point(QQ(Rational(2)), QQ(Rational(0)), budget), context, budget)
    var scaled = poly.scale(two, context, budget).scale(two, context, budget, divide=True)
    c.check(numeric_bounds[D](scaled.evaluate(point, context, budget), budget).contains(dense_value, budget), "successful dense scalar division")
    var sparse_scaled = p.scale(two, context, budget).scale(two, context, budget, divide=True)
    c.check(numeric_bounds[D](sparse_scaled.evaluate(assignments, context, budget), budget).contains(dense_value, budget), "successful sparse scalar division")
    var dense_sum = poly.add(poly, context, budget).sub(poly, context, budget)
    c.check(numeric_bounds[D](dense_sum.neg(context, budget).neg(context, budget).evaluate(point, context, budget), budget).contains(dense_value, budget), "dense combine and negation exact fixture")
    var squared = p.pow(2, context, budget).to_dense(poly.variable, context, budget)
    c.check(numeric_bounds[D](squared.evaluate(point, context, budget), budget).contains(numeric_bounds[D](product.evaluate(point, context, budget), budget), budget), "sparse binary powering")


def main() raises:
    var check = Checker("numerical arithmetic")
    var budget = Budget()
    var fc = ArithmeticContext(format=FloatFormat(128))
    var cc = ComplexContext(fc)
    var bc = BallContext(128)
    var f = NumericDensePoly[FloatDomain].from_coefficients(0, [Float(1), Float(2), Float(1)], budget)
    check.check(f.evaluate(Float(3), fc, budget) == Rational(16), "Float Horner")
    check.check(f.mul(f, fc, budget).evaluate(Float(3), fc, budget) == Rational(256), "Float convolution")
    check.check(f.derivative(fc, budget).evaluate(Float(3), fc, budget) == Rational(8), "Float derivative")
    var c = NumericDensePoly[ComplexDomain].from_coefficients(0, [Complex(1), Complex(0, 1)], budget)
    var cv = c.evaluate(Complex(0, 1), cc, budget)
    check.check(cv.is_zero(), "complex i squared")
    var b = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(1), Ball(2), Ball(1)], budget)
    check.check(contains(b.evaluate(Ball(3), bc, budget), Rational(16)), "Ball Horner")
    check.check(contains(b.pow(2, bc, budget).evaluate(Ball(3), bc, budget), Rational(256)), "Ball powering")
    var cb = NumericDensePoly[ComplexBallDomain].from_coefficients(0, [ComplexBall(1), ComplexBall(0, 1)], budget)
    var cbv = cb.evaluate(ComplexBall(0, 1), bc, budget)
    check.check(contains(cbv.real(), Rational(0)) and contains(cbv.imag(), Rational(0)), "complex Ball i squared")
    sparse(f, Float(3), fc, check, budget)
    sparse(c, Complex(0, 1), cc, check, budget)
    sparse(b, Ball(3), bc, check, budget)
    sparse(cb, ComplexBall(0, 1), bc, check, budget)
    var uncertain = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(1), Ball(0, 1)], budget)
    check.check(uncertain.degree_bound() == 1 and not uncertain.certified_degree(budget), "uncertain leading coefficient retained")
    var trimmed = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(1), Ball(0)], budget)
    check.check(trimmed.degree_bound() == 0 and trimmed.certified_degree(budget).value() == 0, "singleton zero trimmed")
    var uncertain_divisor = Ball(0, 1)
    try:
        _ = b.scale(uncertain_divisor, bc, budget, divide=True)
        check.check(False, "uncertain scalar divisor accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "uncertain scalar divisor rejected")
    var one_ball = Ball(1)
    var implicit_precision = BallContext()
    try:
        _ = b.evaluate(one_ball, implicit_precision, budget)
        check.check(False, "implicit precision accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "explicit precision required")
    var tiny = Budget(0)
    try:
        _ = b.mul(b, bc, tiny)
        check.check(False, "exhausted budget accepted")
    except e:
        check.check(e.kind == PolyErrorKind.LIMIT and b.degree_bound() == 2, "exhaustion preserves input")
    var exact = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1, 3)), QQ(Rational(2, 3))], budget)
    var enclosed = numeric_from_dense[BallDomain, QQ](exact, BallContext(8), budget)
    check.check(contains(enclosed.coefficients[0], Rational(1, 3)), "rational conversion encloses without point rounding")
    var complex_enclosed = convert_numeric_dense[ComplexBallDomain, BallDomain](enclosed, bc, budget)
    check.check(contains(complex_enclosed.coefficients[0].real(), Rational(1, 3)), "explicit real to complex enclosure")
    var bounds = NumericBox.real_interval(QQ(Rational(1)), QQ(Rational(2)), budget)
    try:
        _ = numeric_value[FloatDomain](bounds, fc, budget)
        check.check(False, "uncertain ball midpoint silently converted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "uncertain Float conversion rejected")
    var point_bounds = NumericBox.point(QQ(Rational(9, 8)), QQ(Rational(0)), budget)
    var up = ArithmeticContext(format=FloatFormat(2), rounding=RoundingMode.toward_positive)
    check.check(numeric_value[FloatDomain](point_bounds, up, budget) == Rational(3, 2), "directed arbitrary precision rounding")
    var both = b.evaluate_with_derivative(Ball(3), bc, budget)
    check.check(contains(both[0], Rational(16)) and contains(both[1], Rational(8)), "simultaneous Horner derivative")
    var sparse_uncertain = NumericSparsePoly[BallDomain].from_dense(uncertain, MonomialOrder(), bc, budget)
    check.check(sparse_uncertain.degree_bound(budget) == 1 and not sparse_uncertain.certified_degree(budget), "sparse uncertain degree")
    var m = Monomial([(Int64(0), UInt32(1))], budget)
    var duplicates: List[NumericTerm[FloatDomain]] = [NumericTerm[FloatDomain](m.copy(), Float(1)), NumericTerm[FloatDomain](m.copy(), Float(-1)), NumericTerm[FloatDomain](m.copy(), Float(3))]
    var merged = NumericSparsePoly[FloatDomain].from_terms(duplicates^, MonomialOrder(), fc, budget)
    check.check(len(merged.terms) == 1 and merged.terms[0].coefficient == Rational(3), "stable duplicate cancellation and reappearance")
    var cancellation_group: List[NumericTerm[FloatDomain]] = [NumericTerm[FloatDomain](m.copy(), Float(1)), NumericTerm[FloatDomain](m.copy(), Float(-1)), NumericTerm[FloatDomain](m.copy(), Float(Rational(9, 8)))]
    var low_precision = ArithmeticContext(format=FloatFormat(2))
    var cancellation_sum = NumericSparsePoly[FloatDomain].from_terms(cancellation_group^, MonomialOrder(), low_precision, budget)
    check.check(cancellation_sum.terms[0].coefficient == Rational(1), "intermediate duplicate zero retains the prescribed rounding step")
    var nonfinite = Float.infinity()
    try:
        _ = NumericDensePoly[FloatDomain].from_coefficients(0, [nonfinite.copy()], budget)
        check.check(False, "infinite coefficient accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "infinite coefficient rejected")
    var nan = Float.nan()
    try:
        _ = f.evaluate(nan, fc, budget)
        check.check(False, "NaN argument accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "NaN argument rejected")
    var trapped = ArithmeticContext(format=FloatFormat(2), trap_inexact=True)
    try:
        _ = numeric_value[FloatDomain](point_bounds, trapped, budget)
        check.check(False, "inexact trap ignored")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "context inexact trap mapped")
    var modular_zero = DensePoly[ModInt](0, Modulus(7, budget))
    try:
        _ = numeric_from_dense[BallDomain, ModInt](modular_zero, bc, budget)
        check.check(False, "zero modular polynomial embedded")
    except e:
        check.check(e.kind == PolyErrorKind.UNSUPPORTED, "modular domain rejected even for zero")
    var missing = List[Tuple[Int64, Float]]()
    try:
        _ = merged.evaluate(missing, fc, budget)
        check.check(False, "missing numerical assignment accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "missing assignment rejected")
    var duplicate_assignments: List[Tuple[Int64, Float]] = [(Int64(0), Float(1)), (Int64(0), Float(2))]
    try:
        _ = merged.evaluate(duplicate_assignments, fc, budget)
        check.check(False, "duplicate numerical assignment accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "duplicate assignment rejected")
    var certain_interval = NumericBox.real_interval(QQ(Rational(71, 64)), QQ(Rational(73, 64)), budget)
    check.check(numeric_value[FloatDomain](certain_interval, low_precision, budget) == Rational(1), "public certainty conversion accepts a uniformly rounded interval")
    var zero_crossing = NumericBox.real_interval(QQ(Rational(-1, 1048576)), QQ(Rational(1, 1048576)), budget)
    var underflow_context = ArithmeticContext(format=FloatFormat(2, emin=-2, emax=2))
    try:
        _ = numeric_value[FloatDomain](zero_crossing, underflow_context, budget)
        check.check(False, "opposite signed underflow zeros accepted as certain")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "certainty preserves signed-zero and status distinctions")
    var wide = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(1, Rational(1, 4))], budget)
    var difference = wide.sub(wide, bc, budget)
    check.check(difference.degree_bound() == 0 and contains(difference.coefficients[0], Rational(1, 2)), "uncertain cancellation retains independent coefficient choices")
    check.finish()
