"""Independent series identities, prepared-state failures and numerical schedules."""
from apn_mojo import Rational, ExactComplex, Ball, BallContext, ComplexBall, Float, Complex, ArithmeticContext, FloatFormat, ComplexContext
from polynomials import *
from polynomials._toolkit_exact import _ExactOps
from support.check import Checker


def series_identities[C: CoefficientField](ctx: C.Context, mut check: Checker, mut work: Budget) raises:
    for n in [31, 32, 33, 65, 128, 256]:
        var coefficients = List[C](length=n, fill=C.one(ctx))
        var f = TruncatedSeries[C].from_polynomial(DensePoly[C].from_coefficients(0, ctx, coefficients^, work), n, work)
        var inverse = f.inverse(work)
        var expected = DensePoly[C].from_coefficients(0, ctx, [C.one(ctx), C.one(ctx).neg(work)], work)
        check.check(inverse.polynomial == expected, "dense geometric inverse across dispatch sizes")
        coefficients = List[C]()
        for k in range(n):
            coefficients.append(_ExactOps[C].integer(1 + k % 5, ctx, work))
        var dense = TruncatedSeries[C].from_polynomial(DensePoly[C].from_coefficients(0, ctx, coefficients^, work), n, work)
        var got = dense.inverse(work)
        # Independent ordinary polynomial multiplication, then inspect x^0..x^(n-1).
        var product = dense.polynomial.mul(got.polynomial, work)
        var valid = product.coefficients[0] == C.one(ctx)
        for k in range(1, min(n, len(product.coefficients))):
            valid = valid and product.coefficients[k].is_zero()
        check.check(valid, "dense inverse multiplication identity")
        check.check(dense.divide(dense, work).polynomial == DensePoly[C].from_coefficients(0, ctx, [C.one(ctx)], work), "dense self quotient")
    for n in [8, 17, 33, 65]:
        var coefficients = List[C](length=n, fill=C.one(ctx))
        coefficients[0] = C.zero(ctx)
        var f = TruncatedSeries[C].from_polynomial(DensePoly[C].from_coefficients(0, ctx, coefficients^, work), n, work)
        var reverted = f.revert(work)
        var valid = len(reverted.polynomial.coefficients) == n
        for k in range(n):
            valid = valid and reverted.polynomial.coefficients[k] == _ExactOps[C].integer(0 if k == 0 else (1 if k % 2 else -1), ctx, work)
        check.check(valid, "x/(1-x) reverts to x/(1+x)")
        var tiny = Budget(17)
        var saved = f.polynomial.copy()
        try:
            _ = f.revert(tiny)
            check.check(False, "series reversion limit ignored")
        except e:
            check.check(e.kind == PolyErrorKind.LIMIT and f.polynomial == saved, "series failure preserves borrowed input")


def prepared[C: CoefficientField](ctx: C.Context, mut check: Checker, mut work: Budget) raises:
    var points = List[C]()
    var values = List[C]()
    for k in range(12):
        points.append(_ExactOps[C].integer(k, ctx, work))
        values.append(_ExactOps[C].integer(k * k + 1, ctx, work))
    var plan = MultipointPlan[C].prepare_interpolation(-7, ctx, points, work)
    var saved = plan.copy()
    var tiny = Budget(29)
    try:
        _ = plan.interpolate(values, tiny)
        check.check(False, "interpolation limit ignored")
    except e:
        check.check(e.kind == PolyErrorKind.LIMIT, "interpolation preparation is budgeted")
    var first_work = Budget()
    var first = plan.interpolate(values, first_work)
    var cold_work = Budget()
    var cold = saved.interpolate(values, cold_work)
    var warm_work = Budget()
    var warm = plan.interpolate(values, warm_work)
    var ordinary = MultipointPlan[C].prepare(-7, ctx, points, work)
    var ordinary_work = Budget()
    _ = ordinary.interpolate(values, ordinary_work)
    check.check(warm_work.remaining() > ordinary_work.remaining(), "prepared inverse weights remove repeated work")
    check.check(first == cold and warm == cold and check_interpolation(warm, points, values, work), "prepared interpolation survives failure and copy")
    check.check(plan.evaluate(first, work) == values, "interpolation retains evaluation state")
    values[0] = values[0].add(C.one(ctx), work)
    var changed = plan.interpolate(values, work)
    check.check(check_interpolation(changed, points, values, work), "cached geometry does not cache values")
    var duplicates: List[C] = [C.one(ctx), C.one(ctx)]
    try:
        _ = MultipointPlan[C].prepare_interpolation(0, ctx, duplicates, work)
        check.check(False, "duplicate prepared interpolation accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "prepared interpolation rejects duplicates")
    tiny = Budget(0)
    try:
        _ = plan.interpolate(values, tiny)
        check.check(False, "warm interpolation limit ignored")
    except e:
        check.check(e.kind == PolyErrorKind.LIMIT, "warm plan still charges work")


def main() raises:
    var check = Checker("toolkit_review")
    var work = Budget()
    series_identities[QQ](NoContext(), check, work)
    series_identities[QQi](NoContext(), check, work)
    var prime = Modulus(65521, work)
    series_identities[ModInt](prime, check, work)
    # Newton reversion needs only the linear inverse, even in small characteristic.
    series_identities[ModInt](Modulus(7, work), check, work)
    prepared[QQ](NoContext(), check, work)
    prepared[QQi](NoContext(), check, work)
    prepared[ModInt](prime, check, work)
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1)), QQ(Rational(2)), QQ(Rational(3))], work)
    var ctx = BallContext(16)
    var numeric = numeric_from_dense[BallDomain, QQ](f, ctx, work)
    var basis = NumericBasisPolynomial[BallDomain].bernstein(numeric, 9, ctx, work)
    for k in [-1, 0, 1, 2]:
        var point = Ball(Rational(k, 3), Rational(1, 1000))
        var split = basis.subdivide(point, ctx, work)
        var evaluated = basis.evaluate(point, ctx, work)
        check.check(numeric_bounds[BallDomain](evaluated, work) == numeric_bounds[BallDomain](split[0].coefficients[9], work), "Bernstein evaluation retains subdivision arithmetic schedule")
    check.finish()
