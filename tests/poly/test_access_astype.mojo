"""Explicit numerical destinations and re-rounding retain certified conversion rules."""
from poly import dense, sparse, algebraic, coefficients
from apn_mojo import Float, Complex, Ball, ComplexBall, Rational, ArithmeticContext, ComplexContext, BallContext, FloatFormat, contains, Batch, Integer
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var bc = BallContext(precision=71)
    var fc = ArithmeticContext(format=FloatFormat(precision=17))
    var cc = ComplexContext(fc)
    var z = dense.astype[ComplexBall]([1, 2], context=bc)
    comptime assert type_of(z[0]) == ComplexBall
    var r = dense.astype[Ball](z, context=bc)
    comptime assert type_of(r[0]) == Ball
    assert_true(contains(r[1], 2))
    assert_equal(dense.astype[Float](r, context=fc)[1].to_rational_exact(), Rational(2))
    assert_equal(dense.astype[Complex]([1, 2], context=cc)[1].real().precision(), 17)
    var high = dense.asarray([Rational(1, 3)], context=ArithmeticContext(format=FloatFormat(precision=103)))
    var low = dense.astype[Float](high, context=fc)
    assert_equal(low[0].precision(), 17)
    assert_equal(high[0].precision(), 103)
    assert_equal(dense.asarray(high, context=fc)[0].precision(), 103)
    assert_equal(low[0], Float(high[0].to_rational_exact(), context=fc))
    var mixed = dense.astype[ComplexBall]((1, Rational(1, 3), Complex(0, 1)), context=bc)
    assert_true(contains(mixed[1].real(), Rational(1, 3)))
    assert_true(contains(mixed[2].imag(), 1))
    var sp = sparse.astype[ComplexBall]((((2, 0), 1), ((0, 1), Rational(1, 3)), ((2, 0), Float(2, context=fc))), context=bc)
    assert_equal(len(sp), 2)
    assert_true(contains(sp[1][1].real(), 3))
    var real_sp = sparse.astype[Ball](sp, context=bc)
    comptime assert type_of(real_sp[0][1]) == Ball
    assert_true(contains(real_sp[1][1], 3))
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var algebraic_z = algebraic.ascomplex(a, 0)
    assert_true(contains(dense.astype[Ball]([algebraic_z.copy()], context=bc)[0] * dense.astype[Ball]([a.copy()], context=bc)[0], 2))
    assert_equal(len(dense.astype[ComplexBall](List[Int](), context=bc)), 0)
    assert_equal(len(sparse.astype[Ball](List[Tuple[Int, Int]](), context=bc)), 0)
    var inputs = Batch[Integer]([1, 2, 3])
    assert_equal(dense.astype[Float](inputs, context=fc)[2].to_rational_exact(), Rational(3))
    var shaped = inputs.reshape([1, 3])
    var rank_rejected = False
    try: _ = dense.astype[Float](shaped, context=fc)
    except e: rank_rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rank_rejected)
    var sparse_high = sparse.astype[Float]([(2, Rational(1, 3))], context=ArithmeticContext(format=FloatFormat(precision=103)))
    var sparse_low = sparse.astype[Float](sparse_high, context=fc)
    assert_equal(sparse_low[0][1].precision(), 17)
    assert_equal(sparse_high[0][1].precision(), 103)
    assert_equal(len(dense.astype[Ball](Tuple(), context=bc)), 0)
    var nonreal = Complex(1, 1)
    var uncertain = Ball(0, 1, precision=71)
    var nonreal_ball = ComplexBall(1, 1)
    var automatic = BallContext()
    for mode in range(4):
        var rejected = False
        try:
            if mode == 0: _ = dense.astype[Ball]([nonreal.copy()], context=bc)
            elif mode == 1: _ = dense.astype[Float]([uncertain.copy()], context=fc)
            elif mode == 2: _ = sparse.astype[Ball]([(0, nonreal_ball.copy())], context=bc)
            else: _ = dense.astype[Ball]([1], context=automatic)
        except e: rejected = e.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    var budget = Budget()
    _ = dense.astype[ComplexBall]([1, 2], budget, context=bc)
    var used = UInt64.MAX - budget.remaining()
    var exact = Budget(used)
    _ = dense.astype[ComplexBall]([1, 2], exact, context=bc)
    assert_equal(exact.remaining(), UInt64(0))
    var limited = False
    try: _ = dense.astype[ComplexBall]([1, 2], context=bc, max_steps=used - 1)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    print("access explicit dtype checks passed")
