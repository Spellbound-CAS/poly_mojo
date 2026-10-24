"""Certainly rounded conversions agree with independent square-root arithmetic."""
from poly import dense, sparse, algebraic, coefficients
from apn_mojo import Integer, Float, Ball, BallContext, ArithmeticContext, ComplexContext, FloatFormat, Rational, RoundingMode, sqrt, contains
from polynomials import Budget, PolyErrorKind, QQ
from std.testing import assert_equal, assert_true


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    for mode in [RoundingMode.nearest_even, RoundingMode.toward_zero, RoundingMode.toward_positive, RoundingMode.toward_negative]:
        var context = ArithmeticContext(format=FloatFormat(precision=24), rounding=mode)
        var converted = dense.asarray([a.copy()], context=context)
        assert_equal(converted[0], sqrt(Float(2), context=context))
        assert_equal(converted[0].precision(), 24)
    var bc = BallContext(precision=80)
    var fc = ArithmeticContext(format=FloatFormat(precision=80))
    var p = dense.asarray(coefficients((a.copy(), 0, 1)))
    var ball = dense.asarray(p, context=bc)
    assert_true(contains(ball[0] * ball[0], Rational(2)))
    var z = algebraic.ascomplex(a, 1)
    var cc = ComplexContext(real=fc, imag=ArithmeticContext(format=FloatFormat(precision=32)))
    var complex = dense.asarray([z.copy()], context=cc)
    assert_equal(complex[0].real(), sqrt(Float(2), context=fc))
    assert_equal(complex[0].imag().precision(), 32)
    var complex_ball = dense.asarray((1, z.copy()), context=bc)
    assert_true(contains(complex_ball[1].imag(), Rational(1)))
    var sp = sparse.asarray([((1, 0), z.copy())], context=bc)
    assert_equal(len(sp), 1)
    var empty = dense.asarray(List[algebraic.Complex](), context=bc)
    assert_equal(len(empty), 0)
    var real_complex = algebraic.ascomplex(a, 0)
    assert_equal(dense.asarray([real_complex.copy()], context=fc)[0], sqrt(Float(2), context=fc))
    var bad = False
    try:
        _ = dense.asarray([z.copy()], context=fc)
    except e:
        bad = e.kind == PolyErrorKind.DOMAIN
    assert_true(bad)
    # The exact value 1 can have a nonminimal irrational-style root witness.
    var budget = Budget()
    var witness = a.witness().copy()
    witness.polynomial = dense.to_core([Rational(-1), Rational(1), Rational(-1), Rational(1)])
    witness.low = QQ(Rational(0))
    witness.high = QQ(Rational(2))
    var exact_one = algebraic.asreal(witness)
    var trapped = ArithmeticContext(format=FloatFormat(precision=24), trap_inexact=True)
    assert_equal(dense.asarray([exact_one.copy()], context=trapped)[0], Float(1))
    bad = False
    try:
        _ = dense.asarray([a.copy()], context=trapped)
    except e:
        bad = e.kind == PolyErrorKind.DOMAIN
    assert_true(bad)
    var ceiling_context = ArithmeticContext(format=FloatFormat(precision=80), max_precision=2)
    bad = False
    try:
        _ = dense.asarray([a.copy()], context=ceiling_context)
    except e:
        bad = e.kind == PolyErrorKind.DOMAIN
    assert_true(bad)
    var work = Budget()
    _ = dense.asarray([a.copy()], work, context=fc)
    var spent = UInt64.MAX - work.remaining()
    var exact = Budget(spent)
    _ = dense.asarray([a.copy()], exact, context=fc)
    assert_equal(exact.remaining(), UInt64(0))
    bad = False
    try:
        _ = dense.asarray([a.copy()], context=fc, max_steps=spent - 1)
    except e:
        bad = e.kind == PolyErrorKind.LIMIT
    assert_true(bad)
    assert_true(algebraic.check(a))
    var scale = Integer(1) << 200
    var tiny_square = Rational(2, scale * scale)
    var tiny = dense.root_between([-tiny_square, Rational(0), Rational(1)], 0, 1)
    var narrow = ArithmeticContext(format=FloatFormat(precision=40), max_precision=64)
    assert_equal(dense.asarray([tiny.copy()], context=narrow)[0], sqrt(Float(tiny_square, context=ArithmeticContext(format=FloatFormat(precision=512))), context=narrow))
    var tiny_ball = dense.asarray([tiny.copy()], context=BallContext(precision=40))
    assert_true(not contains(tiny_ball[0], Rational(0)))
    var negative = algebraic.negative(tiny)
    assert_equal(dense.asarray([negative.copy()], context=narrow)[0], -dense.asarray([tiny.copy()], context=narrow)[0])
    print("array algebraic numerical conversion checks passed")
