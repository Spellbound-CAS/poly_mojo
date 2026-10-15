"""Numerical fraction outputs remain finite and honor shared work limits."""
from poly import rational
from apn_mojo import Rational, FloatFormat, ArithmeticContext, BallContext, ComplexContext, ExactComplex, contains
from polynomials import Budget, PolyErrorKind
from std.testing import assert_true, assert_equal

def main() raises:
    var context = ArithmeticContext(format=FloatFormat(precision=8, emin=-2, emax=2))
    var f = rational.asarray([3], [Rational(1, 2)], context=context)
    var failed = False
    try:
        _ = rational.evaluate(f, 0, context=context)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var balls = BallContext(precision=64)
    var normal = rational.asarray([1], [1, 1], context=balls)
    var budget = Budget()
    _ = rational.evaluate(normal, 1, budget, context=balls)
    var spent = UInt64.MAX - budget.remaining()
    var enough = Budget(spent)
    _ = rational.evaluate(normal, 1, enough, context=balls)
    assert_equal(enough.remaining(), UInt64(0))
    failed = False
    try:
        _ = rational.evaluate(normal, 1, context=balls, max_steps=spent-1)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var complex_context = ComplexContext(ArithmeticContext(format=FloatFormat(precision=64)))
    var complex_fraction = rational.asarray([ExactComplex(1)], [ExactComplex(0, 1)], context=complex_context)
    var complex_value = rational.evaluate(complex_fraction, 0, context=complex_context)
    assert_equal(complex_value.real().to_rational_exact(), Rational(0))
    assert_equal(complex_value.imag().to_rational_exact(), Rational(-1))
    var complex_balls = rational.asarray([ExactComplex(1)], [ExactComplex(0, 1)], context=balls)
    var ball_value = rational.evaluate(complex_balls, 0, context=balls)
    assert_true(contains(ball_value.real(), Rational(0)))
    assert_true(contains(ball_value.imag(), Rational(-1)))
    print("completion numerical limit checks passed")
