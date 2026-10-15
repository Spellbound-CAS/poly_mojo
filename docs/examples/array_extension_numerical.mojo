"""Recover numerical coefficients and inspect checked division."""
from apn_mojo import ArithmeticContext, FloatFormat, Rational
from poly import dense, sparse
from std.testing import assert_equal, assert_true


def main() raises:
    var context = ArithmeticContext(format=FloatFormat(precision=96))
    var p = dense.interpolate([-1, 0, 1], [0, 1, 4], context=context)
    assert_equal(p[0].to_rational_exact(), Rational(1))
    assert_equal(p[1].to_rational_exact(), Rational(2))
    assert_equal(p[2].to_rational_exact(), Rational(1))
    var samples = dense.evaluate_many(p, [-1, 0, 1, -1], context=context)
    assert_equal(len(samples), 4)
    assert_equal(samples[3].to_rational_exact(), Rational(0))
    var divisor = dense.asarray([1, 1], context=context)
    var division = dense.divide_checked(p, divisor, context=context)
    assert_true(dense.check_division(division))
    assert_equal(division.quotient()[1].to_rational_exact(), Rational(1))
    assert_equal(len(division.remainder()), 0)
    var surface = sparse.asarray([((1, 1), 2)], context=context)
    var integrated = sparse.definite_integral(surface, 0, 2, axis=0, context=context)
    assert_equal(sparse.evaluate(integrated, (0, 3), context=context).to_rational_exact(), Rational(12))
    print("Recovered (x + 1)^2 and retained four samples")
    print("Numerical division checked; integrated surface is 4*y")
