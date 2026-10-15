"""Enclose formal coefficients and work directly in a Bernstein basis."""
from apn_mojo import BallContext, Rational, contains
from poly import series, basis
from std.testing import assert_equal, assert_true


def main() raises:
    var context = BallContext(precision=96)
    var x = series.asarray([0, 1], order=5, context=context)
    var exponential = series.exp(x, context=context)
    assert_equal(exponential.order(), 5)
    assert_true(contains(exponential[4], Rational(1, 24)))
    var curve = basis.from_polynomial([0, 0, 1], basis="bernstein", degree=3, context=context)
    assert_equal(len(curve), 4)
    assert_true(contains(basis.evaluate(curve, Rational(1, 3), context=context), Rational(1, 9)))
    var halves = basis.subdivide(curve, Rational(1, 2), context=context)
    assert_true(contains(basis.evaluate(halves[0], 1, context=context), Rational(1, 4)))
    print("exp(x): five known coefficients with enclosing balls")
    print("Bernstein evaluation and subdivision enclose the expected values")
