"""Evaluate a polynomial with rounded arithmetic and check an interval enclosure."""
from std.testing import assert_equal, assert_true
from apn_mojo import Rational, Float, FloatFormat, ArithmeticContext, Ball, BallContext
from polynomials import Budget, QQ, NoContext, DensePoly, NumericDensePoly, FloatDomain, BallDomain, enclose_dense, check_enclose_dense


def main() raises:
    var work = Budget()
    var context = ArithmeticContext(format=FloatFormat(128))
    var polynomial = NumericDensePoly[FloatDomain].from_coefficients(0, [Float(-2), Float(0), Float(1)], work)
    var value = polynomial.evaluate(Float("1.5"), context, work)
    assert_equal(value, Float("0.25"))
    print("Rounded value at 1.5:", value)

    var exact = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], work)
    var argument = Ball(Rational(3, 2), Rational(1, 16))
    var enclosure = enclose_dense[BallDomain](exact, argument, BallContext(128), work)
    var verified = check_enclose_dense[BallDomain](exact, argument, enclosure, work)
    assert_true(verified)
    print("Enclosure checked:", verified)
