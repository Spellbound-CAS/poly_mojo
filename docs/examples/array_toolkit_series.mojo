"""Keep a formal truncation order through reciprocals and exponentials."""
from poly import series
from apn_mojo import Rational
from std.testing import assert_equal


def main() raises:
    var denominator = series.asarray([1, -1], order=6)
    var geometric = series.reciprocal(denominator)
    var product = series.multiply(denominator, geometric)
    assert_equal(product[0], Rational(1))
    assert_equal(product[5], Rational(0))
    var exponential = series.exp([0, 1], order=6)
    assert_equal(exponential[5], Rational(1, 120))
    print("Known coefficients:", len(geometric))
    print("Coefficient of x^5:", exponential[5])
    print("Derivative order:", series.derivative(exponential).order())
