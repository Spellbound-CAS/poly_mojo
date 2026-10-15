"""Convert x squared, evaluate it, and subdivide its Bernstein representation."""
from poly import basis
from apn_mojo import Rational
from std.testing import assert_equal


def main() raises:
    var chebyshev = basis.from_polynomial([0, 0, 1], basis="chebyshev")
    assert_equal(chebyshev[0], Rational(1, 2))
    assert_equal(chebyshev[2], Rational(1, 2))
    var bernstein = basis.asarray([0, 0, 1], basis="bernstein")
    var halves = basis.subdivide(bernstein, Rational(1, 2))
    assert_equal(basis.evaluate(halves[0], 1), Rational(1, 4))
    assert_equal(basis.evaluate(halves[1], 0), Rational(1, 4))
    var slope = basis.derivative(bernstein)
    assert_equal(basis.evaluate(slope, 1), Rational(2))
    print("Value at one half:", basis.evaluate(bernstein, Rational(1, 2)))
    print("Subdivision degree:", len(halves[0]) - 1)
    print("Derivative at one:", basis.evaluate(slope, 1))
