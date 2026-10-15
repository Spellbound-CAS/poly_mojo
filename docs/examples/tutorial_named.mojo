"""Build a named polynomial and evaluate its value and derivative exactly."""
from std.testing import assert_equal
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, PolynomialRing, evaluate, derivative


def main() raises:
    var work = Budget(1000000)
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var polynomial = ring.parse("(x + 1)^3", work)
    var point: List[Tuple[Int64, QQ]] = [(Int64(0), QQ(Rational(3)))]
    var value = evaluate(polynomial, point, work)
    var slope = evaluate(derivative(polynomial, 0, work), point, work)
    assert_equal(value.value, Rational(64))
    assert_equal(slope.value, Rational(48))
    print("Value at 3:", value.value)
    print("Derivative at 3:", slope.value)
