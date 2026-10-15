"""Integrate a polynomial exactly and evaluate a definite integral."""
from std.testing import assert_equal
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, integrate_dense, derivative_dense, definite_integral, compose_dense


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var polynomial = DensePoly[QQ].from_sparse(ring.parse("3*x^2 + 2*x", work), work, variable=Int64(0))
    var primitive = integrate_dense(polynomial, work)
    assert_equal(derivative_dense(primitive, work), polynomial)
    var area = definite_integral(polynomial, QQ(Rational(0)), QQ(Rational(1)), work)
    assert_equal(area.value, Rational(2))
    var inner = DensePoly[QQ].from_sparse(ring.parse("x + 1", work), work, variable=Int64(0))
    var shifted = compose_dense(polynomial, inner, work)
    var expected = DensePoly[QQ].from_sparse(ring.parse("3*x^2 + 8*x + 5", work), work, variable=Int64(0))
    assert_equal(shifted, expected)
    print("Integral from 0 to 1:", area.value)
    print("Substitution checked:", shifted == expected)
