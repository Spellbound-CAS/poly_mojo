"""Express a quadratic in Bernstein and Chebyshev bases."""
from std.testing import assert_equal
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, BasisPolynomial


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var polynomial = DensePoly[QQ].from_sparse(ring.parse("x^2", work), work, variable=Int64(0))
    var bernstein = BasisPolynomial[QQ].bernstein(polynomial, 2, work)
    var chebyshev = BasisPolynomial[QQ].chebyshev(polynomial, work)
    assert_equal(bernstein.to_polynomial(work), polynomial)
    assert_equal(chebyshev.to_polynomial(work), polynomial)
    assert_equal(bernstein.coefficients[0].value, Rational(0))
    assert_equal(bernstein.coefficients[1].value, Rational(0))
    assert_equal(bernstein.coefficients[2].value, Rational(1))
    assert_equal(chebyshev.coefficients[0].value, Rational(1, 2))
    assert_equal(chebyshev.coefficients[2].value, Rational(1, 2))
    var value = bernstein.evaluate(QQ(Rational(1, 2)), work)
    assert_equal(value.value, Rational(1, 4))
    print("Bernstein coefficients:", bernstein.coefficients[0].value, bernstein.coefficients[1].value, bernstein.coefficients[2].value)
    print("Value at 1/2:", value.value)
