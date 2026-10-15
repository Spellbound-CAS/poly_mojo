"""Invert a formal series modulo x to the fourth power."""
from std.testing import assert_equal
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, TruncatedSeries


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var polynomial = DensePoly[QQ].from_sparse(ring.parse("1 - x", work), work, variable=Int64(0))
    var series = TruncatedSeries[QQ].from_polynomial(polynomial, 4, work)
    var reciprocal = series.inverse(work)
    var expected = DensePoly[QQ].from_sparse(ring.parse("1 + x + x^2 + x^3", work), work, variable=Int64(0))
    assert_equal(reciprocal.polynomial, expected)
    var product = series.mul(reciprocal, work)
    var one = DensePoly[QQ].from_sparse(ring.parse("1", work), work, variable=Int64(0))
    assert_equal(product.polynomial, one)
    print("Reciprocal through x^3:", reciprocal.polynomial == expected)
    print("Product modulo x^4 is 1:", product.polynomial == one)
