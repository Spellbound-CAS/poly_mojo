"""Factor a rational polynomial and verify its reconstruction."""
from std.testing import assert_equal, assert_true
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, factor, check_factorization


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var polynomial = DensePoly[QQ].from_sparse(ring.parse("(x - 1)^2*(x + 2)", work), work, variable=Int64(0))
    var result = factor(polynomial, work)
    var verified = check_factorization(polynomial, result, work)
    assert_equal(len(result.factors), 2)
    assert_true(verified)
    print("Distinct factors:", len(result.factors))
    print("Product checked:", verified)
