"""Divide two rational polynomials and check the quotient and remainder."""
from std.testing import assert_equal, assert_true
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, divmod, check_division


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var numerator = DensePoly[QQ].from_sparse(ring.parse("x^2 + 1", work), work, variable=Int64(0))
    var divisor = DensePoly[QQ].from_sparse(ring.parse("x - 1", work), work, variable=Int64(0))
    var result = divmod(numerator, divisor, work)
    var expected = DensePoly[QQ].from_sparse(ring.parse("x + 1", work), work, variable=Int64(0))
    assert_equal(result.quotient, expected)
    assert_equal(result.remainder.coefficients[0].value, Rational(2))
    assert_true(check_division(numerator, divisor, result, work))
    print("Quotient is x + 1:", result.quotient == expected)
    print("Remainder:", result.remainder.coefficients[0].value)
