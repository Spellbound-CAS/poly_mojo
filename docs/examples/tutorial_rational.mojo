"""Decompose a rational function and compute its residue at a pole."""
from std.testing import assert_equal, assert_true
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, PolynomialRing, RationalFunction, partial_fractions, check_partial_fractions, rational_evaluate, rational_residue


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var normalized = RationalFunction.normalize(ring.parse("x^2", work), ring.parse("x - 1", work), work)
    var parts = partial_fractions(normalized[0], 0, work)
    assert_true(check_partial_fractions(normalized[0], 0, parts, work))
    var point: List[Tuple[Int64, QQ]] = [(Int64(0), QQ(Rational(2)))]
    var value = rational_evaluate(normalized[0], point, work)
    var residue = rational_residue(normalized[0], 0, QQ(Rational(1)), work)
    assert_equal(value.value, Rational(4))
    assert_equal(residue.value, Rational(1))
    assert_equal(len(parts.terms), 1)
    print("Value at 2:", value.value)
    print("Residue at 1:", residue.value)
