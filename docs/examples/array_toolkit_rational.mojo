"""Normalize a rational expression and inspect a pole using exact coefficients."""
from poly import rational
from apn_mojo import Rational
from std.testing import assert_equal, assert_true


def main() raises:
    var simplified = rational.asarray([-1, 0, 1], [-1, 1])
    assert_true(rational.check(simplified))
    assert_equal(rational.evaluate(simplified, 1), Rational(2))
    var function = rational.asarray([1], [0, 1, -1])
    var parts = rational.partial_fractions(function)
    assert_true(rational.check(parts))
    var local = rational.laurent(function, 0, terms=4)
    assert_true(rational.check(local))
    assert_equal(local.valuation, -1)
    assert_equal(local.coefficients[0], Rational(1))
    var quartic = rational.partial_fractions(rational.asarray([1], [1, 0, 0, 0, 1]))
    assert_true(rational.check(quartic))
    assert_equal(len(quartic.terms), 1)
    assert_equal(quartic.terms[0][1].to_list(), [Rational(1), 0, 0, 0, 1])
    print("Reduced value at one:", rational.evaluate(simplified, 1))
    print("Proper fraction terms:", len(parts.terms))
    print("Lowest Laurent exponent:", local.valuation)
    print("Residue at zero:", rational.residue(function, 0))
    print("Irreducible quartic denominator checked")
