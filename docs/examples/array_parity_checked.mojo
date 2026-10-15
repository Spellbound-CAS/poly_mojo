"""Verify an exact division and a GCD from ascending coefficients."""
from apn_mojo import Rational
from poly import dense, exact
from std.testing import assert_equal, assert_true


def main() raises:
    var division = dense.divide_checked([1, 0, 1], [1, 2])
    assert_true(exact.check(division))
    var quotient = exact.quotient(division)
    assert_equal(quotient[0], Rational(-1, 4))
    assert_equal(quotient[1], Rational(1, 2))
    assert_equal(exact.remainder(division)[0], Rational(5, 4))
    var gcd = dense.gcd_checked([-1, 0, 1], [-1, 1])
    assert_true(exact.check(gcd))
    assert_equal(exact.value(gcd).to_list(), dense.asarray([-1, 1]).to_list())
    print("Exact division and GCD verified")
