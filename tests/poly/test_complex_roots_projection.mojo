"""Nonquadratic elimination and roots sharing real or imaginary coordinates."""
from apn_mojo import Rational
from poly import dense, roots, algebraic
from std.testing import assert_equal, assert_true


def main() raises:
    var s3 = dense.root_between([-3, 0, 1], 1, 2)
    var height = algebraic.divide(s3, 2)
    var cubic = dense.complex_roots([-1, 0, 0, 1])
    assert_equal(len(cubic), 3)
    assert_true(cubic[0] == algebraic.ascomplex(Rational(-1, 2), algebraic.negative(height)))
    assert_true(cubic[1] == algebraic.ascomplex(Rational(-1, 2), height))
    assert_true(cubic[2] == algebraic.ascomplex(1))
    assert_true(roots.check(cubic))
    var quartic = dense.complex_roots([-1, 0, 0, 0, 1])
    assert_equal(len(quartic), 4)
    assert_true(quartic[0] == algebraic.ascomplex(-1))
    assert_true(quartic[1] == algebraic.ascomplex(0, -1))
    assert_true(quartic[2] == algebraic.ascomplex(0, 1))
    assert_true(quartic[3] == algebraic.ascomplex(1))
    assert_true(roots.check(quartic))
    var irreducible = dense.complex_roots([1, 1, 0, 1])
    assert_equal(len(irreducible), 3)
    for root in irreducible:
        assert_true(dense.evaluate([1, 1, 0, 1], root).is_zero())
    assert_true(roots.check(irreducible))
    print("array complex root projection checks passed")
