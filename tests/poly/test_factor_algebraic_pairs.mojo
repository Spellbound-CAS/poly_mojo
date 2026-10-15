"""Conjugates sharing a real coordinate and factors with irrational coefficients."""
from poly import dense, algebraic, coefficients, exact
from std.testing import assert_equal, assert_true


def main() raises:
    # All four roots have real coordinate zero; sorted neighbors are not pairs.
    var p = dense.asarray(coefficients((algebraic.asreal(4), 0, 5, 0, 1)))
    var f = dense.factor_checked(p)
    assert_true(exact.check(f))
    var fs = exact.factors(f)
    assert_equal(len(fs), 2)
    assert_equal(fs[0][0].to_list(), dense.asarray(coefficients((algebraic.asreal(1), 0, 1))).to_list())
    assert_equal(fs[1][0].to_list(), dense.asarray(coefficients((algebraic.asreal(4), 0, 1))).to_list())
    var a = dense.root_between([-2, 0, 1], 1, 2)
    # (x - sqrt(2))^2 + 1 has discriminant -4 and stays irreducible over AA.
    var q = dense.asarray(coefficients((3, algebraic.multiply(-2, a), 1)))
    assert_equal(dense.factor(q)[1][0][0].to_list(), q.to_list())
    # A rational cubic splits further over AA; the quadratic has discriminant -3.
    var cubic = dense.factor_checked(coefficients((algebraic.asreal(-1), 0, 0, 1)))
    assert_true(exact.check(cubic))
    fs = exact.factors(cubic)
    assert_equal(len(fs), 2)
    assert_equal(fs[0][0].to_list(), dense.asarray(coefficients((algebraic.asreal(-1), 1))).to_list())
    assert_equal(fs[1][0].to_list(), dense.asarray(coefficients((algebraic.asreal(1), 1, 1))).to_list())
    var complex_cubic = dense.factor_checked(coefficients((algebraic.ascomplex(-1), 0, 0, 1)))
    assert_true(exact.check(complex_cubic))
    var cs = exact.factors(complex_cubic)
    assert_equal(len(cs), 3)
    for entry in cs:
        assert_equal(len(entry[0]), 2)
        assert_equal(entry[1], 1)
    # An irreducible rational cubic acquires genuinely irrational AA factors.
    var rho = dense.root_between([1, 1, 0, 1], -1, 0)
    var general = dense.factor_checked(coefficients((algebraic.asreal(1), 1, 0, 1)))
    assert_true(exact.check(general))
    fs = exact.factors(general)
    assert_equal(len(fs), 2)
    assert_equal(fs[0][0].to_list(), dense.asarray(coefficients((algebraic.negative(rho), 1))).to_list())
    assert_equal(fs[1][0].to_list(), dense.asarray(coefficients((algebraic.add(algebraic.multiply(rho, rho), 1), rho.copy(), 1))).to_list())
    print("array algebraic conjugate factor checks passed")
