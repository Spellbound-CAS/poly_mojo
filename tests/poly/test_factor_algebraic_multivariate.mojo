"""Exact sparse factorization over real algebraic coefficients."""
from poly import sparse, dense, algebraic, exact
from std.testing import assert_equal, assert_true


def main() raises:
    var root2 = dense.root_between([-2, 0, 1], 1, 2)
    # x^2-2*y^2 splits over AA, although it is irreducible over QQ.
    var p = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 2), -2)))
    var result = sparse.factor_checked(p)
    assert_true(exact.check(result))
    var factors = exact.factors(result)
    assert_equal(len(factors), 2)
    # Storage-order monicity makes the x coefficient one.
    var negative = sparse.asarray((((1, 0), 1), ((0, 1), algebraic.negative(root2))))
    var positive = sparse.asarray((((1, 0), 1), ((0, 1), root2.copy())))
    assert_equal(factors[0][0].to_list(), negative.to_list())
    assert_equal(factors[1][0].to_list(), positive.to_list())
    assert_equal(factors[0][1], 1)
    assert_equal(factors[1][1], 1)
    assert_true(exact.unit(result) == algebraic.asreal(1))
    var plain = sparse.factor(p)
    assert_equal(plain[1][0][0].to_list(), factors[0][0].to_list())
    # This real quadratic remains irreducible even though its univariate
    # Kronecker image has several real irreducible factors.
    var sum_squares = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 2), 1)))
    var irreducible = sparse.factor_checked(sum_squares)
    assert_true(exact.check(irreducible))
    assert_equal(len(exact.factors(irreducible)), 1)
    assert_equal(exact.factors(irreducible)[0][0].to_list(), sum_squares.to_list())
    # Repeated factors and an irrational scalar keep their multiplicities.
    var repeated = sparse.multiply(sparse.power(negative, 3), sparse.power(positive, 2))
    repeated = sparse.multiply(repeated, sparse.asarray((((0, 0), algebraic.negative(root2)),)))
    var rr = sparse.factor_checked(repeated)
    assert_true(exact.check(rr))
    var rf = exact.factors(rr)
    assert_equal(len(rf), 2)
    assert_equal(rf[0][0].to_list(), negative.to_list())
    assert_equal(rf[0][1], 3)
    assert_equal(rf[1][0].to_list(), positive.to_list())
    assert_equal(rf[1][1], 2)
    assert_true(exact.unit(rr) == algebraic.negative(root2))
    # A multivariate result may contain an irreducible univariate quadratic.
    var quadratic = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 0), 1)))
    var shifted_y = sparse.asarray((((0, 1), algebraic.asreal(1)), ((0, 0), 1)))
    var grouped = sparse.factor_checked(sparse.multiply(quadratic, shifted_y))
    assert_true(exact.check(grouped))
    assert_equal(exact.factors(grouped)[0][0].to_list(), shifted_y.to_list())
    assert_equal(exact.factors(grouped)[1][0].to_list(), quadratic.to_list())
    # Coefficient-content extraction works with algebraic coefficients too.
    var shifted_x = sparse.asarray((((1, 0), 1), ((0, 0), root2.copy())))
    var algebraic_input = sparse.factor_checked(sparse.multiply(shifted_x, shifted_y))
    assert_true(exact.check(algebraic_input))
    assert_equal(exact.factors(algebraic_input)[0][0].to_list(), shifted_y.to_list())
    assert_equal(exact.factors(algebraic_input)[1][0].to_list(), shifted_x.to_list())
    # Inactive axes and more than two active variables survive reconstruction.
    var a = sparse.asarray((((1, 0, 0, 0), algebraic.asreal(1)), ((0, 0, 1, 0), 1)))
    var b = sparse.asarray((((0, 0, 0, 1), algebraic.asreal(1)), ((0, 0, 0, 0), 1)))
    var many = sparse.factor_checked(sparse.multiply(a, b))
    assert_true(exact.check(many))
    var mf = exact.factors(many)
    assert_equal(len(mf), 2)
    assert_equal(mf[0][0].to_list(), b.to_list())
    assert_equal(mf[1][0].to_list(), a.to_list())
    # Duplicate input terms cancel before the support and degree box are built.
    var cancelled = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 2), -1), ((0, 3), 7), ((0, 3), -7)))
    assert_equal(len(sparse.factor(cancelled)[1]), 2)
    print("array multivariate real algebraic factor checks passed")
