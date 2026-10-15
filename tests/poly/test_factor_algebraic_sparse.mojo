"""Univariate factorization preserves sparse shapes, axes and coefficient fields."""
from poly import dense, sparse, algebraic, exact
from std.testing import assert_equal, assert_true


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = sparse.asarray(((2, algebraic.asreal(1)), (0, -2)))
    var factorization = sparse.factor(p)
    assert_equal(len(factorization[1]), 2)
    assert_equal(factorization[1][0][0].to_list(), sparse.asarray(((1, 1), (0, algebraic.negative(a)))).to_list())
    var q = sparse.asarray((((0, 0, 2), algebraic.ascomplex(1)), ((0, 0, 0), 1)))
    var result = sparse.factor_checked(q)
    assert_true(exact.check(result))
    var factors = exact.factors(result)
    assert_equal(len(factors), 2)
    assert_equal(factors[0][0].to_list(), sparse.asarray((((0, 0, 1), 1), ((0, 0, 0), algebraic.ascomplex(0, -1)))).to_list())
    assert_equal(factors[1][0].to_list(), sparse.asarray((((0, 0, 1), 1), ((0, 0, 0), algebraic.ascomplex(0, 1)))).to_list())
    assert_equal(sparse.factor(q)[1][0][0].to_list(), factors[0][0].to_list())
    var tuple_one = sparse.asarray((((2,), algebraic.asreal(1)), ((0,), 1)))
    assert_equal(sparse.factor(tuple_one)[1][0][0].to_list(), tuple_one.to_list())
    # Active axes are determined after duplicate terms cancel.
    var cancelled = sparse.asarray((((1, 0), a.copy()), ((1, 0), algebraic.negative(a)), ((0, 2), 1), ((0, 0), -2)))
    var cancellation = sparse.factor_checked(cancelled)
    assert_true(exact.check(cancellation))
    assert_equal(exact.factors(cancellation)[0][0].to_list(), sparse.asarray((((0, 1), 1), ((0, 0), algebraic.negative(a)))).to_list())
    var constant = sparse.asarray((((0, 0), a.copy()),))
    assert_true(exact.check(sparse.factor_checked(constant)))
    assert_equal(len(sparse.factor(constant)[1]), 0)
    var zero = sparse.subtract(constant, constant)
    assert_true(exact.check(sparse.factor_checked(zero)))
    assert_true(sparse.factor(zero)[0].is_zero())
    var many = sparse.asarray((((1, 0), a.copy()), ((0, 1), 1)))
    assert_equal(len(sparse.factor(many)[1]), 1)
    assert_true(exact.check(sparse.factor_checked(many)))
    var specialized = sparse.specialize(many, [(0, algebraic.asreal(1))])
    assert_true(exact.check(sparse.factor_checked(specialized)))
    print("array sparse algebraic factorization checks passed")
