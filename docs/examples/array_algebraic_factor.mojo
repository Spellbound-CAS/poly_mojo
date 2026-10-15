"""Choose the factorization field through ordinary coefficient values."""
from poly import dense, sparse, algebraic, coefficients, exact
from std.testing import assert_equal, assert_true


def main() raises:
    # x^2 - 2 splits once its coefficients belong to the real algebraic field.
    var p = coefficients((-2, 0, algebraic.asreal(1)))
    var real_result = dense.factor_checked(p)
    var real_factors = exact.factors(real_result)
    var sqrt_two = dense.root_between([-2, 0, 1], 1, 2)
    assert_true(real_factors[0][0][0] == algebraic.negative(sqrt_two))
    assert_true(real_factors[1][0][0] == sqrt_two)
    assert_true(exact.check(real_result))
    print("AA: x^2 - 2 = (x - sqrt(2))(x + sqrt(2))")

    # Over AA, x^2 + 1 stays quadratic. Over QQbar it has two linear factors.
    var real_pair = dense.factor(coefficients((1, 0, algebraic.asreal(1))))
    var complex_result = dense.factor_checked(coefficients((1, 0, algebraic.ascomplex(1))))
    assert_equal(len(real_pair[1]), 1)
    assert_equal(len(exact.factors(complex_result)), 2)
    assert_true(exact.check(complex_result))
    print("QQbar: x^2 + 1 = (x - i)(x + i)")

    # The active variable can be any axis; output terms keep the same shape.
    var terms = sparse.asarray((((0, 2), algebraic.asreal(1)), ((0, 0), -2)))
    var sparse_result = sparse.factor_checked(sparse.power(terms, 3))
    assert_true(exact.check(sparse_result))
    for factor, multiplicity in exact.factors(sparse_result):
        assert_equal(multiplicity, 3)
        assert_equal(factor.domain(), "AA")
    print("Sparse factors preserve axis 1 and multiplicity 3")
