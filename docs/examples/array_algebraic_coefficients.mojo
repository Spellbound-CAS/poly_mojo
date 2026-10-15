"""Build and calculate with exact real and complex algebraic coefficients."""
from poly import dense, sparse, algebraic, coefficients, exact
from std.testing import assert_true, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var left = dense.asarray(coefficients((algebraic.negative(a), 1)))
    var right = dense.asarray(coefficients((a.copy(), 1)))
    var product = dense.multiply(left, right)
    assert_true(product[0] == algebraic.asreal(-2))
    assert_true(product[1].is_zero())
    assert_true(product[2] == algebraic.asreal(1))
    print(product.domain(), ": (x - sqrt(2))(x + sqrt(2)) = x^2 - 2", sep="")

    var division = dense.divide_checked(product, left)
    assert_true(exact.check(division))
    assert_equal(exact.quotient(division).to_list(), right.to_list())

    var z = algebraic.ascomplex(a, 1)
    var complex_polynomial = dense.asarray(coefficients((algebraic.negative(z), 1)))
    assert_true(dense.evaluate(complex_polynomial, z).is_zero())
    print(complex_polynomial.domain(), ": x - (sqrt(2) + i) vanishes exactly at sqrt(2) + i", sep="")

    var multivariate = sparse.asarray((((1, 0), a.copy()), ((0, 1), 1)))
    var specialized = sparse.specialize(multivariate, [(0, a.copy())])
    assert_true(sparse.evaluate(specialized, (0, 3)) == algebraic.asreal(5))
    print("Sparse specialization: sqrt(2)*x + y becomes 2 + y")
