"""Factor repeated finite-field and Gaussian-rational polynomials."""
import poly
from apn_mojo import ExactComplex
from std.testing import assert_equal, assert_true


def main() raises:
    # In GF(3), x^3 + 1 = (x + 1)^3, even though its derivative is zero.
    var finite = poly.dense.asarray([1, 0, 0, 1], modulus=3)
    var factors = poly.dense.factor(finite)
    assert_equal(factors[1][0][1], 3)
    print("GF(3) factor multiplicity:", factors[1][0][1])

    # ExactComplex coefficients select QQ(i), so x^2 + 1 splits.
    var complex_p = poly.coefficients((1, 0, ExactComplex(1)))
    var checked = poly.dense.factor_checked(complex_p)
    assert_true(poly.exact.check(checked))
    print("Gaussian-rational factors:", len(poly.exact.factors(checked)))

    var a = poly.sparse.asarray([((1, 0), 1), ((0, 1), 1)])
    var b = poly.sparse.asarray([((0, 1), 1), ((0, 0), 1)])
    var product = poly.sparse.multiply(poly.sparse.power(a, 2), poly.sparse.power(b, 3))
    var blocks = poly.sparse.squarefree_checked(product)
    assert_true(poly.exact.check(blocks))
    var values = poly.exact.factors(blocks)
    print("Sparse multiplicities:", values[0][1], values[1][1])
