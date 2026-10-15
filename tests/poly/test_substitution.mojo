"""Ordered polynomial reduction and simultaneous multivariate changes."""
import poly
from apn_mojo import Rational
from polynomials import PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var dense_divisors: List[List[Int]] = [[1, 2], [1]]
    var reduction = poly.dense.divide_many([1, 0, 1], dense_divisors)
    var restored = poly.dense.add(poly.dense.multiply(reduction[0][0], [1, 2]), poly.dense.multiply(reduction[0][1], [1]))
    assert_equal(poly.dense.add(restored, reduction[1]).to_list(), poly.coefficients((Rational(1), 0, 1)).to_list())
    var empty = poly.dense.divide_many([1, 2], List[List[Int]]())
    assert_equal(len(empty[0]), 0)
    assert_equal(empty[1].to_list(), poly.coefficients((Rational(1), 2)).to_list())
    var x = poly.sparse.asarray([((1, 0), 1)])
    var y = poly.sparse.asarray([((0, 1), 1)])
    var p = poly.sparse.asarray([((2, 1), 1), ((0, 1), 1)])
    var swapped = poly.sparse.substitute(p, [y.copy(), x.copy()])
    assert_equal(swapped.to_list(), poly.sparse.asarray([((1, 2), 1), ((1, 0), 1)]).to_list())
    var replaced = poly.sparse.substitute(p, [y.copy()], axes=[0])
    assert_equal(replaced.to_list(), poly.sparse.asarray([((0, 3), 1), ((0, 1), 1)]).to_list())
    var divisors = [poly.sparse.asarray([((1, 1), 1), ((0, 0), 1)]), y.copy()]
    var reduced = poly.sparse.divide_many(p, divisors)
    var restored_sparse = reduced[1].copy()
    for i in range(len(divisors)):
        restored_sparse = poly.sparse.add(restored_sparse, poly.sparse.multiply(reduced[0][i], divisors[i]))
    assert_equal(poly.sparse.evaluate(restored_sparse, (2, 3)), Rational(15))
    var transformed = poly.sparse.affine_transform(p, 2, 1, axis=0)
    assert_equal(poly.sparse.evaluate(transformed, (1, 2)), poly.sparse.evaluate(p, (3, 2)))
    var shifted = poly.sparse.taylor_shift(p, 1, axis=1)
    assert_equal(poly.sparse.evaluate(shifted, (2, 2)), poly.sparse.evaluate(p, (2, 3)))
    var integral = poly.sparse.definite_integral(p, 0, 1, axis=0)
    assert_equal(integral.to_list(), poly.sparse.asarray([((0, 1), Rational(4, 3))]).to_list())
    var failed = False
    try:
        _ = poly.sparse.substitute(p, [x.copy(), y.copy()], axes=[0, 0])
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    failed = False
    try:
        _ = poly.sparse.divide_many(p, divisors, max_steps=0)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    print("array substitution checks passed")
