"""Exact roots, arithmetic, comparisons, and explicit approximation."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind
from poly import dense, sparse, roots, algebraic
from std.testing import assert_equal, assert_true


def main() raises:
    var p = [-2, 0, 1]
    var result = dense.real_roots(p)
    assert_equal(len(result), 2)
    assert_true(roots.check(result))
    assert_equal(algebraic.compare(result[0], 0), -1)
    assert_equal(algebraic.compare(result[1], 0), 1)
    var square = algebraic.multiply(result[1], result[1])
    assert_true(algebraic.check(square))
    assert_equal(algebraic.compare(square, 2), 0)
    var sum_value = algebraic.add(result[0], result[1])
    assert_true(algebraic.check(sum_value))
    assert_equal(algebraic.compare(sum_value, 0), 0)
    assert_equal(algebraic.compare(algebraic.subtract(result[1], result[1]), 0), 0)
    assert_equal(algebraic.compare(algebraic.negative(result[0]), result[1]), 0)
    var quotient = algebraic.divide(6, 3)
    assert_equal(algebraic.compare(quotient, 2), 0)
    assert_true(algebraic.check(quotient))
    assert_true(Bool(quotient.division_evidence()))
    var inverse = algebraic.reciprocal(result[1])
    assert_true(algebraic.check(inverse))
    assert_equal(algebraic.compare(algebraic.multiply(inverse, result[1]), 1), 0)
    var selected = dense.root_between(p, 1, 2)
    assert_equal(algebraic.compare(selected, result[1]), 0)
    var counted = dense.count_roots(p, -2, 2)
    assert_equal(counted.count(), 2)
    assert_true(roots.check(counted))
    var narrow = algebraic.refine(selected, Rational(1, 1024))
    var interval = algebraic.bounds(narrow)
    assert_true(interval[1] - interval[0] < Rational(1, 1024))
    var approximation = algebraic.approximate(selected, Rational(1, 1024))
    assert_true(approximation > Rational(1) and approximation < Rational(2))
    var count = 0
    for root in result:
        assert_true(algebraic.check(root))
        count += 1
    assert_equal(count, 2)
    var repeated = dense.real_roots([1, -2, 1])
    assert_equal(len(repeated), 1)
    assert_equal(algebraic.compare(repeated[0], 1), 0)
    assert_true(roots.check(repeated))
    assert_equal(len(dense.real_roots([7])), 0)
    var rational = dense.real_roots([-1, 3])
    assert_equal(algebraic.bounds(rational[0])[0], Rational(1, 3))
    var sp = [(0, -2), (2, 1)]
    assert_equal(sparse.count_roots(sp, -2, 2).count(), 2)
    assert_equal(len(sparse.real_roots(sp)), 2)
    assert_equal(algebraic.compare(sparse.root_between(sp, 1, 2), selected), 0)
    var mixed_roots = dense.real_roots([2, -6, -1, 3])
    assert_true(roots.check(mixed_roots))
    assert_equal(len(mixed_roots), 3)
    assert_equal(algebraic.compare(mixed_roots[0], result[0]), 0)
    assert_equal(algebraic.compare(mixed_roots[1], Rational(1, 3)), 0)
    assert_equal(algebraic.compare(mixed_roots[2], result[1]), 0)
    print("array root value checks passed")
