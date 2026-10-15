"""Algebraic interval queries and checked numerical diagnostics."""
from poly import dense, sparse, algebraic, roots, coefficients
from apn_mojo import BallContext, Rational, ComplexContext, ArithmeticContext, FloatFormat
from polynomials import Budget, PolyError, PolyErrorKind, NumericOptions
from std.testing import assert_equal, assert_true, assert_false

def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray(coefficients((algebraic.negative(a), 1)))
    var count = dense.count_roots(p, 1, 2)
    assert_equal(count.count(), 1)
    assert_true(roots.check(count))
    count._count = 2
    assert_false(roots.check(count))
    assert_equal(algebraic.compare(dense.root_between(p, 1, 2), a), 0)
    assert_equal(algebraic.compare(dense.root_between(p, a, a), a), 0)
    assert_equal(dense.sign_on_interval(p, 2, 3).value(), 1)
    assert_false(Bool(dense.sign_on_interval(p, 1, a)))
    assert_equal(dense.sign_on_interval(p, 0, 1).value(), -1)
    assert_true(roots.check(dense.count_roots([-2, 0, 1], 0, algebraic.add(a, 1))))
    var repeated = dense.multiply(p, p)
    assert_equal(dense.count_roots(repeated, 1, 2).count(), 1)
    assert_false(Bool(dense.sign_on_interval(repeated, 1, 2)))
    var sp = sparse.asarray([((1, 0), algebraic.asreal(1)), ((0, 1), algebraic.negative(a))])
    var result = sparse.roots_at(sp, [(1, 1)], axis=0)
    assert_equal(len(result), 1)
    assert_equal(algebraic.compare(result[0], a), 0)
    assert_true(roots.check(result))
    for mode in range(3):
        var domain_rejected = False
        try:
            if mode == 0: _ = dense.count_roots(p, a, 2)
            elif mode == 1: _ = dense.count_roots(p, 2, 1)
            elif mode == 2: _ = dense.root_between(p, 2, 3)
        except e: domain_rejected = e.kind == PolyErrorKind.DOMAIN
        assert_true(domain_rejected)
    var rejected = False
    try:
        _ = dense.count_roots(p, 1, 2, max_steps=0)
    except e: rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    print("completion algebraic root checks passed")
