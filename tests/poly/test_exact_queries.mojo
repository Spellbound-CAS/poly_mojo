"""Closed-interval signs and full determinant certificates with independent answers."""
from std.testing import assert_equal, assert_true, assert_false
from apn_mojo import Rational
from polynomials import Budget, PolyError, PolyErrorKind, QQ
from poly import dense, sparse, exact


def main() raises:
    assert_equal(dense.sign_on_interval([-2, 0, 1], 2, 3).value(), 1)
    assert_equal(dense.sign_on_interval([-2, 0, 1], 0, Rational(1, 2)).value(), -1)
    assert_false(Bool(dense.sign_on_interval([-2, 0, 1], -2, 2)))
    assert_false(Bool(dense.sign_on_interval([1, -2, 1], 0, 2)))
    assert_false(Bool(dense.sign_on_interval([-1, 1], 1, 2)))
    assert_false(Bool(dense.sign_on_interval([-1, 1], 0, 1)))
    assert_false(Bool(dense.sign_on_interval([0], 0, 0)))
    assert_equal(dense.sign_on_interval([-2, 0, 1], 0, 0).value(), -1)
    assert_equal(dense.sign_on_interval([3], -2, 2).value(), 1)
    var budget = Budget()
    assert_equal(sparse.sign_on_interval([(0, -2), (2, 1)], 2, 3, budget).value(), 1)
    assert_equal(sparse.sign_on_interval([((0,), -3)], 0, 1).value(), -1)
    var linear = dense.resultant_certified([-2, 1], [Rational(-3), Rational(1)])
    assert_true(exact.check(linear))
    assert_equal(exact.value(linear).domain(), "QQ")
    assert_equal(exact.value(linear).to_list(), [Rational(-1)])
    var quadratic = dense.resultant_certified([1, 0, 1], [-2, 0, 1], budget)
    assert_true(exact.check(quadratic, budget))
    assert_equal(exact.value(quadratic, budget).to_list(), [Rational(9)])
    var common = dense.resultant_certified([-1, 0, 1], [-1, 1])
    assert_true(exact.check(common))
    assert_equal(len(exact.value(common)), 0)
    var constant = dense.resultant_certified([2], [1, 0, 1])
    assert_true(exact.check(constant))
    assert_equal(exact.value(constant)[0], Rational(4))
    assert_equal(exact.value(dense.resultant_certified([2], [3]))[0], Rational(1))
    var zero = dense.resultant_certified(List[Int](), [1, 2])
    assert_true(exact.check(zero))
    assert_equal(len(exact.value(zero)), 0)
    # A valid transcript for another axis cannot replace the requested claim.
    zero._data[].proof.variable = 1
    assert_false(exact.check(zero))
    var f = sparse.asarray([((1, 0), 1), ((0, 1), -1)])
    var rx = sparse.resultant_certified(f, [((1, 0), 1), ((0, 0), -2)], budget)
    assert_true(exact.check(rx))
    assert_equal(exact.value(rx).to_list(), sparse.asarray([((0, 1), Rational(1)), ((0, 0), Rational(-2))]).to_list())
    var ry = sparse.resultant_certified(f, [((0, 1), 1), ((0, 0), -2)], axis=1)
    assert_true(exact.check(ry))
    assert_equal(exact.value(ry).to_list(), sparse.asarray([((1, 0), Rational(-1)), ((0, 0), Rational(2))]).to_list())
    var copy = dense.to_core(exact.value(quadratic))
    copy.coefficients[0] = QQ(Rational(123))
    assert_true(exact.check(quadratic))
    quadratic._data[].proof.steps[0].pivot_row = -1
    assert_false(exact.check(quadratic))
    linear._data[].proof.resultant.terms.clear()
    assert_false(exact.check(linear))
    for scenario in range(5):
        var failed = False
        try:
            if scenario == 0:
                _ = dense.sign_on_interval([1], 2, 1)
            elif scenario == 1:
                _ = sparse.sign_on_interval([((1, 0), 1)], 0, 1)
            elif scenario == 2:
                _ = sparse.resultant_certified(f, f, axis=2)
            elif scenario == 3:
                _ = dense.sign_on_interval([1], 0, 1, max_steps=0)
            else:
                _ = dense.resultant_certified([1, 2], [3, 4], max_steps=0)
        except e:
            failed = e.kind == (PolyErrorKind.LIMIT if scenario >= 3 else PolyErrorKind.DOMAIN)
        assert_true(failed)
    print("array exact query checks passed")
