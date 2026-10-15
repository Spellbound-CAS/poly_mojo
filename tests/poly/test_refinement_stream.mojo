"""Algebraic coefficient and coordinate refinement in bounded-history streams."""
from poly import dense, sparse, algebraic, evaluation, coefficients
from apn_mojo import BallContext, Rational, Integer, contains
from polynomials import Budget, NumericDomain, BallDomain, ComplexBallDomain, numeric_bounds, QQ, PolyErrorKind
from std.testing import assert_true, assert_false, assert_equal


def consume[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](var stream: evaluation.AlgebraicEnclosureStream[C, D, numerical, sparse, False], var checker: evaluation.AlgebraicEnclosureStream[C, D, numerical, sparse, True], mut budget: Budget) raises -> D.Scalar:
    assert_equal(stream.step_count(), checker.step_count())
    var rejected = False
    try:
        _ = stream.finish(budget)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    var count = 0
    var last_kind = 0
    while True:
        var step = stream.next(budget)
        if not step:
            break
        assert_true(step.value().kind >= last_kind)
        last_kind = step.value().kind
        assert_true(checker.accept(step.take(), budget))
        count += 1
    assert_false(Bool(stream.next(budget)))
    var receipt = stream.finish(budget)
    assert_equal(count, stream.step_count())
    assert_equal(count, receipt.steps)
    assert_true(checker.finish(receipt, budget))
    return receipt.value(budget)


def narrow[D: NumericDomain](value: D.Scalar, bits: Int, mut budget: Budget) raises:
    var box = numeric_bounds[D](value, budget)
    var target = QQ(Rational(Integer(1), Integer(1) << bits))
    assert_true(box.real_high.sub(box.real_low, budget).sub(target, budget).value.sign() <= 0)
    assert_true(box.imag_high.sub(box.imag_low, budget).sub(target, budget).value.sign() <= 0)


def main() raises:
    var budget = Budget()
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var bc = BallContext(precision=80)
    var p = dense.asarray(coefficients((a.copy(), 1)))
    var value = consume(dense.enclose_stream(p, 1, context=bc), dense.enclose_checker(p, 1, context=bc), budget)
    narrow[BallDomain](value, 70, budget)
    assert_false(contains(value, Rational(5, 2)))
    var low = consume(dense.enclose_stream([a.copy()], 0, context=BallContext(precision=12)), dense.enclose_checker([a.copy()], 0, context=BallContext(precision=12)), budget)
    var high = consume(dense.enclose_stream([a.copy()], 0, context=bc), dense.enclose_checker([a.copy()], 0, context=bc), budget)
    var lb = numeric_bounds[BallDomain](low, budget)
    var hb = numeric_bounds[BallDomain](high, budget)
    assert_true(hb.real_high.sub(hb.real_low, budget).sub(lb.real_high.sub(lb.real_low, budget), budget).value.sign() < 0)
    # Independent rational bounds locate sqrt(2), rather than comparing two producers.
    assert_true(hb.real_low.mul(hb.real_low, budget).sub(QQ(Rational(2)), budget).value.sign() < 0)
    assert_true(hb.real_high.mul(hb.real_high, budget).sub(QQ(Rational(2)), budget).value.sign() > 0)
    var product = consume(dense.enclose_stream(coefficients((0, a.copy())), a, context=bc), dense.enclose_checker(coefficients((0, a.copy())), a, context=bc), budget)
    assert_true(contains(product, Rational(2)))
    narrow[BallDomain](product, 68, budget)
    var z = algebraic.ascomplex(a, a)
    var complex = consume(dense.enclose_stream([z.copy()], 1, context=bc), dense.enclose_checker([z.copy()], 1, context=bc), budget)
    narrow[ComplexBallDomain](complex, 70, budget)
    var point = algebraic.ascomplex(a, 1)
    var at_point = consume(dense.enclose_stream([0, 0, 1], point, context=bc), dense.enclose_checker([0, 0, 1], point, context=bc), budget)
    assert_true(contains(at_point.real(), Rational(1)))
    narrow[ComplexBallDomain](at_point, 68, budget)
    var sp = sparse.asarray([((1, 0), z.copy())])
    var sv = consume(sparse.enclose_stream(sp, (a.copy(), 1), context=bc), sparse.enclose_checker(sp, (a.copy(), 1), context=bc), budget)
    assert_true(contains(sv.real(), Rational(2)))
    assert_true(contains(sv.imag(), Rational(2)))
    narrow[ComplexBallDomain](sv, 68, budget)
    var unit = consume(sparse.enclose_stream([(1, a.copy())], a, context=bc), sparse.enclose_checker([(1, a.copy())], a, context=bc), budget)
    assert_true(contains(unit, Rational(2)))
    var constant = consume(sparse.enclose_stream([((), a.copy())], (), context=bc), sparse.enclose_checker([((), a.copy())], (), context=bc), budget)
    narrow[BallDomain](constant, 70, budget)
    var numeric = dense.asarray([0, 0, 1], context=bc)
    var nv = consume(dense.enclose_stream(numeric, a, context=bc), dense.enclose_checker(numeric, a, context=bc), budget)
    assert_true(contains(nv, Rational(2)))
    var ns = sparse.asarray([((2, 0), 1)], context=bc)
    var nsv = consume(sparse.enclose_stream(ns, (a.copy(), 0), context=bc), sparse.enclose_checker(ns, (a.copy(), 0), context=bc), budget)
    assert_true(contains(nsv, Rational(2)))
    var un = sparse.asarray([(2, 1)], context=bc)
    var uns = consume(sparse.enclose_stream(un, a, context=bc), sparse.enclose_checker(un, a, context=bc), budget)
    assert_true(contains(uns, Rational(2)))
    var empty = List[algebraic.Real]()
    var ev = consume(dense.enclose_stream(empty, 0, context=bc), dense.enclose_checker(empty, 0, context=bc), budget)
    assert_true(contains(ev, Rational(0)))
    var zero = consume(dense.enclose_stream(empty, a, context=bc), dense.enclose_checker(empty, a, context=bc), budget)
    assert_true(contains(zero, Rational(0)))
    var exact = consume(dense.enclose_stream([algebraic.multiply(a, a)], 0, context=bc), dense.enclose_checker([algebraic.multiply(a, a)], 0, context=bc), budget)
    assert_true(contains(exact, Rational(2)))
    assert_equal(a.witness().low, QQ(Rational(1)))
    assert_equal(a.witness().high, QQ(Rational(2)))
    assert_true(algebraic.check(a))
    print("array algebraic refinement stream checks passed")
