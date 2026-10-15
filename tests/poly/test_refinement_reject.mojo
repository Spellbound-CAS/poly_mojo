"""Refinement evidence is bound to its source, position, and complete arithmetic trace."""
from poly import dense, sparse, algebraic, evaluation
from apn_mojo import BallContext, Rational, Ball
from polynomials import Budget, QQ, NumericBox, BallDomain, numeric_value
from std.testing import assert_true, assert_false


def main() raises:
    var budget = Budget()
    var bc = BallContext(precision=48)
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = [a.copy(), algebraic.asreal(1)]
    var stream = dense.enclose_stream(p, 1, context=bc)
    var pending_1 = stream.next(budget)
    var first = pending_1.take()
    var bad = first.copy()
    bad.refinement.value().real = dense.root_between([-2, 0, 1], -2, -1).witness().copy()
    var checker = dense.enclose_checker(p, 1, context=bc)
    assert_false(checker.accept(bad^, budget))
    assert_false(checker.accept(first.copy(), budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.refinement.value().real.low = QQ(Rational(3))
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.refinement.value().real.sequence.polynomials.clear()
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.refinement.value().real.polynomial.coefficients.clear()
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.refinement.value().value = numeric_value[BallDomain](NumericBox.point(QQ(Rational(9)), QQ(Rational(0)), budget), bc, budget)
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.refinement.value().imag = algebraic.asreal(1).witness().copy()
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.index += 1
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.kind = 1
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.refinement = None
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    bad = first.copy()
    bad.bounds = NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), budget)
    assert_false(checker.accept(bad^, budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    assert_true(checker.accept(first.copy(), budget))
    assert_false(checker.accept(first.copy(), budget))
    checker = dense.enclose_checker(p, 1, context=bc)
    var pending_2 = stream.next(budget)
    var second = pending_2.take()
    assert_false(checker.accept(second.copy(), budget))
    var different = dense.enclose_checker([algebraic.asreal(3), algebraic.asreal(1)], 1, context=bc)
    assert_false(different.accept(first.copy(), budget))
    # A separately validated refinement may be tighter than the producer's;
    # the checker verifies evidence, not a replay of its bisection schedule.
    var independent = dense.enclose_stream(p, 1, context=BallContext(precision=80))
    var pending_3 = independent.next(budget)
    var tighter = pending_3.take()
    checker = dense.enclose_checker(p, 1, context=bc)
    assert_true(checker.accept(tighter^, budget))
    assert_true(checker.accept(second.copy(), budget))
    var pending_4 = stream.next(budget)
    var arithmetic = pending_4.take()
    bad = arithmetic.copy()
    bad.index = 999
    assert_false(checker.accept(bad^, budget))
    # Coordinate evidence has its own source binding, including mixed tuples.
    var points = sparse.enclose_stream([((1, 1), 1)], (2, a.copy()), context=bc)
    var pc = sparse.enclose_checker([((1, 1), 1)], (2, a.copy()), context=bc)
    var pending_5 = points.next(budget)
    var point = pending_5.take()
    point.index = 0
    assert_false(pc.accept(point^, budget))
    points = sparse.enclose_stream([((1, 1), 1)], (2, a.copy()), context=bc)
    pc = sparse.enclose_checker([((1, 1), 1)], (2, a.copy()), context=bc)
    var pending_6 = points.next(budget)
    point = pending_6.take()
    point.refinement.value().real = algebraic.asreal(2).witness().copy()
    assert_false(pc.accept(point^, budget))
    var z = algebraic.ascomplex(a, a)
    var cs = dense.enclose_stream([z.copy()], 0, context=bc)
    var cc = dense.enclose_checker([z.copy()], 0, context=bc)
    var pending_7 = cs.next(budget)
    var complex = pending_7.take()
    complex.refinement.value().imag = dense.root_between([-2, 0, 1], -2, -1).witness().copy()
    assert_false(cc.accept(complex^, budget))
    # Completion requires all events, the right precision, and the actual result.
    stream = dense.enclose_stream(p, 1, context=bc)
    checker = dense.enclose_checker(p, 1, context=bc)
    var truncated = dense.enclose_checker(p, 1, context=bc)
    var last = Optional[evaluation.EnclosureStep[BallDomain]]()
    while True:
        var step = stream.next(budget)
        if not step:
            break
        last = step.value().copy()
        assert_true(checker.accept(step.take(), budget))
    var receipt = stream.finish(budget)
    assert_true(checker.finish(receipt, budget))
    assert_false(truncated.finish(receipt, budget))
    receipt.steps -= 1
    assert_false(checker.finish(receipt, budget))
    receipt.steps += 1
    receipt.precision += 1
    assert_false(checker.finish(receipt, budget))
    receipt.precision -= 1
    var original = receipt._value[0].copy()
    receipt._value[0] = numeric_value[BallDomain](NumericBox.point(QQ(Rational(99)), QQ(Rational(0)), budget), bc, budget)
    assert_false(checker.finish(receipt, budget))
    receipt._value[0] = original^
    assert_false(checker.accept(last.take(), budget))
    assert_false(checker.finish(receipt, budget))
    print("array algebraic refinement rejection checks passed")
