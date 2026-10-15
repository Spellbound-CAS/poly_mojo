"""Reject malformed receipts and accept independently refined input witnesses."""
from poly import dense
from apn_mojo import BallContext
from polynomials import Budget
from std.testing import assert_true, assert_false


def main() raises:
    var budget = Budget()
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = [a.copy()]
    var bc = BallContext(precision=48)
    var stream = dense.enclose_stream(p, 0, context=bc)
    var checker = dense.enclose_checker(p, 0, context=bc)
    var truncated = dense.enclose_checker(p, 0, context=bc)
    var wrong_context = dense.enclose_checker(p, 0, context=BallContext(precision=80))
    var alternative = dense.enclose_stream(p, 0, context=BallContext(precision=96))
    # Supply a valid tighter input witness, then the 48-bit arithmetic stream.
    # A checker that reran bisection and demanded identical witnesses would fail.
    var tighter = alternative.next(budget)
    assert_true(checker.accept(tighter.take(), budget))
    var first = stream.next(budget)
    assert_true(truncated.accept(first.value().copy(), budget))
    assert_true(wrong_context.accept(first.take(), budget))
    while True:
        var step = stream.next(budget)
        if not step:
            break
        assert_true(wrong_context.accept(step.value().copy(), budget))
        assert_true(checker.accept(step.take(), budget))
    var receipt = stream.finish(budget)
    assert_true(checker.finish(receipt, budget))
    assert_false(truncated.finish(receipt, budget))
    assert_false(wrong_context.finish(receipt, budget))
    var malformed = receipt.copy()
    malformed._value.clear()
    assert_false(checker.finish(malformed, budget))
    malformed = receipt.copy()
    malformed._value.append(receipt._value[0].copy())
    assert_false(checker.finish(malformed, budget))
    assert_true(checker.finish(receipt, budget))
    print("array algebraic refinement receipt checks passed")
