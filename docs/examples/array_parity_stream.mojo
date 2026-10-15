"""Check an enclosure without keeping the full witness trace."""
from apn_mojo import BallContext, Rational, contains
from polynomials import Budget
from poly import dense
from std.testing import assert_true


def main() raises:
    var budget = Budget()
    var context = BallContext(precision=80)
    var stream = dense.enclose_stream([1, 2, 3], Rational(1, 3), budget, context=context)
    var checker = dense.enclose_checker([1, 2, 3], Rational(1, 3), budget, context=context)
    while True:
        var step = stream.next(budget)
        if not step:
            break
        assert_true(checker.accept(step.take(), budget))
    var receipt = stream.finish(budget)
    assert_true(checker.finish(receipt, budget))
    assert_true(contains(receipt.value(budget), Rational(2)))
    print("Streamed enclosure verified")
