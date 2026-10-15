"""Refine and check algebraic inputs without retaining the event history."""
from apn_mojo import BallContext, Rational, contains
from polynomials import Budget
from poly import dense, coefficients
from std.testing import assert_true


def main() raises:
    var budget = Budget()
    var root = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray(coefficients((0, root.copy())))
    var context = BallContext(precision=80)
    var stream = dense.enclose_stream(p, root, budget, context=context)
    var checker = dense.enclose_checker(p, root, budget, context=context)
    while True:
        var step = stream.next(budget)
        if not step:
            break
        assert_true(checker.accept(step.take(), budget))
    var receipt = stream.finish(budget)
    assert_true(checker.finish(receipt, budget))
    assert_true(contains(receipt.value(budget), Rational(2)))
    print("Algebraic stream verified")
