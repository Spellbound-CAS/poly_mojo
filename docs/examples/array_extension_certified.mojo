"""Check evaluation bounds separately from requested accuracy."""
from apn_mojo import Rational, Ball, BallContext, contains
from polynomials import NumericOptions, NumericOutcome
from poly import dense, evaluation
from std.testing import assert_equal, assert_true


def main() raises:
    var fixed = dense.enclose([1, 2, 3], Rational(1, 3), context=BallContext(precision=80))
    assert_true(evaluation.check(fixed))
    assert_true(contains(fixed.value(), Rational(2)))
    var options = NumericOptions(initial_precision=16, max_precision=128, absolute_bits=64)
    var refined = dense.adaptive_evaluate([1, 2, 3], Rational(1, 3), options=options)
    assert_true(evaluation.check(refined))
    assert_true(refined.accurate())
    assert_true(contains(refined.value(), Rational(2)))
    var uncertain = dense.adaptive_evaluate([0, 1], Ball(Rational(0), Rational(1), precision=80), options=options)
    assert_true(evaluation.check(uncertain))
    assert_true(not uncertain.accurate())
    assert_equal(uncertain.outcome().value(), NumericOutcome.PRECISION_CEILING)
    print("Evaluation enclosure checked; requested accuracy reached")
    print("Uncertain input: valid enclosure, precision ceiling reached")
