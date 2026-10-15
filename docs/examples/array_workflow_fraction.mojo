"""Request quotient accuracy without choosing the final working precision."""
from poly import rational
from polynomials import Budget, NumericOptions, BallDomain, numeric_bounds
from apn_mojo import Rational
from std.testing import assert_true

def main() raises:
    var f = rational.asarray([1], [3, 1])
    var result = rational.adaptive_evaluate(
        f, 0, options=NumericOptions(initial_precision=8, max_precision=64,
                                    absolute_bits=30))
    assert_true(result.has_value() and result.accurate())
    assert_true(rational.check(result))
    var work = Budget()
    var bounds = numeric_bounds[BallDomain](result.value(), work)
    assert_true(bounds.real_low.value <= Rational(1, 3))
    assert_true(bounds.real_high.value >= Rational(1, 3))
    print("1/3 enclosed with certified accuracy")
