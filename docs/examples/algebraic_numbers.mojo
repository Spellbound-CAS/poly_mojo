"""Construct and independently check exact algebraic arithmetic."""
from std.testing import assert_true
from apn_mojo import Rational
from polynomials import Budget, DensePoly, QQ, NoContext, real_root, algebraic_multiply, check_algebraic_arithmetic


def main() raises:
    var budget = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], budget)
    var sqrt2 = real_root(f, QQ(Rational(1)), QQ(Rational(2)), budget)
    var square = algebraic_multiply(sqrt2, sqrt2, budget)
    assert_true(check_algebraic_arithmetic(sqrt2, sqrt2.copy(), square, budget))
    assert_true(square.value.rational and square.value.low.value == Rational(2))
    print("sqrt(2) squared: 2")
    print("Arithmetic certificate checked")
