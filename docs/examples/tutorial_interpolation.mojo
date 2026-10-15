"""Recover a quadratic from exact values at three distinct points."""
from std.testing import assert_equal, assert_true
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, MultipointPlan, check_interpolation


def main() raises:
    var work = Budget()
    var nodes: List[QQ] = [QQ(Rational(0)), QQ(Rational(1)), QQ(Rational(2))]
    var values: List[QQ] = [QQ(Rational(1)), QQ(Rational(4)), QQ(Rational(9))]
    var plan = MultipointPlan[QQ].prepare_interpolation(0, NoContext(), nodes, work)
    var polynomial = plan.interpolate(values, work)
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var expected = DensePoly[QQ].from_sparse(ring.parse("(x + 1)^2", work), work, variable=Int64(0))
    assert_equal(polynomial, expected)
    assert_true(check_interpolation(polynomial, nodes, values, work))
    var recovered = plan.evaluate(polynomial, work)
    for i in range(len(values)):
        assert_equal(recovered[i], values[i])
    print("Recovered (x + 1)^2:", polynomial == expected)
    print("Value at 2:", recovered[2].value)
