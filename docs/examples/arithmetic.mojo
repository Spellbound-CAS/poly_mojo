"""Compute exact powers, derivatives, and kernel bytes with one shared budget."""
from std.testing import assert_equal
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, MonomialOrder, SparsePoly, DensePoly, evaluate, evaluate_dense, derivative, to_kernel_bytes, from_kernel_bytes


def main() raises:
    var budget = Budget(1000000)
    var order = MonomialOrder()
    var x = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var one = SparsePoly[QQ].constant(QQ(Rational(1)), order, budget)
    var p = x.add(one, budget).pow(3, budget)
    var point: List[Tuple[Int64, QQ]] = [(Int64(0), QQ(Rational(3)))]
    var value = evaluate(p, point, budget)
    var slope = evaluate(derivative(p, 0, budget), point, budget)
    assert_equal(value, QQ(Rational(64)))
    assert_equal(slope, QQ(Rational(48)))

    var dense = DensePoly[QQ].from_sparse(p, budget)
    assert_equal(evaluate_dense(dense, point[0][1], budget), value)
    var bytes = to_kernel_bytes(p, budget)
    assert_equal(from_kernel_bytes(bytes[:], budget), p)
    assert_equal(len(bytes), 100)
    print("Value at 3:", value.value)
    print("Derivative at 3:", slope.value)
    print("Kernel bytes:", len(bytes))
