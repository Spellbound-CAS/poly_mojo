"""Scalar arithmetic keeps coefficient input concise and sparse powers sparse."""
from poly import dense, sparse
from apn_mojo import ArithmeticContext, FloatFormat, Rational
from std.testing import assert_equal


def main() raises:
    var p = dense.asarray([1, 2, 3])
    var shifted = dense.subtract(4, p)
    var half = dense.divide(p, 2)
    assert_equal(dense.evaluate(shifted, 2), -13)
    assert_equal(half[0], Rational(1, 2))
    print("Integer division promotes to", half.domain())
    var context = ArithmeticContext(format=FloatFormat(precision=80))
    var high_power = sparse.asarray([(1000000000, 6)], context=context)
    var scaled = sparse.divide(high_power, 2, context=context, max_steps=1000000)
    assert_equal(scaled[0][1].to_rational_exact(), Rational(3))
    print("Stored terms after scaling:", len(scaled))
