"""Expand (1 + x)^3 and evaluate it and its derivative exactly."""
import poly
from apn_mojo import Integer
from std.testing import assert_equal


def main() raises:
    var polynomial = poly.dense.power([1, 1], 3)
    var value = poly.dense.evaluate(polynomial, 3)
    var slope = poly.dense.evaluate(poly.dense.derivative(polynomial), 3)
    assert_equal(value, Integer(64))
    assert_equal(slope, Integer(48))
    print("Value at 3:", value)
    print("Derivative at 3:", slope)
