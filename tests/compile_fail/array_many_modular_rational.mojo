# expect: modular coefficients cannot mix
import poly
from apn_mojo import Rational


def main() raises:
    var p = poly.dense.asarray([1, 2], modulus=7)
    _ = poly.dense.evaluate_many(p, [Rational(1, 2)])
