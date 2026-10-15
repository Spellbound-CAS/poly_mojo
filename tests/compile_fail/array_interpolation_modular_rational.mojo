# expect: modular coefficients cannot mix
import poly
from apn_mojo import Rational, ExactComplex


def main() raises:
    _ = poly.dense.interpolate([poly.modint(0, 5)], [Rational(1, 2)])
