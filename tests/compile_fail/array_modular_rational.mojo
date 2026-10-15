# expect: modular coefficients cannot mix
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.dense.add([poly.modint(1, 7)], [Rational(1, 2)])
