# expect: modular coefficients cannot mix
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.asarray(((0, poly.modint(1, 7)), (2, ExactComplex(0, 1))))
