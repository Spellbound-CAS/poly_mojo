# expect: modular coefficients cannot mix
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.add([(0, poly.modint(1, 7))], [(0, ExactComplex(0, 1))])
