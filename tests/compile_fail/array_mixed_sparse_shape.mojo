# expect: same power shape
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.asarray(((0, 1), ((1, 2), Rational(1, 2))))
