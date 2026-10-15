# expect: wrong number of coordinates
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.evaluate([((1, 0), 1)], (1, 2, 3))
