# expect: sparse tuples need
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.asarray(((0, 1), (2, Float64(1.5))))
