# expect: convert floats explicitly
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.asarray([(0, 1.5)])
