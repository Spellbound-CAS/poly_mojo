# expect: each power must be Int
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.asarray([((Rational(1, 2), Int(0)), 1)])
