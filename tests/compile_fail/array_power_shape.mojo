# expect: failed to infer
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.sparse.add([((1, 0), 1)], [(1, 1)])
