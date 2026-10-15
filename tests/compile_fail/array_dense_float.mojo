# expect: convert floats explicitly
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    _ = poly.dense.asarray([1.0, 2.0])
