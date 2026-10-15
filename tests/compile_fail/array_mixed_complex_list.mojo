# expect: error:
import poly
from apn_mojo import Rational, ExactComplex

def main() raises:
    var n = Int(1)
    var z = ExactComplex(0, 1)
    var values = [n, z]
    _ = poly.dense.asarray(values)
