# expect: sequence entries must be exact numbers
import poly
from apn_mojo import Rational, ExactComplex


def main() raises:
    _ = poly.dense.interpolate([0.0, 1.0], [1, 2])
