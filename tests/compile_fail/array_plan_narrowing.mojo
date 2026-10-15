# expect: input exceeds the plan field
import poly
from apn_mojo import Rational, ExactComplex


def main() raises:
    var plan = poly.dense.prepare_points([0, 1])
    _ = plan.evaluate([ExactComplex(0, 1)])
