# expect: CAD requires exact real coefficients
from poly import cad
from apn_mojo import ExactComplex

def main() raises:
    _ = cad.negative([(0, ExactComplex(0, 1))])
