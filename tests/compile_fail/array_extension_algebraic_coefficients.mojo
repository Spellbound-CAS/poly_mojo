# expect: real roots require exact integer or rational coefficients
from poly import dense
from apn_mojo import ExactComplex

def main() raises:
    _ = dense.real_roots([ExactComplex(1, 1), ExactComplex(1, 0)])
