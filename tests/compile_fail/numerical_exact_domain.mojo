# expect: CoefficientRing
"""Rounded scalars cannot satisfy the exact polynomial ring contract."""
from apn_mojo import Float
from polynomials import DensePoly


def reject(value: DensePoly[Float]):
    _ = value


def main():
    pass
