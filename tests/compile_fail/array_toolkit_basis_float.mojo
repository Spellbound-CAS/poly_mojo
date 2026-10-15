# expect: basis coefficients must be exact numbers
from poly import basis

def main() raises:
    _ = basis.asarray([Float64(1)], basis="chebyshev")
