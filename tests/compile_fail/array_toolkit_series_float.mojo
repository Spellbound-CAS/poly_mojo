# expect: exact coefficients must be integers, Rational, ExactComplex or ModInt
from poly import series

def main() raises:
    _ = series.asarray([Float64(1)], order=2)
