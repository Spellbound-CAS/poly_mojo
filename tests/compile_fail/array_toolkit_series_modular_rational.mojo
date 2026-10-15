# expect: modular coefficients cannot mix implicitly with rational or complex coefficients
from poly import series, dense
from apn_mojo import Rational

def main() raises:
    var finite = series.asarray(dense.asarray([1], modulus=5), order=3)
    _ = series.add(finite, [Rational(1, 2)])
