# expect: modular coefficients cannot mix implicitly with rational or complex coefficients
from poly import basis, dense
from apn_mojo import Rational

def main() raises:
    var p = basis.asarray(dense.asarray([1], modulus=5), basis="chebyshev")
    _ = basis.evaluate(p, Rational(1, 2))
