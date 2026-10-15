# expect: modular coefficients cannot mix implicitly with rational or complex coefficients
from poly import sparse
from apn_mojo import Rational

def main() raises:
    var p = sparse.asarray([(1, 1)], modulus=5)
    _ = sparse.groebner((p.copy(), [(0, Rational(1, 2))]))
