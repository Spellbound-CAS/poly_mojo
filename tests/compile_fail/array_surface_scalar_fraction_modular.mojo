# expect: modular coefficients cannot mix implicitly with rational or complex coefficients
from poly import sparse, modint
from apn_mojo import Rational

def main() raises:
    _ = sparse.multiply([(1, modint(2, 7))], Rational(1, 2))
