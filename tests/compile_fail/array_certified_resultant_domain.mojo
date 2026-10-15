# expect: modular coefficients cannot mix implicitly with rational or complex coefficients
from poly import dense, modint
from apn_mojo import Rational

def main() raises:
    _ = dense.resultant_certified([modint(1, 5)], [Rational(1, 2)])
