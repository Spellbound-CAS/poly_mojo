# expect: expression must be mutable in assignment
from poly import dense
from polynomials import QQ
from apn_mojo import Rational

def main() raises:
    var result = dense.divide_checked([1, 0, 1], [1, 2])
    result.evidence().remainder.coefficients[0] = QQ(Rational(1))
