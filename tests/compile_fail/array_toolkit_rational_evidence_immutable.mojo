# expect: expression must be mutable in assignment
from poly import rational
from polynomials import QQ
from apn_mojo import Rational

def main() raises:
    var p = rational.asarray([1], [1])
    p.evidence().common.terms[0].coefficient = QQ(Rational(2))
