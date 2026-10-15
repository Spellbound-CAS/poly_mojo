# expect: expression must be mutable in assignment
from poly import roots
from polynomials import QQ
from apn_mojo import Rational

def main() raises:
    var result = roots.diagnostics([1, 2], 0)
    result.evidence().backward_error_high = QQ(Rational(0))
