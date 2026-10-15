# expect: expression must be mutable in assignment
from poly import dense, algebraic, coefficients

def main() raises:
    var result = dense.factor_checked(coefficients((algebraic.asreal(-2), 0, 1)))
    result.evidence().factors[0][1] = 2
