# expect: expression must be mutable in assignment
from poly import dense

def main() raises:
    var values = dense.complex_roots([1, 0, 1])
    values.evidence().multiplicities[0] = 2
