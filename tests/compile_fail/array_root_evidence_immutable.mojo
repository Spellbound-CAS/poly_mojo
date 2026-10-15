# expect: expression must be mutable in assignment
from poly import dense

def main() raises:
    var r = dense.count_roots([-2, 0, 1], -2, 2)
    r.evidence().count = 99
