# expect: expression must be mutable in assignment
from poly import rational


def main() raises:
    var p = rational.asarray([1], [1])
    p.coprimality_evidence().root = -1
