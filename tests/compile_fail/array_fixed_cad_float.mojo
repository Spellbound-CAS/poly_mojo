# expect: real algebraic values require exact integers or rationals
from poly import cad

def main() raises:
    _ = cad.decompose([[(1, 1)]], fixed=[(0, Float64(0.5))])
