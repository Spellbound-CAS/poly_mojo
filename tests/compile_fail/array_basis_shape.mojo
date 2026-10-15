# expect: every polynomial must use the same power shape and exact coefficients
from poly import sparse

def main() raises:
    _ = sparse.groebner(([(1, 1)], [((1, 0), 1)]))
