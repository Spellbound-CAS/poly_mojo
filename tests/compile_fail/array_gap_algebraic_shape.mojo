# expect: evaluation point has the wrong number of coordinates
from poly import sparse, algebraic

def main() raises:
    _ = sparse.evaluate([((1, 0), 1)], (algebraic.asreal(1),))
