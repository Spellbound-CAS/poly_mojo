"""Array algebraic points."""
from poly import dense, sparse, algebraic


def main() raises:
    var root = dense.root_between([-2, 0, 1], 1, 2)
    var value = dense.evaluate([0, 0, 1], root)
    print("Square equals two:", algebraic.compare(value, 2) == 0)
    var p = [((2, 0), 1), ((0, 1), 3)]
    var specialized = sparse.specialize(p, [(0, root.copy())])
    var remaining = algebraic.coefficients(specialized, axis=1)
    print("Value at y=2 equals eight:", algebraic.compare(dense.evaluate(remaining, 2), 8) == 0)
    print("Specialization checked:", algebraic.check(specialized))
