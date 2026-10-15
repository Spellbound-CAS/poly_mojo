from poly import dense, rational, algebraic


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var f = rational.asarray([1], [-2, 0, 1])
    var expansion = rational.laurent(f, a, terms=3)
    print("valuation:", expansion.valuation)
    print("checked:", rational.check(expansion))
    print("residue:", rational.residue(f, a) == algebraic.reciprocal(algebraic.multiply(a, 2)))
