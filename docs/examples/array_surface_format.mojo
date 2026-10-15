"""Give variables names when displaying coefficient-based polynomials."""
from poly import dense, sparse


def main() raises:
    print(dense.format([1, -2, 3], variable="t"))
    var p = sparse.asarray([((2, 1), 3), ((0, 0), -1)])
    print(sparse.format(p, variables=["u", "v"]))
