"""Array serialization."""
from apn_mojo import Integer, Rational
from poly import dense, sparse


def main() raises:
    var p = dense.asarray([1, -2, 3])
    var saved = dense.dumps(p)
    var restored = dense.loads[Integer](saved)
    print("Restored value at two:", dense.evaluate(restored, 2))
    var q = sparse.asarray([((1, 0), Rational(1, 3)), ((0, 2), Rational(2))])
    var terms = sparse.loads[Rational, ndim=2](sparse.dumps(q))
    print("Sparse value:", sparse.evaluate(terms, (3, 2)))
