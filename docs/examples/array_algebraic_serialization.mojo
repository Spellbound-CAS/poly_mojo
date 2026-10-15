from poly import dense, sparse, algebraic, coefficients


def main() raises:
    var root = dense.root_between([-2, 0, 1], 1, 2)
    var original = dense.asarray(coefficients((root.copy(), 0, 1)))
    var restored = dense.loads[algebraic.Real](dense.dumps(original))
    print(restored.domain(), restored[0] == root, algebraic.check(restored[0]))
    var z = algebraic.ascomplex(root, 1)
    var terms = sparse.asarray((((1, 0), z.copy()), ((0, 2), 1)))
    var saved = sparse.dumps(terms)
    var loaded = sparse.loads[algebraic.Complex, ndim=2](saved)
    print(loaded.domain(), sparse.dumps(loaded) == saved)
