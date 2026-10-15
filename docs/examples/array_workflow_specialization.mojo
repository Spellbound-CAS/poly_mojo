"""Retain and archive the full calculation from AA specialization to exact roots."""
from poly import dense, sparse, algebraic, roots, io
from std.testing import assert_equal, assert_true

def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = sparse.asarray([((1, 0), a.copy()), ((0, 1), a.copy())])
    var result = sparse.roots_at(p, [(1, a.copy())], axis=0)
    assert_true(roots.check(result))
    assert_true(algebraic.check(result.specialization()))
    assert_equal(algebraic.compare(result[0], algebraic.negative(a)), 0)
    var restored = io.loads[type_of(result)](io.dumps(result))
    assert_true(roots.check(restored))
    print("sqrt(2)*x + 2 has root -sqrt(2); specialization and roots checked")
