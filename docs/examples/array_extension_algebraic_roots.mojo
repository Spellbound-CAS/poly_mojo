"""Find exact roots with an algebraic coefficient or parameter."""
from poly import dense, sparse, algebraic, roots
from std.testing import assert_equal, assert_true


def main() raises:
    var sqrt2 = dense.root_between([-2, 0, 1], 1, 2)
    var coefficients = [algebraic.negative(sqrt2), algebraic.asreal(1)]
    var answer = dense.real_roots(coefficients)
    assert_equal(len(answer), 1)
    assert_equal(algebraic.compare(answer[0], sqrt2), 0)
    assert_true(roots.check(answer))
    var surface = [((1, 0), 1), ((0, 1), -1)]
    var specialized = sparse.specialize(surface, [(1, sqrt2.copy())])
    var remaining = sparse.real_roots(specialized, axis=0)
    assert_true(roots.check(remaining))
    assert_equal(algebraic.compare(remaining[0], sqrt2), 0)
    var direct = sparse.roots_at(surface, [(1, sqrt2.copy())], axis=0)
    assert_true(roots.check(direct))
    print("x - sqrt(2): one exact root")
    print("x - y at y = sqrt(2): specialization and complete roots checked")
