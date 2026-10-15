"""Root finding from exact algebraic coefficients and specialization."""
from std.testing import assert_equal, assert_true
from poly import dense, sparse, algebraic, roots


def main() raises:
    var sqrt2 = dense.root_between([-2, 0, 1], 1, 2)
    var coefficients = [algebraic.negative(sqrt2), algebraic.asreal(1)]
    var result = dense.real_roots(coefficients)
    assert_equal(len(result), 1)
    assert_equal(algebraic.compare(result[0], sqrt2), 0)
    assert_true(roots.check(result))
    var p = [((1, 0), 1), ((0, 1), -1)]
    var specialized = sparse.specialize(p, [(1, sqrt2.copy())])
    var after = sparse.real_roots(specialized, axis=0)
    assert_equal(len(after), 1)
    assert_equal(algebraic.compare(after[0], sqrt2), 0)
    assert_true(roots.check(after))
    assert_true(algebraic.check(after.specialization()))
    for value in after:
        assert_equal(algebraic.compare(value, sqrt2), 0)
    var direct = sparse.roots_at(p, [(1, sqrt2.copy())], axis=0)
    assert_true(roots.check(direct))
    assert_equal(algebraic.compare(direct[0], sqrt2), 0)
    var no_roots = dense.real_roots([algebraic.asreal(1)])
    assert_equal(len(no_roots), 0)
    assert_true(roots.check(no_roots))
    print("array extension algebraic_roots checks passed")
