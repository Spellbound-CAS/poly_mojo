"""Exact polynomial values and partial assignments at algebraic points."""
from std.testing import assert_equal, assert_true
from poly import dense, sparse, algebraic


def main() raises:
    var x = dense.root_between([-2, 0, 1], 1, 2)
    var result = dense.evaluate([0, 0, 1], x)
    assert_equal(algebraic.compare(result, 2), 0)
    assert_true(algebraic.check(result))
    var sp = [((2, 0), 1), ((0, 1), 3)]
    var evaluated = sparse.evaluate(sp, (x.copy(), 2))
    assert_equal(algebraic.compare(evaluated, 8), 0)
    assert_true(algebraic.check(evaluated))
    var specialized = sparse.specialize(sp, [(0, x.copy())])
    assert_equal(len(specialized), 2)
    assert_true(algebraic.check(specialized))
    for term in specialized:
        var expected = 2 if term[0] == (0, 0) else 3
        assert_equal(algebraic.compare(term[1], expected), 0)
    var remaining = algebraic.coefficients(specialized, axis=1)
    assert_equal(algebraic.compare(dense.evaluate(remaining, 2), 8), 0)
    var sparse_value = sparse.evaluate([(2, 1)], x)
    assert_equal(algebraic.compare(sparse_value, 2), 0)
    assert_true(algebraic.check(sparse_value))
    print("array gap algebraic checks passed")
