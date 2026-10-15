"""Complete root solving over real and complex algebraic coefficient fields."""
from poly import dense, sparse, roots, algebraic, coefficients
from std.testing import assert_equal, assert_true


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var linear = dense.complex_roots(coefficients((algebraic.negative(a), 1)))
    assert_equal(len(linear), 1)
    assert_true(linear[0] == algebraic.ascomplex(a))
    assert_true(roots.check(linear))
    var p = dense.asarray(coefficients((a.copy(), 0, 1)))
    var r = dense.complex_roots(p)
    var fourth_root = dense.root_between([-2, 0, 0, 0, 1], 1, 2)
    assert_equal(len(r), 2)
    assert_true(r[0] == algebraic.ascomplex(0, algebraic.negative(fourth_root)))
    assert_true(r[1] == algebraic.ascomplex(0, fourth_root))
    assert_true(roots.check(r))
    var z = algebraic.ascomplex(a, 1)
    var w = algebraic.ascomplex(a, -2)
    var complex_p = dense.from_roots([z.copy(), w.copy()])
    var complex_r = dense.complex_roots(complex_p)
    assert_equal(len(complex_r), 2)
    assert_true(complex_r[0] == w)
    assert_true(complex_r[1] == z)
    assert_true(roots.check(complex_r))
    var terms = sparse.asarray(((2, 1), (0, a.copy())))
    var sr = sparse.complex_roots(terms)
    assert_true(sr[0] == r[0])
    assert_true(sr[1] == r[1])
    print("array algebraic complex root checks passed")
