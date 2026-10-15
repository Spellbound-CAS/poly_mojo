"""Exact roots, multiplicity blocks, and Gaussian rational inputs."""
from apn_mojo import Rational, ExactComplex
from poly import dense, sparse, roots, algebraic, coefficients
from std.testing import assert_equal, assert_true


def main() raises:
    var unit = dense.complex_roots([7])
    assert_equal(len(unit), 0)
    assert_true(roots.check(unit))
    var linear = dense.complex_roots([6, 3, 0])
    assert_equal(len(linear), 1)
    assert_true(linear[0] == algebraic.ascomplex(-2))
    assert_equal(linear.multiplicity(0), 1)
    var r = dense.complex_roots([1, 0, 1])
    assert_equal(len(r), 2)
    assert_true(r[0] == algebraic.ascomplex(0, -1))
    assert_true(r[1] == algebraic.ascomplex(0, 1))
    assert_true(roots.check(r))
    var repeated = dense.multiply(dense.power([1, 0, 1], 3), [-2, 1])
    var rr = dense.complex_roots(repeated)
    assert_equal(len(rr), 3)
    assert_true(rr[0] == r[0])
    assert_true(rr[1] == r[1])
    assert_true(rr[2] == algebraic.ascomplex(2))
    assert_equal(rr.multiplicity(0), 3)
    assert_equal(rr.multiplicity(1), 3)
    assert_equal(rr.multiplicity(2), 1)
    assert_true(roots.check(rr))
    var p = [ExactComplex(Rational(-1), Rational(5)), ExactComplex(Rational(-3), Rational(-4)), ExactComplex(1)]
    var gaussian = dense.complex_roots(p)
    assert_equal(len(gaussian), 2)
    assert_true(gaussian[0] == algebraic.ascomplex(1, 1))
    assert_true(gaussian[1] == algebraic.ascomplex(2, 3))
    assert_true(roots.check_complex(p, gaussian.evidence()))
    var rebuilt = dense.from_roots(gaussian)
    for i in range(len(p)):
        assert_true(rebuilt[i] == algebraic.ascomplex(p[i]))
    var sparse_roots = sparse.complex_roots([((0, 0), 1), ((0, 2), 1)], axis=1)
    assert_equal(len(sparse_roots), 2)
    assert_true(sparse_roots[0] == r[0])
    assert_true(roots.check(sparse_roots))
    var mult_zero = dense.complex_roots([0, 0, 0, 1])
    assert_equal(len(mult_zero), 1)
    assert_true(mult_zero[0].is_zero())
    assert_equal(mult_zero.multiplicity(0), 3)
    print("array exact complex root checks passed")
