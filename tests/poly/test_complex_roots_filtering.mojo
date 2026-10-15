"""Filter extraneous projections, including fibers with one zero equation."""
from apn_mojo import ExactComplex
from poly import dense, roots, algebraic
from std.testing import assert_equal, assert_true


def main() raises:
    # The projection also contains pairwise averages -1/2, 1/2, and 2.
    # Their imaginary-coordinate equations have no real solution.
    var p = [6, -5, -2, 1]
    var values = dense.complex_roots(p)
    assert_equal(len(values), 3)
    assert_true(values[0] == algebraic.ascomplex(-2))
    assert_true(values[1] == algebraic.ascomplex(1))
    assert_true(values[2] == algebraic.ascomplex(3))
    assert_true(roots.check(values))
    # At real coordinate zero, Re(p(x+iy)) vanishes identically in y.
    var axis = dense.complex_roots([0, -1, 0, 1])
    assert_equal(len(axis), 3)
    assert_true(axis[0] == algebraic.ascomplex(-1))
    assert_true(axis[1].is_zero())
    assert_true(axis[2] == algebraic.ascomplex(1))
    assert_true(roots.check(axis))
    var scaled = dense.multiply([ExactComplex(2, 3)], [0, -1, 0, 1])
    var scaled_values = dense.complex_roots(scaled)
    assert_equal(len(scaled_values), 3)
    for i in range(3):
        assert_true(scaled_values[i] == axis[i])
        assert_equal(scaled_values.multiplicity(i), 1)
    assert_true(roots.check_complex(scaled, scaled_values.evidence()))
    var constant = dense.complex_roots([ExactComplex(2, 3)])
    assert_equal(len(constant), 0)
    assert_true(roots.check(constant))
    print("array complex root filtering checks passed")
