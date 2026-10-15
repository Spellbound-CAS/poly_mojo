"""Selected CAD workflows bind axis order, query points and Lazard orders."""
from std.testing import assert_true, assert_false, assert_equal
from apn_mojo import Rational
from polynomials import PolyErrorKind
from poly import cad, dense, sparse, algebraic, roots


def main() raises:
    var x = [((1, 0), 1), ((0, 0), -2)]
    var y = [((0, 1), 1), ((0, 0), 3)]
    var inputs = [x.copy(), y.copy()]
    var projection = cad.project(inputs, axes=[1, 0])
    assert_true(cad.check(projection))
    assert_true(len(cad.projection_polynomials(projection)) >= 2)
    assert_true(len(cad.projection_polynomials(projection, axis=0)) > 0)
    var selected = cad.cell_at(inputs, (2, Rational(-3)), axes=[1, 0])
    assert_true(cad.check(selected))
    assert_equal(cad.signs(selected), [0, 0])
    var point = cad.sample(selected)
    assert_equal(algebraic.compare(point[0], 2), 0)
    assert_equal(algebraic.compare(point[1], -3), 0)
    var interior = cad.cell_at(inputs, [3, -2], axes=[1, 0], lazard=False)
    assert_true(cad.check(interior))
    assert_equal(cad.signs(interior), [1, 1])
    var sqrt2 = dense.real_roots([-2, 0, 1])[1].copy()
    var root_cell = cad.cell_at([[(0, -2), (2, 1)]], [sqrt2.copy()])
    assert_true(cad.check(root_cell))
    assert_equal(cad.signs(root_cell), [0])
    # x^2(y-1) vanishes at x=0. The least nonzero derivative is 2(y-1).
    var lazard = cad.lazard_evaluate([((2, 1), 1), ((2, 0), -1)], [(0, 0)], axis=1)
    assert_true(cad.check(lazard))
    assert_equal(cad.vanishing_orders(lazard), [(0, 2)])
    var coefficients = cad.lifting_coefficients(lazard)
    assert_equal(algebraic.compare(coefficients[0], -2), 0)
    assert_equal(algebraic.compare(coefficients[1], 2), 0)
    var lifted = dense.real_roots(coefficients)
    assert_true(roots.check(lifted))
    assert_equal(len(lifted), 1)
    assert_equal(algebraic.compare(lifted[0], 1), 0)
    var base = cad.lazard_evaluate([(0, -2), (2, 1)], List[Tuple[Int, Int]]())
    assert_true(cad.check(base))
    assert_equal(len(cad.vanishing_orders(base)), 0)
    assert_equal(len(cad.lifting_coefficients(base)), 3)
    var p = [(0, -1), (2, 1)]
    var feasible = cad.cover([cad.negative(p)])
    assert_true(feasible.satisfiable() and cad.check(feasible))
    assert_equal(len(cad.cells(feasible)), 5)
    assert_equal(cad.cells(feasible)[2].signs, [-1])
    var sample = cad.sample(feasible)
    assert_equal(algebraic.compare(sample.value()[0], -1), 1)
    assert_equal(algebraic.compare(sample.value()[0], 1), -1)
    var impossible = cad.cover([cad.negative(p), cad.nonnegative(p)])
    assert_false(impossible.satisfiable())
    assert_true(cad.check(impossible))
    assert_false(Bool(cad.sample(impossible)))
    projection._data[].proof.variables[0] = 0
    assert_false(cad.check(projection))
    selected._data[].proof.signs[0].relation = 1
    assert_false(cad.check(selected))
    lazard._data[].proof.orders[0] = 1
    assert_false(cad.check(lazard))
    var failed = False
    try:
        _ = cad.cell_at(inputs, [1])
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    failed = False
    try:
        _ = cad.lazard_evaluate(x, [(0, 1)], axis=0)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    print("array parity CAD checks passed")
