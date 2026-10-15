"""Sign constraints, exact samples, and complete one-variable cells."""
from poly import cad, algebraic
from std.testing import assert_equal, assert_true


def main() raises:
    var p = [(0, -1), (2, 1)]
    var result = cad.decompose([p.copy()])
    assert_true(cad.check(result))
    var cells = cad.cells(result)
    assert_equal(len(cells), 5)
    assert_equal(cells[0].signs[0], 1)
    assert_equal(cells[1].signs[0], 0)
    assert_equal(cells[2].signs[0], -1)
    assert_equal(cells[3].signs[0], 0)
    assert_equal(cells[4].signs[0], 1)
    assert_equal(cells[1].dimension, 0)
    assert_equal(cells[2].dimension, 1)
    var inside = cad.solve([cad.negative(p)])
    assert_true(inside.satisfiable())
    assert_true(cad.check(inside))
    var point = cad.sample(inside)
    assert_true(Bool(point))
    assert_equal(algebraic.compare(point.value()[0], -1), 1)
    assert_equal(algebraic.compare(point.value()[0], 1), -1)
    var impossible = cad.solve([cad.negative(p), cad.nonnegative(p)])
    assert_true(not impossible.satisfiable())
    assert_true(cad.check(impossible))
    assert_true(not cad.sample(impossible))
    var boundary = cad.solve([cad.zero(p)])
    assert_true(boundary.satisfiable() and cad.check(boundary))
    var boundary_point = cad.sample(boundary)
    var boundary_root = boundary_point.value()[0].copy()
    assert_equal(algebraic.compare(algebraic.multiply(boundary_root, boundary_root), 1), 0)
    var nonzero = cad.solve([cad.nonzero(p), cad.nonpositive(p)])
    assert_true(nonzero.satisfiable() and cad.check(nonzero))
    var positive = cad.solve([cad.positive(p)])
    assert_true(positive.satisfiable() and cad.check(positive))
    var x = [((1, 0), 1), ((0, 0), -2)]
    var y = [((0, 1), 1), ((0, 0), 3)]
    var ordered = cad.solve([cad.zero(x), cad.zero(y)], axes=[1, 0])
    assert_true(cad.check(ordered))
    var ordered_point = cad.sample(ordered)
    assert_equal(algebraic.compare(ordered_point.value()[0], 2), 0)
    assert_equal(algebraic.compare(ordered_point.value()[1], -3), 0)
    var circle = [((2, 0), 1), ((0, 2), 1), ((0, 0), -1)]
    var plane = cad.decompose([circle.copy()])
    assert_equal(len(cad.cells(plane)), 13)
    assert_true(cad.check(plane))
    var empty = List[cad.Constraint[Int]]()
    var unconstrained = cad.solve(empty)
    assert_true(unconstrained.satisfiable() and cad.check(unconstrained))
    assert_equal(len(cad.sample(unconstrained).value()), 1)
    print("array CAD checks passed")
