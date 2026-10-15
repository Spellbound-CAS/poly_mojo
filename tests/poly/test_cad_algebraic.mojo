"""Algebraic CAD cells, mixed constraints, witnesses, and original coordinate shapes."""
from poly import cad, dense, sparse, algebraic
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var one = algebraic.asreal(1)
    var p = [(0, algebraic.negative(a)), (1, one.copy())]
    var decomposition = cad.decompose([p.copy()])
    assert_true(cad.check(decomposition))
    var cells = cad.cells(decomposition)
    assert_equal(len(cells), 3)
    assert_equal(cells[0].signs, [-1])
    assert_equal(cells[1].signs, [0])
    assert_equal(cells[2].signs, [1])
    assert_equal(cells[0].dimension, 1)
    assert_equal(cells[1].dimension, 0)
    assert_equal(len(cells[1].sample), 1)
    assert_true(cells[1].sample[0] == a)
    var exact = cad.solve([cad.zero(p), cad.positive([(1, 1)])])
    assert_true(exact.satisfiable() and cad.check(exact))
    assert_true(cad.sample(exact).value()[0] == a)
    var impossible = cad.solve([cad.negative(p), cad.nonnegative(p)])
    assert_false(impossible.satisfiable())
    assert_true(cad.check(impossible))
    assert_false(Bool(cad.sample(impossible)))
    var covering = cad.cover([cad.nonzero(p), cad.zero(p)])
    assert_false(covering.satisfiable())
    assert_true(cad.check(covering))
    assert_equal(len(cad.cells(covering)), 3)
    var feasible = cad.cover([cad.positive(p)])
    assert_true(feasible.satisfiable() and cad.check(feasible))
    assert_equal(algebraic.compare(cad.sample(feasible).value()[0], a), 1)
    var x = [((1, 0), one.copy()), ((0, 0), algebraic.asreal(-1))]
    var y = [((0, 1), one.copy()), ((0, 0), algebraic.negative(a))]
    var reordered = cad.solve([cad.zero(x), cad.zero(y)], axes=[1, 0])
    assert_true(cad.check(reordered))
    var point = cad.sample(reordered).value().copy()
    assert_equal(len(point), 2)
    assert_true(point[0] == one and point[1] == a)
    var constant = [((), a.copy())]
    var scalar = cad.decompose([constant.copy()])
    assert_true(cad.check(scalar))
    assert_equal(len(cad.cells(scalar)), 1)
    assert_equal(len(cad.cells(scalar)[0].sample), 0)
    assert_equal(cad.cells(scalar)[0].dimension, 0)
    assert_equal(cad.cells(scalar)[0].signs, [1])
    var zero = List[Tuple[Int, algebraic.Real]]()
    var zero_result = cad.decompose([zero.copy()])
    assert_true(cad.check(zero_result))
    assert_equal(cad.cells(zero_result)[0].signs, [0])
    print("array algebraic CAD checks passed")
