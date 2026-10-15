"""Algebraic coefficient prefixes, vanishing specialization, and empty conjunctions."""
from poly import cad, dense, sparse, algebraic
from polynomials import Budget
from std.testing import assert_true, assert_equal

def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = sparse.asarray((((1, 0), 1), ((0, 0), algebraic.negative(a))))
    var budget = Budget()
    var check = cad.checker([p.copy()], fixed=[(1, 0)])
    var receipt = cad.decompose_to([p.copy()], check, budget, fixed=[(1, 0)])
    assert_true(check.finish(receipt, budget))
    assert_equal(receipt.cells, 3)
    var solve = cad.solve([cad.zero(p)], fixed=[(1, 0)])
    assert_true(solve.satisfiable() and cad.check(solve))
    var cover = cad.cover([cad.zero(p)], fixed=[(1, 0)])
    assert_true(cover.satisfiable() and cad.check(cover))
    var nullified = [((2, 1), 1)]
    var result = cad.decompose([nullified.copy()], fixed=[(0, 0)])
    assert_true(cad.check(result))
    var cells = cad.cells(result)
    assert_equal(len(cells), 3)
    assert_equal(cells[0].dimension, 1)
    assert_equal(cells[1].dimension, 0)
    assert_equal(cells[2].dimension, 1)
    for cell in cells:
        assert_equal(cell.signs, [0])
    var empty = List[cad.Constraint[Tuple[Int, Int]]]()
    var point = cad.solve(empty, fixed=[(1, 1), (0, 2)])
    assert_true(point.satisfiable() and cad.check(point))
    var sample = cad.sample(point)
    assert_true(sample.value()[0] == algebraic.asreal(2))
    assert_true(sample.value()[1] == algebraic.asreal(1))
    print("fixed CAD edge checks passed")
