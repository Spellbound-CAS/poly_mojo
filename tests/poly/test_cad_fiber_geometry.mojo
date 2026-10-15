"""Exact coefficient fibers, independent fields, dependencies, and axis order."""
from poly import cad, dense, sparse, algebraic
from polynomials import Budget
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var b = dense.root_between([-3, 0, 1], 1, 2)
    var p = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    var q = [(0, algebraic.negative(b)), (1, algebraic.asreal(1))]
    var result = cad.decompose([p.copy(), q.copy()])
    assert_true(cad.check(result))
    var cells = cad.cells(result)
    # sqrt(2) < sqrt(3): five cells on the line, with no conjugate branches.
    assert_equal(len(cells), 5)
    assert_equal(cells[0].signs, [-1, -1])
    assert_equal(cells[1].signs, [0, -1])
    assert_equal(cells[2].signs, [1, -1])
    assert_equal(cells[3].signs, [1, 0])
    assert_equal(cells[4].signs, [1, 1])
    assert_true(cells[1].sample[0] == a and cells[3].sample[0] == b)
    assert_equal(len(result.evidence().stacks), 3)
    assert_equal(len(result.evidence().cells), 8)
    var budget = Budget()
    var checker = cad.checker([p.copy(), q.copy()], budget)
    var receipt = cad.decompose_to([p.copy(), q.copy()], checker, budget)
    assert_true(checker.finish(receipt, budget))
    assert_equal(receipt.stacks, 3)
    assert_equal(receipt.cells, 5)
    # Opposite and repeated values share an auxiliary coordinate.
    var opposite = [(0, a.copy()), (1, algebraic.asreal(1))]
    var repeated = cad.decompose([p.copy(), opposite.copy(), p.copy()])
    assert_true(cad.check(repeated))
    assert_equal(len(repeated.evidence().projection.variables), 2)
    assert_equal(len(cad.cells(repeated)), 5)
    # A nonminimal witness and its dependent double must select compatible roots.
    var witness = dense.root_between([-2, 0, -1, 0, 1], 1, 2)
    var double = algebraic.multiply(witness, 2)
    var dependent = [(0, algebraic.negative(double)), (1, algebraic.asreal(1))]
    var related = cad.decompose([p.copy(), dependent.copy()])
    assert_true(cad.check(related))
    assert_equal(len(related.evidence().stacks), 3)
    assert_equal(len(cad.cells(related)), 5)
    var x = [((1, 0), algebraic.asreal(1)), ((0, 0), algebraic.negative(a))]
    var y = [((0, 1), algebraic.asreal(1)), ((0, 0), algebraic.negative(b))]
    var plane = cad.decompose([x.copy(), y.copy()], axes=[1, 0])
    assert_true(cad.check(plane))
    var plane_cells = cad.cells(plane)
    assert_equal(len(plane_cells), 9)
    assert_equal(len(plane.evidence().stacks), 6)
    assert_equal(len(plane.evidence().cells), 15)
    var point_count = 0
    for cell in plane_cells:
        if cell.dimension == 0:
            point_count += 1
            assert_true(cell.sample[0] == a and cell.sample[1] == b)
            assert_equal(cell.signs, [0, 0])
    assert_equal(point_count, 1)
    var rational = [(0, algebraic.asreal(-1)), (1, algebraic.asreal(1))]
    var empty_prefix = cad.decompose([rational.copy()], lazard=False)
    assert_true(cad.check(empty_prefix))
    assert_equal(len(empty_prefix.evidence().stacks), 1)
    var mccallum = cad.decompose([p.copy()], lazard=False)
    assert_true(cad.check(mccallum))
    assert_equal(len(cad.cells(mccallum)), 3)
    print("array CAD fiber geometry checks passed")
