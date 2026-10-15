"""Decompose, solve, stream, and save a CAD slice at exact coordinates."""
from poly import cad, dense, algebraic, io
from polynomials import Budget
from std.testing import assert_true, assert_equal


def main() raises:
    # x - y, with y fixed at sqrt(2).
    var p = [((1, 0), 1), ((0, 1), -1)]
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var fixed = [(1, a.copy())]
    var result = cad.decompose([p.copy()], fixed=fixed)
    assert_true(cad.check(result))
    var cells = cad.cells(result)
    assert_equal(len(cells), 3)
    assert_equal(cells[1].signs, [0])
    assert_true(cells[1].sample[0] == a and cells[1].sample[1] == a)
    print("Cells on y = sqrt(2):", len(cells))

    var solved = cad.solve([cad.zero(p)], fixed=fixed)
    assert_true(solved.satisfiable() and cad.check(solved))
    var point = cad.sample(solved)
    assert_true(point.value()[0] == a)
    print("Exact satisfying point checked")

    var budget = Budget()
    var checker = cad.checker([p.copy()], budget, fixed=fixed)
    var receipt = cad.decompose_to([p.copy()], checker, budget, fixed=fixed)
    assert_true(checker.finish(receipt, budget))
    assert_equal(receipt.cells, 3)
    print("Restricted stream checked")

    var restored = io.loads[type_of(result)](io.dumps(result))
    assert_true(cad.check(restored))
    assert_equal(cad.fixed_coordinates(restored)[0][0], 1)
    assert_true(cad.fixed_coordinates(restored)[0][1] == a)
    print("Saved coordinate binding checked")
