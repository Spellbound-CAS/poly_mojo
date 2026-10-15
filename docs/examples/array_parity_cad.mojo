"""Select a CAD cell and recover a nullified Lazard lifting polynomial."""
from poly import cad, algebraic, dense
from std.testing import assert_equal, assert_true


def main() raises:
    var inputs = [[((1, 0), 1), ((0, 0), -2)], [((0, 1), 1), ((0, 0), 3)]]
    var projection = cad.project(inputs, axes=[1, 0])
    assert_true(cad.check(projection))
    var cell = cad.cell_at(inputs, (2, -3), axes=[1, 0])
    assert_true(cad.check(cell))
    assert_equal(cad.signs(cell), [0, 0])
    var point = cad.sample(cell)
    assert_equal(algebraic.compare(point[0], 2), 0)
    assert_equal(algebraic.compare(point[1], -3), 0)
    var lifting = cad.lazard_evaluate([((2, 1), 1), ((2, 0), -1)], [(0, 0)], axis=1)
    assert_true(cad.check(lifting))
    assert_equal(cad.vanishing_orders(lifting), [(0, 2)])
    var roots = dense.real_roots(cad.lifting_coefficients(lifting))
    assert_equal(len(roots), 1)
    assert_equal(algebraic.compare(roots[0], 1), 0)
    print("Selected cell and Lazard lifting verified")
