"""Find a point inside a circle and rule out contradictory constraints."""
from poly import sparse, cad, algebraic
from std.testing import assert_equal, assert_true


def main() raises:
    var circle = sparse.asarray([((2, 0), 1), ((0, 2), 1), ((0, 0), -1)])
    var interior = cad.solve([cad.negative(circle)])
    assert_true(cad.check(interior))
    var point = cad.sample(interior)
    assert_true(Bool(point))
    var x = point.value()[0].copy()
    var y = point.value()[1].copy()
    var radius_squared = algebraic.add(algebraic.multiply(x, x), algebraic.multiply(y, y))
    assert_equal(algebraic.compare(radius_squared, 1), -1)
    print("Interior point found:", interior.satisfiable())

    var impossible = cad.solve([cad.negative(circle), cad.nonnegative(circle)])
    assert_true(cad.check(impossible))
    print("Contradictory constraints satisfiable:", impossible.satisfiable())
    var decomposition = cad.decompose([circle.copy()])
    assert_true(cad.check(decomposition))
    print("Circle decomposition cells:", len(cad.cells(decomposition)))
