"""CAD outputs refine every irrational coordinate to width at most 2^-830."""
from apn_mojo import Integer, Rational
from poly import cad, algebraic
from std.testing import assert_true


def main() raises:
    var floor = Rational(Integer(1), Integer(1) << Integer(830))
    var circle = [((2, 0), 1), ((0, 2), 1), ((0, 0), -2)]
    var cubic = [((0, 2), 1), ((3, 0), -1), ((0, 0), -2)]
    var irrational = 0
    var one = cad.decompose([circle.copy()])
    var two = cad.decompose([circle.copy(), cubic.copy()])
    assert_true(cad.check(one) and cad.check(two))
    for cell in cad.cells(one) + cad.cells(two):
        for value in cell.sample:
            var b = algebraic.bounds(value)
            if b[0] != b[1]:
                irrational += 1
                assert_true(b[1] - b[0] <= floor)
    assert_true(irrational > 0)
    var solved = cad.solve([cad.zero(circle), cad.positive(cubic)])
    assert_true(cad.check(solved))
    for value in cad.sample(solved).value():
        var b = algebraic.bounds(value)
        assert_true(b[1] - b[0] <= floor)
    print("CAD output precision checks passed")
