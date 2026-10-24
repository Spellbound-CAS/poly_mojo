"""Section sign witnesses survive an archive round trip and still check."""
from poly import cad, io
from polynomials import Budget
from std.testing import assert_true, assert_equal


def main() raises:
    var budget = Budget()
    var circle = [((2, 0), 1), ((0, 2), 1), ((0, 0), -1)]
    var parabola = [((0, 1), 2), ((2, 0), -2), ((0, 0), 1)]
    var inputs = [circle.copy(), parabola.copy()]
    for mode in range(2):
        var result = cad.decompose(inputs, lazard=Bool(mode))
        assert_true(cad.check(result))
        var sections = 0
        if result.evidence().is_compact():
            for cell in result.evidence().cells(budget):
                if cell.position % 2 and len(cell.signs):
                    for sign in cell.signs:
                        sections += Int(sign == 0)
        else:
            for cell in result.evidence().legacy_decomposition().value().cells:
                for w in cell.witnesses:
                    sections += Int(w.kind == 4)
        assert_true(sections > 0)
        var text = io.dumps(result)
        var restored = io.loads[type_of(result)](text)
        assert_equal(io.dumps(restored), text)
        assert_true(cad.check(restored))
        assert_equal(len(cad.cells(restored)), len(cad.cells(result)))
    print("section witness archive checks passed")
