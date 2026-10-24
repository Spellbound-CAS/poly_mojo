"""Coverings evidence of cad.solve: archive tag, round trips, and revision-9 archives."""
from poly import cad, io
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_true, assert_equal


def main() raises:
    var p = [((2, 0), 1), ((0, 2), 1), ((0, 0), -1)]
    var line = [((1, 0), 1), ((0, 0), -2)]
    var inside = cad.solve([cad.negative(p)])
    var outside = cad.solve([cad.zero(p), cad.zero(line)])
    assert_true(inside.satisfiable() and not outside.satisfiable())
    assert_true(cad.check(inside) and cad.check(outside))
    assert_equal(len(inside.evidence().witness), 2)
    assert_true(len(outside.evidence().nodes) > 0)
    for solution in [inside.copy(), outside.copy()]:
        var text = io.dumps(solution)
        assert_true(text.find("\nCHECKED-2:poly.cad.Solution[") >= 0)
        var restored = io.loads[type_of(solution)](text)
        assert_equal(io.dumps(restored), text)
        # An archive under the revision-9 tag fails cleanly.
        var old = text.replace("\nCHECKED-2:", "\nCHECKED-1:")
        var kind = PolyErrorKind.LIMIT
        try:
            _ = io.loads[type_of(solution)](old)
        except e:
            kind = e.kind
        assert_true(kind == PolyErrorKind.ENCODING)
    print("coverings archive checks passed")
