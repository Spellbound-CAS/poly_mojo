"""Compact CAD archives and genuine revision-14 archives share accepted semantics."""
from poly import cad, io, algebraic
from polynomials import Budget, PolyErrorKind
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var quadratic = io.loads[cad.Decomposition[Int]](open("tests/fixtures/cad-revision14/quadratic.txt", "r").read())
    assert_true(cad.check(quadratic))
    assert_false(quadratic.evidence().is_compact())
    assert_equal(len(cad.cells(quadratic)), 5)
    var impossible = io.loads[cad.Covering[Int]](open("tests/fixtures/cad-revision14/unsatisfiable.txt", "r").read())
    assert_true(cad.check(impossible))
    assert_false(impossible.satisfiable())
    var fixed = io.loads[cad.Decomposition[Tuple[Int, Int]]](open("tests/fixtures/cad-revision14/fixed.txt", "r").read())
    assert_true(cad.check(fixed))
    assert_equal(len(cad.cells(fixed)), 3)
    assert_true(cad.fixed_coordinates(fixed)[0][1] == algebraic.asreal(2))
    var aa = io.loads[cad.Decomposition[Int]](open("tests/fixtures/cad-revision14/algebraic.txt", "r").read())
    assert_true(cad.check(aa))
    assert_equal(len(cad.cells(aa)), 3)
    # Old evidence is retained, checked, and written with the new explicit tag.
    var migrated = io.dumps(quadratic)
    assert_true(migrated.find("CHECKED-3:") >= 0)
    var again = io.loads[type_of(quadratic)](migrated)
    assert_equal(io.dumps(again), migrated)
    var current = cad.decompose([[(0, -2), (2, 1)]])
    assert_true(current.evidence().is_compact())
    var text = io.dumps(current)
    var restored = io.loads[type_of(current)](text)
    assert_true(restored.evidence().is_compact() and cad.check(restored))
    assert_equal(io.dumps(restored), text)
    var old_cells = cad.cells(quadratic)
    var new_cells = cad.cells(restored)
    for i in range(len(old_cells)):
        assert_equal(old_cells[i].dimension, new_cells[i].dimension)
        assert_equal(old_cells[i].signs, new_cells[i].signs)
        assert_true(old_cells[i].sample[0] == new_cells[i].sample[0])
    var covering = cad.cover([cad.zero([(0, -2), (2, 1)])])
    var restored_cover = io.loads[type_of(covering)](io.dumps(covering))
    assert_true(restored_cover.satisfiable() and cad.check(restored_cover))
    assert_true(cad.sample(restored_cover).value()[0] == cad.sample(covering).value()[0])
    # The new tag cannot be used as an old archive envelope.
    var bad = text.replace("CHECKED-3:", "CHECKED-1:")
    var rejected = False
    try:
        _ = io.loads[type_of(current)](bad, max_steps=1000000)
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING or e.kind == PolyErrorKind.LIMIT or e.kind == PolyErrorKind.CAPACITY
    assert_true(rejected)
    print("final CAD archive compatibility checks passed")
