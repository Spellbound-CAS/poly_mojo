"""Reject missing branches, wrong coefficient roots, and incomplete conflict covers."""
from poly import cad, dense, algebraic
from polynomials import Budget, check_cad, CADFiber, RealAlgebraic, PolyErrorKind
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    var original = cad.decompose([p.copy()])
    var budget = Budget()
    # The raw restricted tree is not evidence for all auxiliary coefficient values.
    assert_false(check_cad(original._data[].inputs, original.evidence(), budget))
    var missing = cad.decompose([p.copy()])
    _ = missing._data[].proof.cells.pop()
    assert_false(cad.check(missing))
    var missing_stack = cad.decompose([p.copy()])
    _ = missing_stack._data[].proof.stacks.pop()
    assert_false(cad.check(missing_stack))
    var wrong = cad.decompose([p.copy()])
    wrong._data[].proof.cells[1].sample = a.witness().copy()
    assert_false(cad.check(wrong))
    var sector = cad.decompose([p.copy()])
    sector._data[].proof.cells[1].position = 0
    assert_false(cad.check(sector))
    var duplicate = cad.decompose([p.copy()])
    duplicate._data[].proof.cells.append(duplicate._data[].proof.cells[2].copy())
    assert_false(cad.check(duplicate))
    var bad_lift = cad.decompose([p.copy()])
    _ = bad_lift._data[].proof.stacks[0].ordering.pop()
    assert_false(cad.check(bad_lift))
    var alien = CADFiber([a.witness().copy()], budget)
    assert_false(alien.check(original._data[].inputs, original.evidence(), budget))
    var too_long = CADFiber([a.witness().copy(), a.witness().copy(), a.witness().copy()], budget)
    assert_false(too_long.check(original._data[].inputs, original.evidence(), budget))
    var absent = CADFiber([algebraic.asreal(100).witness().copy()], budget)
    var rejected = False
    try:
        _ = absent.decompose(original._data[].inputs, original.evidence().projection.variables, budget, lazard=True)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    var impossible = cad.solve([cad.negative(p), cad.nonnegative(p)])
    assert_false(impossible.satisfiable())
    assert_true(cad.check(impossible))
    for conflict in impossible.evidence().conflicts:
        assert_true(impossible.evidence().cells[conflict.parent].level >= 0)
    _ = impossible._data[].proof.conflicts.pop()
    assert_false(cad.check(impossible))
    var covering = cad.cover([cad.negative(p), cad.nonnegative(p)])
    assert_true(cad.check(covering))
    assert_equal(len(covering.evidence().violations), 3)
    _ = covering._data[].proof.violations.pop()
    assert_false(cad.check(covering))
    # A contradiction at a coefficient section can close the fiber immediately.
    var constant = [((), a.copy())]
    var scalar = cad.solve([cad.zero(constant)])
    assert_false(scalar.satisfiable())
    assert_true(cad.check(scalar))
    assert_equal(len(scalar.evidence().conflicts), 1)
    scalar._data[].proof.conflicts[0].lower = 0
    assert_false(cad.check(scalar))
    print("array CAD fiber evidence checks passed")
