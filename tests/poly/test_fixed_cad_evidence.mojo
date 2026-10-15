"""Reject proofs for another fixed value, axis, or incomplete slice."""
from poly import cad, dense, algebraic
from polynomials import Budget, PolyErrorKind, check_cad
from std.testing import assert_true, assert_false


def main() raises:
    var p = [((1, 0), 1), ((0, 1), -1)]
    var original = cad.decompose([p.copy()], fixed=[(1, 2)])
    var budget = Budget()
    assert_false(check_cad(original._data[].inputs, original.evidence(), budget))
    var wrong_value = cad.decompose([p.copy()], fixed=[(1, 2)])
    wrong_value._data[].fixed[0][1] = algebraic.asreal(3)
    assert_false(cad.check(wrong_value))
    var wrong_axis = cad.decompose([p.copy()], fixed=[(1, 2)])
    wrong_axis._data[].fixed[0][0] = 0
    assert_false(cad.check(wrong_axis))
    var missing = cad.decompose([p.copy()], fixed=[(1, 2)])
    _ = missing._data[].proof.cells.pop()
    assert_false(cad.check(missing))
    var unbound = cad.decompose([p.copy()], fixed=[(1, 2)])
    unbound._data[].fixed.clear()
    assert_false(cad.check(unbound))
    var witness = cad.decompose([p.copy()], fixed=[(1, 2)])
    _ = witness._data[].inputs.pop()
    assert_false(cad.check(witness))
    var masks = cad.cover([cad.zero(p)], fixed=[(1, 2)])
    masks._data[].masks[len(masks._data[].masks)-1] = 7
    assert_false(cad.check(masks))
    # Opposite conjugates have the same defining polynomial; the exact binding matters.
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var conjugate = cad.decompose([p.copy()], fixed=[(1, a.copy())])
    conjugate._data[].fixed[0][1] = algebraic.negative(a)
    assert_false(cad.check(conjugate))
    var invalid = a.witness().copy()
    invalid.high = invalid.low.copy()
    var rejected = False
    try:
        _ = cad.decompose([p.copy()], fixed=[(1, invalid.copy())])
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    for mode in range(5):
        rejected = False
        try:
            if mode == 0:
                _ = cad.decompose([p.copy()], fixed=[(1, 2), (1, 2)])
            elif mode == 1:
                _ = cad.decompose([p.copy()], fixed=[(-1, 2)])
            elif mode == 2:
                _ = cad.decompose([p.copy()], fixed=[(2, 2)])
            elif mode == 3:
                _ = cad.decompose([p.copy()], fixed=[(1, 2)], axes=[0, 1])
            else:
                _ = cad.checker([p.copy()], fixed=[(1, 2)], axes=[0, 1])
        except e:
            rejected = e.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    assert_true(cad.check(original))
    print("fixed CAD evidence checks passed")
