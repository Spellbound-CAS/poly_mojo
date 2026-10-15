"""Streaming completion, augmented axis order, budgets, and witness rejection."""
from poly import cad, dense, algebraic
from polynomials import Budget, PolyErrorKind
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    var budget = Budget()
    var checker = cad.checker([p.copy()], budget)
    var receipt = cad.decompose_to([p.copy()], checker, budget)
    assert_true(checker.finish(receipt, budget))
    receipt.cells += 1
    assert_false(checker.finish(receipt, budget))
    var projection = cad.project([p.copy()])
    var wrong = projection.evidence().copy()
    var first_axis = wrong.variables[0]
    wrong.variables[0] = wrong.variables[1]
    wrong.variables[1] = first_axis
    var rejected = False
    var bad_checker = cad.checker([p.copy()])
    try:
        bad_checker.begin(wrong, budget)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    assert_false(bad_checker.finish(receipt, budget))
    # Budget accounting includes preparation and the retained embedding.
    var small = [(0, algebraic.asreal(-1)), (1, algebraic.asreal(1))]
    _ = cad.project([small.copy()])
    budget = Budget()
    var result = cad.project([small.copy()], budget)
    var used = UInt64.MAX - budget.remaining()
    assert_true(cad.check(cad.project([small.copy()], max_steps=used)))
    rejected = False
    try:
        _ = cad.project([small.copy()], max_steps=used-1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    budget = Budget()
    assert_true(cad.check(result, budget))
    used = UInt64.MAX - budget.remaining()
    assert_true(cad.check(result, max_steps=used))
    rejected = False
    try:
        _ = cad.check(result, max_steps=used-1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    var corrupt = dense.root_between([-3, 0, 1], 1, 2)
    corrupt._data[].plain.high = corrupt._data[].plain.low.copy()
    rejected = False
    var invalid = [(0, corrupt.copy())]
    try:
        _ = cad.project([invalid.copy()])
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    # Retained inputs and isolating constraints are part of the checked claim.
    var changed = cad.project([p.copy()])
    _ = changed._data[].inputs.pop()
    assert_false(cad.check(changed))
    var weakened = cad.solve([cad.zero(p)])
    weakened._data[].masks[len(weakened._data[].masks)-1] = 7
    assert_false(cad.check(weakened))
    print("array algebraic CAD stream checks passed")
