"""Preparation, traversal, access, and checking share deterministic work budgets."""
from poly import cad
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_true, assert_equal


def exercise(mode: Int, mut budget: Budget) raises PolyError:
    var p = [(0, -1), (1, 1)]
    if mode == 0:
        _ = cad.decompose([p.copy()], budget, fixed=[(0, 1)])
    elif mode == 1:
        _ = cad.solve([cad.zero(p, budget)], budget, fixed=[(0, 1)])
    elif mode == 2:
        _ = cad.cover([cad.zero(p, budget)], budget, fixed=[(0, 1)])
    else:
        var checker = cad.checker([p.copy()], budget, fixed=[(0, 1)])
        var receipt = cad.decompose_to([p.copy()], checker, budget, fixed=[(0, 1)])
        if not checker.finish(receipt, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "incomplete fixed-coordinate stream")


def main() raises:
    for mode in range(4):
        var budget = Budget()
        exercise(mode, budget)
        budget = Budget()
        exercise(mode, budget)
        var used = UInt64.MAX - budget.remaining()
        budget = Budget(used)
        exercise(mode, budget)
        assert_equal(budget.remaining(), UInt64(0))
        budget = Budget(used-1)
        var rejected = False
        try:
            exercise(mode, budget)
        except e:
            rejected = e.kind == PolyErrorKind.LIMIT
        assert_true(rejected)
    var result = cad.decompose([[(0, -1), (1, 1)]], fixed=[(0, 1)])
    var budget = Budget()
    assert_true(cad.check(result, budget))
    var used = UInt64.MAX - budget.remaining()
    assert_true(cad.check(result, max_steps=used))
    var rejected = False
    try:
        _ = cad.check(result, max_steps=used-1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    rejected = False
    try:
        _ = cad.fixed_coordinates(result, max_steps=0)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    assert_true(cad.check(result))
    print("fixed CAD budget checks passed")
