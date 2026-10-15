"""Shared work boundaries for retained, selective, and streamed fiber CAD."""
from poly import cad, dense, algebraic
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_true, assert_equal


def exercise(mode: Int, p: List[Tuple[Int, algebraic.Real]], mut budget: Budget) raises PolyError:
    if mode == 0:
        _ = cad.decompose([p.copy()], budget)
    elif mode == 1:
        _ = cad.solve([cad.negative(p, budget), cad.nonnegative(p, budget)], budget)
    elif mode == 2:
        _ = cad.cover([cad.zero(p, budget)], budget)
    else:
        var checker = cad.checker([p.copy()], budget)
        var receipt = cad.decompose_to([p.copy()], checker, budget)
        if not checker.finish(receipt, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "incomplete fiber stream")


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p: List[Tuple[Int, algebraic.Real]] = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    for mode in range(4):
        var budget = Budget()
        exercise(mode, p, budget)
        budget = Budget()
        exercise(mode, p, budget)
        var used = UInt64.MAX - budget.remaining()
        budget = Budget(used)
        exercise(mode, p, budget)
        assert_equal(budget.remaining(), UInt64(0))
        budget = Budget(used - 1)
        var rejected = False
        try:
            exercise(mode, p, budget)
        except e:
            rejected = e.kind == PolyErrorKind.LIMIT
        assert_true(rejected)
        assert_true(p[0][1] == algebraic.negative(a))
    var result = cad.decompose([p.copy()])
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
    assert_true(cad.check(result))
    print("array CAD fiber budget checks passed")
