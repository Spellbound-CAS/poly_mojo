"""Budgets cover sparse preparation, image storage, search and checking."""
from poly import sparse, algebraic, exact
from poly import _algebraic_multivariate as implementation
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var p = sparse.asarray((((1, 0), algebraic.asreal(1)), ((0, 1), 1)))
    var before = p.to_list()
    for allowance in [0, 1, 50, 100]:
        var limited = False
        try:
            _ = sparse.factor(p, max_steps=UInt64(allowance))
        except e:
            limited = e.kind == PolyErrorKind.LIMIT
        assert_true(limited)
        assert_equal(p.to_list(), before)
    _ = sparse.factor(p)
    var budget = Budget()
    _ = sparse.factor(p, budget)
    var used = UInt64.MAX - budget.remaining()
    budget = Budget(used)
    _ = sparse.factor(p, budget)
    assert_equal(budget.remaining(), UInt64(0))
    var limited = False
    try:
        _ = sparse.factor(p, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    budget = Budget()
    var retained = sparse.factor_checked(p, budget)
    used = UInt64.MAX - budget.remaining()
    budget = Budget(used)
    _ = sparse.factor_checked(p, budget)
    assert_equal(budget.remaining(), UInt64(0))
    limited = False
    try:
        _ = sparse.factor_checked(p, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    budget = Budget()
    assert_true(exact.check(retained, budget))
    used = UInt64.MAX - budget.remaining()
    budget = Budget(used)
    assert_true(exact.check(retained, budget))
    assert_equal(budget.remaining(), UInt64(0))
    limited = False
    try:
        _ = exact.check(retained, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    # Coordinate factors need no degree-sized buffer even at large exponents.
    var monomial = sparse.asarray((((1000000000, 2000000000), algebraic.ascomplex(2, 1)),))
    var mr = sparse.factor_checked(monomial, max_steps=100000000)
    var mf = exact.factors(mr)
    assert_equal(len(mf), 2)
    assert_equal(mf[0][1], 2000000000)
    assert_equal(mf[1][1], 1000000000)
    assert_true(exact.check(mr, max_steps=100000000))
    # A dense degree box is charged before allocation; overflow is distinct.
    var large = sparse.asarray((((1000000000, 0), algebraic.asreal(1)), ((0, 1000000000), 1), ((0, 0), 1)))
    var core = sparse.to_core(large)
    budget = Budget()
    var box = implementation._box(core, budget)
    limited = False
    budget = Budget(10000)
    try:
        _ = implementation._encode(core, box[0], box[1], box[2], budget)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT or e.kind == PolyErrorKind.CAPACITY
    assert_true(limited)
    # Three representable variable degrees can overflow the dense extent.
    var overflow = sparse.asarray((((2000000000, 0, 0), algebraic.asreal(1)), ((0, 2000000000, 0), 1), ((0, 0, 2000000000), 1)))
    budget = Budget()
    var capacity = False
    try:
        _ = implementation._box(sparse.to_core(overflow), budget)
    except e:
        capacity = e.kind == PolyErrorKind.CAPACITY
    assert_true(capacity)
    assert_equal(len(large), 3)
    assert_equal(p.to_list(), before)
    print("array multivariate algebraic budget checks passed")
