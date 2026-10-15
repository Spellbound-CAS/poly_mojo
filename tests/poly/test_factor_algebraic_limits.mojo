"""Shared budgets cover preparation, solving, checking and output materialization."""
from poly import dense, sparse, algebraic, coefficients, exact
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var p = dense.asarray(coefficients((algebraic.asreal(-2), 0, 1)))
    var before = p.to_list()
    for steps in [0, 1, 100, 1000]:
        var limited = False
        var budget = Budget(UInt64(steps))
        try:
            _ = dense.factor(p, budget)
        except e:
            limited = e.kind == PolyErrorKind.LIMIT
        assert_true(limited)
        assert_equal(p.to_list(), before)
    # Prime discovery is cached by the core. Measure a repeatable warm call.
    _ = dense.factor(p)
    var measured = Budget()
    _ = dense.factor(p, measured)
    var used = UInt64.MAX - measured.remaining()
    var exact_budget = Budget(used)
    _ = dense.factor(p, exact_budget)
    assert_equal(exact_budget.remaining(), UInt64(0))
    var short_budget = Budget(used - 1)
    var limited = False
    try:
        _ = dense.factor(p, short_budget)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    measured = Budget()
    var retained = dense.factor_checked(p, measured)
    used = UInt64.MAX - measured.remaining()
    limited = False
    try:
        _ = dense.factor_checked(p, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(exact.check(retained))
    measured = Budget()
    assert_true(exact.check(retained, measured))
    used = UInt64.MAX - measured.remaining()
    assert_true(exact.check(retained, max_steps=used))
    limited = False
    try:
        _ = exact.check(retained, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    var s = sparse.asarray((((0, 2), algebraic.ascomplex(1)), ((0, 0), 1)))
    measured = Budget()
    _ = sparse.factor(s, measured)
    used = UInt64.MAX - measured.remaining()
    exact_budget = Budget(used)
    _ = sparse.factor(s, exact_budget)
    assert_equal(exact_budget.remaining(), UInt64(0))
    limited = False
    try:
        _ = sparse.factor(s, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(len(s), 2)
    # Charge the full dense image before allocating sparse high-degree inputs.
    var large = sparse.asarray(((1000000000, algebraic.asreal(1)), (0, 1)))
    limited = False
    try:
        _ = sparse.factor(large, max_steps=10000)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(len(large), 2)
    print("array algebraic factor budget checks passed")
