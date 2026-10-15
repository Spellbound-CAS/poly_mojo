"""Work is charged before state advances; retries preserve complete evidence."""
from poly import dense, algebraic
from apn_mojo import BallContext
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = [a.copy()]
    var bc = BallContext(precision=64)
    var budget = Budget()
    _ = dense.enclose_stream(p, 0, budget, context=bc)
    var spent = UInt64.MAX - budget.remaining()
    var exact = Budget(spent)
    var stream = dense.enclose_stream(p, 0, exact, context=bc)
    assert_equal(exact.remaining(), UInt64(0))
    var rejected = False
    try:
        _ = dense.enclose_stream(p, 0, context=bc, max_steps=spent - 1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    budget = Budget()
    var pending_1 = stream.next(budget)
    _ = pending_1.take()
    spent = UInt64.MAX - budget.remaining()
    stream = dense.enclose_stream(p, 0, context=bc)
    var short = Budget(spent - 1)
    rejected = False
    try:
        _ = stream.next(short)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    exact = Budget(spent)
    var pending_2 = stream.next(exact)
    var first = pending_2.take()
    assert_equal(exact.remaining(), UInt64(0))
    assert_equal(first.kind, 0)
    assert_equal(first.index, 0)
    budget = Budget()
    _ = dense.enclose_checker(p, 0, budget, context=bc)
    spent = UInt64.MAX - budget.remaining()
    exact = Budget(spent)
    var checker = dense.enclose_checker(p, 0, exact, context=bc)
    assert_equal(exact.remaining(), UInt64(0))
    rejected = False
    try:
        _ = dense.enclose_checker(p, 0, context=bc, max_steps=spent - 1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    budget = Budget()
    assert_true(checker.accept(first.copy(), budget))
    spent = UInt64.MAX - budget.remaining()
    checker = dense.enclose_checker(p, 0, context=bc)
    short = Budget(spent - 1)
    rejected = False
    try:
        _ = checker.accept(first.copy(), short)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    exact = Budget(spent)
    assert_true(checker.accept(first.copy(), exact))
    assert_equal(exact.remaining(), UInt64(0))
    # The first arithmetic call rebuilds the graph from checked inputs. A failed
    # rebuild must leave the event available on both producer and checker.
    short = Budget(0)
    rejected = False
    try:
        _ = stream.next(short)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    budget = Budget()
    var pending_3 = stream.next(budget)
    var step = pending_3.take()
    short = Budget(0)
    rejected = False
    try:
        _ = checker.accept(step.copy(), short)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    assert_true(checker.accept(step^, budget))
    while True:
        var next = stream.next(budget)
        if not next:
            break
        assert_true(checker.accept(next.take(), budget))
    var receipt = stream.finish(budget)
    short = Budget(0)
    rejected = False
    try:
        _ = stream.finish(short)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    assert_true(checker.finish(receipt, budget))
    short = Budget(0)
    rejected = False
    try:
        _ = checker.finish(receipt, short)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    assert_true(checker.finish(receipt, budget))
    var missing = BallContext()
    rejected = False
    try:
        _ = dense.enclose_stream(p, 0, context=missing)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    assert_true(algebraic.check(a))
    assert_equal(p[0], a)
    print("array algebraic refinement budget checks passed")
