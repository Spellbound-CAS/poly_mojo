"""Exact complex streams reconstruct independently and reject damaged events."""
from poly import dense, sparse, roots, algebraic
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_equal, assert_true, assert_false

@fieldwise_init
struct Counter(roots.ExactComplexRootSink):
    var items: Int
    var degree: Int
    def write(mut self, var event: roots.ExactComplexRootEvent, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.items += 1
        self.degree += event.multiplicity

struct Reject(roots.ExactComplexRootSink):
    def __init__(out self):
        pass
    def write(mut self, var event: roots.ExactComplexRootEvent, mut budget: Budget) raises PolyError:
        raise PolyError(PolyErrorKind.DOMAIN, "consumer rejected event")

def main() raises:
    var p = dense.asarray([1, 0, 2, 0, 1])
    var budget = Budget()
    var checker = dense.complex_roots_checker(p, budget)
    var receipt = dense.complex_roots_to(p, checker, budget)
    assert_equal(receipt.items, 2)
    assert_equal(receipt.degree, 4)
    assert_true(checker.finish(receipt, budget))
    receipt.items += 1
    assert_false(checker.finish(receipt, budget))
    var counter = Counter(0, 0)
    var sparse_receipt = sparse.complex_roots_to([((0, 0), 1), ((0, 2), 1)], counter, axis=1)
    assert_equal(counter.items, 2)
    assert_equal(counter.degree, 2)
    var sparse_checker = sparse.complex_roots_checker([((0, 0), 1), ((0, 2), 1)], axis=1)
    var i = algebraic.ascomplex(0, 1)
    assert_true(sparse_checker.accept(roots.ExactComplexRootEvent(i.copy(), 1), budget))
    assert_true(sparse_checker.accept(roots.ExactComplexRootEvent(algebraic.negative(i), 1), budget))
    assert_true(sparse_checker.finish(sparse_receipt, budget))
    var duplicate = dense.complex_roots_checker([1, 0, 2, 0, 1])
    assert_true(duplicate.accept(roots.ExactComplexRootEvent(i.copy(), 2), budget))
    assert_false(duplicate.accept(roots.ExactComplexRootEvent(i.copy(), 2), budget))
    assert_false(duplicate.finish(roots.ExactComplexRootReceipt(2, 4), budget))
    var wrong = dense.complex_roots_checker([1, 0, 1])
    assert_false(wrong.accept(roots.ExactComplexRootEvent(algebraic.ascomplex(1), 1), budget))
    var incomplete = dense.complex_roots_checker([1, 0, 1])
    assert_false(incomplete.finish(roots.ExactComplexRootReceipt(0, 0), budget))
    var constant = dense.complex_roots_checker([2])
    var none = dense.complex_roots_to([2], constant)
    assert_true(constant.finish(none, budget))
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var ap = [algebraic.negative(a), algebraic.asreal(1)]
    var ac = dense.complex_roots_checker(ap)
    var ar = dense.complex_roots_to(ap, ac)
    assert_true(ac.finish(ar, budget))
    var zp = [algebraic.negative(algebraic.ascomplex(a, 1)), algebraic.ascomplex(1)]
    var zc = dense.complex_roots_checker(zp)
    var zr = dense.complex_roots_to(zp, zc)
    assert_true(zc.finish(zr, budget))
    var work = Budget()
    var counting = Counter(0, 0)
    _ = dense.complex_roots_to([-1, 1], counting, work)
    var spent = UInt64.MAX - work.remaining()
    var exact = Budget(spent)
    var counting2 = Counter(0, 0)
    _ = dense.complex_roots_to([-1, 1], counting2, exact)
    assert_equal(exact.remaining(), UInt64(0))
    for mode in range(3):
        var rejected = False
        try:
            var sink = Counter(0, 0)
            if mode == 0:
                _ = dense.complex_roots_to([-1, 1], sink, max_steps=spent - 1)
            elif mode == 1:
                var rejecting = Reject()
                _ = dense.complex_roots_to([-1, 1], rejecting)
            else:
                _ = dense.complex_roots_checker([0])
        except error:
            rejected = error.kind == (PolyErrorKind.LIMIT if mode == 0 else PolyErrorKind.DOMAIN)
        assert_true(rejected)
    var poisoned = dense.complex_roots_checker([-1, 1])
    var exhausted = Budget(0)
    try:
        _ = poisoned.accept(roots.ExactComplexRootEvent(algebraic.ascomplex(1), 1), exhausted)
    except error:
        assert_equal(error.kind, PolyErrorKind.LIMIT)
    assert_false(poisoned.finish(roots.ExactComplexRootReceipt(0, 0), budget))
    print("workflow exact stream checks passed")
