"""Restricted streams require the requested axes, exact values, and completion."""
from poly import cad, dense, algebraic
from polynomials import Budget, PolyError, PolyErrorKind, CADSink, CADProjectionCertificate, CADStreamEvent
from std.testing import assert_true, assert_false, assert_equal


@fieldwise_init
struct Tape(CADSink, Movable):
    var projection: Optional[CADProjectionCertificate]
    var events: List[CADStreamEvent]

    def begin(mut self, projection: CADProjectionCertificate, mut budget: Budget) raises PolyError:
        self.projection = projection.copy()

    def write(mut self, var event: CADStreamEvent, mut budget: Budget) raises PolyError:
        self.events.append(event^)


def main() raises:
    var p = [((1, 0), 1), ((0, 1), -1)]
    var budget = Budget()
    var tape = Tape(None, List[CADStreamEvent]())
    var receipt = cad.decompose_to([p.copy()], tape, fixed=[(1, 2)])
    assert_equal(receipt.cells, 3)
    assert_equal(receipt.stacks, 2)
    var checker = cad.checker([p.copy()], fixed=[(1, 2)])
    checker.begin(tape.projection.value(), budget)
    for event in tape.events:
        assert_true(checker.accept(event.copy(), budget))
    assert_true(checker.finish(receipt, budget))
    assert_false(checker.accept(tape.events[0].copy(), budget))
    assert_false(checker.finish(receipt, budget))
    var truncated = cad.checker([p.copy()], fixed=[(1, 2)])
    truncated.begin(tape.projection.value(), budget)
    for i in range(len(tape.events)-1):
        assert_true(truncated.accept(tape.events[i].copy(), budget))
    assert_false(truncated.finish(receipt, budget))
    var reordered = cad.checker([p.copy()], fixed=[(1, 2)])
    reordered.begin(tape.projection.value(), budget)
    assert_false(reordered.accept(tape.events[1].copy(), budget))
    assert_false(reordered.finish(receipt, budget))
    var exhausted = cad.checker([p.copy()], fixed=[(1, 2)])
    exhausted.begin(tape.projection.value(), budget)
    var tiny = Budget(0)
    var rejected = False
    try:
        _ = exhausted.accept(tape.events[0].copy(), tiny)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    assert_false(exhausted.accept(tape.events[0].copy(), budget))
    assert_false(exhausted.finish(receipt, budget))
    var axis_checker = cad.checker([p.copy()], fixed=[(0, 2)])
    rejected = False
    try:
        axis_checker.begin(tape.projection.value(), budget)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    assert_false(axis_checker.finish(receipt, budget))
    # Equal defining polynomials cannot hide a different selected conjugate.
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var conjugate_tape = Tape(None, List[CADStreamEvent]())
    var conjugate_receipt = cad.decompose_to([p.copy()], conjugate_tape, fixed=[(1, a.copy())])
    var wrong = cad.checker([p.copy()], fixed=[(1, algebraic.negative(a))])
    wrong.begin(conjugate_tape.projection.value(), budget)
    var accepted = True
    for event in conjugate_tape.events:
        accepted = wrong.accept(event.copy(), budget) and accepted
    assert_false(accepted)
    assert_false(wrong.finish(conjugate_receipt, budget))
    # All-fixed streams and mixed tuples work directly with the checker as sink.
    var all_fixed = cad.checker([p.copy()], fixed=((1, a.copy()), (0, 2)), axes=[1, 0])
    var all_receipt = cad.decompose_to([p.copy()], all_fixed, budget, fixed=((0, 2), (1, a.copy())), axes=[1, 0])
    assert_true(all_fixed.finish(all_receipt, budget))
    assert_equal(all_receipt.cells, 1)
    print("fixed CAD stream checks passed")
