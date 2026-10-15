"""Fiber streams reject truncation, reordering, conjugate headers, and failed reuse."""
from poly import cad, dense, algebraic
from polynomials import Budget, PolyError, PolyErrorKind, CADSink, CADProjectionCertificate, CADStreamEvent, CADStreamReceipt
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
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    var tape = Tape(None, List[CADStreamEvent]())
    var budget = Budget()
    var receipt = cad.decompose_to([p.copy()], tape, budget)
    assert_equal(receipt.stacks, 2)
    assert_equal(receipt.cells, 3)
    assert_equal(len(tape.events), 5)
    var complete = cad.checker([p.copy()])
    complete.begin(tape.projection.value(), budget)
    for event in tape.events:
        assert_true(complete.accept(event.copy(), budget))
    assert_true(complete.finish(receipt, budget))
    assert_false(complete.accept(tape.events[4].copy(), budget))
    assert_false(complete.finish(receipt, budget))
    var truncated = cad.checker([p.copy()])
    truncated.begin(tape.projection.value(), budget)
    for i in range(4):
        assert_true(truncated.accept(tape.events[i].copy(), budget))
    assert_false(truncated.finish(receipt, budget))
    var reordered = cad.checker([p.copy()])
    reordered.begin(tape.projection.value(), budget)
    assert_false(reordered.accept(tape.events[1].copy(), budget))
    assert_false(reordered.accept(tape.events[0].copy(), budget))
    var changed = tape.events[0].copy()
    _ = changed.stack.value().samples.pop()
    var corrupt = cad.checker([p.copy()])
    corrupt.begin(tape.projection.value(), budget)
    assert_false(corrupt.accept(changed^, budget))
    var exhausted = cad.checker([p.copy()])
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
    # An all-fixed, zero-user-axis stream still has exactly one terminal event.
    var constant = [((), a.copy())]
    var scalar = cad.checker([constant.copy()])
    var scalar_receipt = cad.decompose_to([constant.copy()], scalar, budget)
    assert_true(scalar.finish(scalar_receipt, budget))
    assert_equal(scalar_receipt.stacks, 1)
    assert_equal(scalar_receipt.cells, 1)
    print("array CAD fiber stream checks passed")
