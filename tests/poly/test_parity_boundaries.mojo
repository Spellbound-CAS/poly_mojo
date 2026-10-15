"""Shared allowances, single-pass input, exact operand binding and CAD failures."""
from std.testing import assert_true, assert_false, assert_equal
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from apn_mojo import Rational, BallContext
from polynomials import Budget, PolyError, PolyErrorKind
from poly import dense, sparse, exact, cad, algebraic


@fieldwise_init
struct Count(Copyable, Movable):
    var calls: Int
    var position: Int


@fieldwise_init
struct Once(Iterable):
    var state: ArcPointer[Count]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = Cursor

    def __iter__(self) -> Cursor:
        return Cursor(self.state.copy())


@fieldwise_init
struct Cursor(Iterator):
    comptime Element = Int
    var state: ArcPointer[Count]

    def __next__(mut self) raises StopIteration -> Int:
        self.state[].calls += 1
        if self.state[].position == 3:
            raise StopIteration()
        self.state[].position += 1
        return self.state[].position


def workflow(context: BallContext, mut budget: Budget) raises PolyError:
    var division = dense.divide_many_checked([1, 0, 1], List[List[Int]](), budget)
    _ = exact.check(division, budget)
    _ = exact.remainder(division, budget)
    var checked = algebraic.compare_checked(1, 2, budget)
    _ = algebraic.check(checked, budget)
    var inputs = [[(0, -1), (1, 1)]]
    var projection = cad.project(inputs, budget)
    _ = cad.check(projection, budget)
    _ = cad.projection_polynomials(projection, budget)
    var cell = cad.cell_at(inputs, [0], budget)
    _ = cad.check(cell, budget)
    _ = cad.sample(cell, budget)
    _ = cad.signs(cell, budget)
    var stream = dense.enclose_stream([1, 2, 3], 2, budget, context=context)
    var checker = dense.enclose_checker([1, 2, 3], 2, budget, context=context)
    while True:
        var step = stream.next(budget)
        if not step:
            break
        _ = checker.accept(step.take(), budget)
    _ = checker.finish(stream.finish(budget), budget)


def main() raises:
    var context = BallContext(precision=32)
    var measured = Budget()
    workflow(context, measured)
    var cost = UInt64.MAX - measured.remaining()
    var enough = Budget(cost)
    workflow(context, enough)
    assert_equal(enough.remaining(), UInt64(0))
    var short = Budget(cost - 1)
    var failed = False
    try:
        workflow(context, short)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var once = Once(ArcPointer(Count(0, 0)))
    var result = dense.gcd_checked(once, [1])
    assert_true(exact.check(result))
    assert_equal(once.state[].calls, 4)
    var untouched = Once(ArcPointer(Count(0, 0)))
    failed = False
    try:
        _ = dense.gcd_checked(untouched, [1], max_steps=0)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    assert_equal(untouched.state[].calls, 0)
    var preserved = dense.asarray([1, 2, 3])
    failed = False
    try:
        _ = dense.divide_checked(preserved, [0])
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    assert_equal(preserved.to_list(), dense.asarray([1, 2, 3]).to_list())
    var bad = dense.divide_checked([1, 0, 1], [1, 2])
    bad._data[].others.clear()
    assert_false(exact.check(bad))
    var many = sparse.divide_many_checked([(0, 1), (2, 1)], List[List[Tuple[Int, Int]]]())
    many._data[].proof.remainder.terms.clear()
    assert_false(exact.check(many))
    var inputs = [[((1, 0), 1)], [((0, 1), 1)]]
    failed = False
    try:
        _ = cad.project(inputs, axes=[0, 0])
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    failed = False
    try:
        _ = cad.cell_at(inputs, [0, 0, 0])
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var cell = cad.cell_at(inputs, (0, 0), axes=[1, 0])
    cell._data[].variables[0] = 0
    assert_false(cad.check(cell))
    var full = cad.cover([cad.positive([(1, 1)])])
    full._data[].proof.witness_cell = -1
    assert_false(cad.check(full))
    print("array parity boundary checks passed")
