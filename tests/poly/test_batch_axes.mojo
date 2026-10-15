"""Batch-axis reductions preserve shape, fold order, and shared budgets."""
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from std.testing import assert_equal, assert_true
from poly import batch
from polynomials import Budget, PolyError, PolyErrorKind


def plus(a: Int, b: Int, mut budget: Budget) raises PolyError -> Int:
    budget.charge(1)
    return a + b


def digits(a: Int, b: Int, mut budget: Budget) raises PolyError -> Int:
    budget.charge(1)
    return 10 * a + b


def fail_on_five(a: Int, b: Int, mut budget: Budget) raises PolyError -> Int:
    budget.charge(1)
    if b == 5:
        raise PolyError(PolyErrorKind.DOMAIN, "five is rejected")
    return a + b



@fieldwise_init
struct AxisState(Copyable, Movable):
    var calls: Int
    var position: Int


@fieldwise_init
struct AxisSource(Iterable):
    var state: ArcPointer[AxisState]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = AxisCursor

    def __iter__(self) -> AxisCursor:
        return AxisCursor(self.state.copy())


@fieldwise_init
struct AxisCursor(Iterator):
    comptime Element = Int
    var state: ArcPointer[AxisState]

    def __next__(mut self) raises StopIteration -> Int:
        self.state[].calls += 1
        if self.state[].position == 2:
            raise StopIteration()
        var axis = self.state[].position
        self.state[].position += 1
        return axis


def main() raises:
    var grid = batch.asarray([1, 2, 3, 4, 5, 6], shape=List[Int]([2, 3]))
    assert_equal(batch.reduce(plus, grid, 0), 21)
    var scan_rows = batch.accumulate(digits, grid, 0, axis=-1)
    assert_equal(scan_rows.shape(), grid.shape())
    assert_equal(List(scan_rows), List[Int]([1, 12, 123, 4, 45, 456]))
    assert_equal(List(batch.accumulate(digits, grid, 0, axis=0)), List[Int]([1, 2, 3, 14, 25, 36]))
    assert_equal(List(batch.accumulate(plus, grid, 10, axis=1)), List[Int]([11, 13, 16, 14, 19, 25]))
    assert_equal(List(batch.accumulate(digits, grid, 0, axis=None)), List[Int]([1, 12, 123, 1234, 12345, 123456]))
    assert_equal(List(batch.accumulate(plus, grid, 0, axis=Optional[Int](1))), List[Int]([1, 3, 6, 4, 9, 15]))
    assert_equal(batch.accumulate(plus, batch.scalar(2), 3, axis=None).item(), 5)
    for invalid in [2, -3, Int.MIN, Int.MAX]:
        var bad = False
        try:
            _ = batch.accumulate(plus, grid, 0, axis=invalid)
        except error:
            bad = error.kind == PolyErrorKind.DOMAIN
        assert_true(bad)
    var scan_budget = Budget()
    _ = batch.accumulate(plus, grid, 0, scan_budget, axis=0)
    var scan_spent = UInt64.MAX - scan_budget.remaining()
    var scan_exact = Budget(scan_spent)
    _ = batch.accumulate(plus, grid, 0, scan_exact, axis=0)
    assert_equal(scan_exact.remaining(), UInt64(0))
    var scan_limited = False
    try:
        _ = batch.accumulate(plus, grid, 0, axis=0, max_steps=scan_spent - 1)
    except error:
        scan_limited = error.kind == PolyErrorKind.LIMIT
    assert_true(scan_limited)
    var scan_failed = False
    try:
        _ = batch.accumulate(fail_on_five, grid, 0, axis=0)
    except error:
        scan_failed = error.kind == PolyErrorKind.DOMAIN
    assert_true(scan_failed)
    assert_equal(grid[4], 5)
    var rows = batch.reduce(plus, grid, 10, axis=1)
    assert_equal(rows.shape(), List[Int]([2]))
    assert_equal(rows[0], 16)
    assert_equal(rows[1], 25)
    var columns = batch.reduce(plus, grid, 0, axis=0, keepdims=True)
    assert_equal(columns.shape(), List[Int]([1, 3]))
    assert_equal(columns[0], 5)
    assert_equal(columns[1], 7)
    assert_equal(columns[2], 9)
    assert_equal(batch.reduce(digits, grid, 0, axis=-1)[1], 456)
    assert_equal(batch.reduce(digits, grid, 0, axis=0)[1], 25)
    assert_equal(batch.reduce(digits, grid, 0, axis=(1, 0)).item(), 123456)
    assert_equal(batch.reduce(digits, grid, 0, axis=(0, 1)).item(), 123456)
    assert_equal(batch.reduce(plus, grid, 0, axis=range(2)).item(), 21)
    var once = AxisSource(ArcPointer(AxisState(0, 0)))
    assert_equal(batch.reduce(plus, grid, 0, axis=once).item(), 21)
    assert_equal(once.state[].calls, 3)
    var blocked = AxisSource(ArcPointer(AxisState(0, 0)))
    var early_limit = False
    try:
        _ = batch.reduce(plus, grid, 0, axis=blocked, max_steps=0)
    except error:
        early_limit = error.kind == PolyErrorKind.LIMIT
    assert_true(early_limit)
    assert_equal(blocked.state[].calls, 0)
    var total = batch.reduce(plus, grid, 10, axis=None)
    assert_equal(total.ndim(), 0)
    assert_equal(total.item(), 31)
    assert_equal(batch.reduce(plus, grid, 0, keepdims=False).ndim(), 0)
    assert_equal(batch.reduce(plus, grid, 0, keepdims=True).shape(), List[Int]([1, 1]))
    var option = Optional[Int](-1)
    assert_equal(batch.reduce(plus, grid, 0, axis=option)[0], 6)
    option = None
    assert_equal(batch.reduce(plus, grid, 0, axis=option).item(), 21)
    var untouched_axes = batch.reduce(plus, grid, 10, axis=(), keepdims=True)
    assert_equal(untouched_axes.shape(), grid.shape())
    assert_equal(untouched_axes[4], 15)
    var empty_axes = batch.reduce(plus, grid, 10, axis=List[Int]())
    assert_equal(empty_axes[5], 16)
    var cube = batch.asarray([1, 2, 3, 4, 5, 6, 7, 8], shape=List[Int]([2, 2, 2]))
    var slices = batch.reduce(digits, cube, 0, axis=[2, 0])
    assert_equal(slices.shape(), List[Int]([2]))
    assert_equal(slices[0], 1256)
    assert_equal(slices[1], 3478)
    assert_equal(batch.reduce(plus, cube, 0, axis=(-1, 0), keepdims=True).shape(), List[Int]([1, 2, 1]))
    var repeated = batch.broadcast_to(batch.asarray([1, 2, 3], shape=List[Int]([1, 3])), List[Int]([2, 3]))
    assert_equal(batch.reduce(digits, repeated, 0, axis=0)[2], 33)
    assert_equal(batch.reduce(digits, repeated, 0, axis=None).item(), 123123)
    assert_equal(List(batch.accumulate(digits, repeated, 0, axis=0)), List[Int]([1, 2, 3, 11, 22, 33]))
    var scalar = batch.scalar(4)
    assert_equal(batch.reduce(plus, scalar, 1, axis=None, keepdims=True).item(), 5)
    assert_equal(batch.reduce(plus, scalar, 1, axis=()).ndim(), 0)
    var empty = batch.asarray(List[Int](), shape=List[Int]([2, 0, 3]))
    assert_equal(len(batch.accumulate(plus, empty, 7, axis=1)), 0)
    var identities = batch.reduce(plus, empty, 7, axis=1)
    assert_equal(identities.shape(), List[Int]([2, 3]))
    assert_equal(len(identities), 6)
    assert_equal(identities[5], 7)
    assert_equal(len(batch.reduce(plus, empty, 7, axis=0)), 0)
    assert_equal(batch.reduce(plus, empty, 7, axis=None).item(), 7)
    var huge_empty = batch.asarray(List[Int](), shape=List[Int]([0, Int.MAX, 0, Int.MAX]))
    assert_equal(len(batch.reduce(plus, huge_empty, 0, axis=2)), 0)
    assert_equal(len(batch.reduce(plus, huge_empty, 0, axis=(1, 3))), 0)
    for invalid in [List[Int]([0, -2]), List[Int]([2]), List[Int]([-3]), List[Int]([Int.MIN]), List[Int]([Int.MAX])]:
        var rejected = False
        try:
            _ = batch.reduce(plus, grid, 0, axis=invalid)
        except error:
            rejected = error.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    var rejected = False
    try:
        _ = batch.reduce(plus, scalar, 0, axis=0)
    except error:
        rejected = error.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    rejected = False
    try:
        _ = batch.reduce(plus, huge_empty, 0, axis=(0, 2))
    except error:
        rejected = error.kind == PolyErrorKind.CAPACITY
    assert_true(rejected)
    var work = Budget()
    _ = batch.reduce(plus, repeated, 0, work, axis=-1, keepdims=True)
    var spent = UInt64.MAX - work.remaining()
    var exact = Budget(spent)
    _ = batch.reduce(plus, repeated, 0, exact, axis=-1, keepdims=True)
    assert_equal(exact.remaining(), UInt64(0))
    rejected = False
    try:
        _ = batch.reduce(plus, repeated, 0, axis=-1, keepdims=True, max_steps=spent - 1)
    except error:
        rejected = error.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    rejected = False
    try:
        _ = batch.reduce(fail_on_five, grid, 0, axis=0)
    except error:
        rejected = error.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    assert_equal(grid[4], 5)
    rejected = False
    try:
        _ = grid.item()
    except error:
        rejected = error.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    rejected = False
    try:
        _ = empty.item()
    except error:
        rejected = error.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    var item_budget = Budget(1)
    assert_equal(total.item(item_budget), 31)
    assert_equal(item_budget.remaining(), UInt64(0))
    print("batch axis checks passed")
