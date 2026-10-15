"""Axis reductions keep whole polynomial values and serializable batch shapes."""
from poly import batch, dense, sparse, io
from polynomials import Budget, PolyError
from std.testing import assert_equal


def plus(a: dense.Array[0], b: dense.Array[0], mut budget: Budget) raises PolyError -> dense.Array[0]:
    return dense.add(a, b, budget)


def minus(a: dense.Array[0], b: dense.Array[0], mut budget: Budget) raises PolyError -> dense.Array[0]:
    return dense.subtract(a, b, budget)


def sparse_plus(a: sparse.Array[Tuple[Int, Int], 0], b: sparse.Array[Tuple[Int, Int], 0], mut budget: Budget) raises PolyError -> sparse.Array[Tuple[Int, Int], 0]:
    return sparse.add(a, b, budget)


def main() raises:
    var zero = dense.asarray(List[Int]())
    var grid = batch.asarray([dense.asarray([1, 1]), dense.asarray([1, -1]), dense.asarray([2]), dense.asarray([-2])], shape=List[Int]([2, 2]))
    var scans = batch.accumulate(plus, grid, zero, axis=1)
    assert_equal(scans[0][1], 1)
    assert_equal(scans[1][0], 2)
    assert_equal(len(scans[3]), 0)
    var rows = batch.reduce(plus, grid, zero, axis=1, keepdims=True)
    assert_equal(rows.shape(), List[Int]([2, 1]))
    assert_equal(rows[0][0], 2)
    assert_equal(len(rows[1]), 0)
    var centered = batch.vmap(minus, grid, rows)
    assert_equal(centered.shape(), grid.shape())
    assert_equal(centered[0][0], -1)
    assert_equal(centered[0][1], 1)
    var total = batch.reduce(plus, grid, zero, axis=None)
    assert_equal(total.ndim(), 0)
    assert_equal(total.item()[0], 2)
    var saved = io.loads[type_of(total)](io.dumps(total))
    assert_equal(saved.ndim(), 0)
    assert_equal(saved.item()[0], 2)
    var with_axes = io.loads[type_of(rows)](io.dumps(rows))
    assert_equal(with_axes.shape(), rows.shape())
    assert_equal(len(with_axes[1]), 0)
    var work = Budget()
    var all_kept = batch.reduce(plus, grid, zero, work, keepdims=True)
    assert_equal(all_kept.shape(), List[Int]([1, 1]))
    assert_equal(all_kept.item()[0], 2)
    var p = sparse.asarray([((1, 0), 1), ((0, 1), 2)])
    var q = sparse.asarray([((0, 1), -1)])
    var sparse_zero = sparse.asarray(List[Tuple[Tuple[Int, Int], Int]]())
    var items = batch.asarray([p.copy(), q.copy()], shape=List[Int]([1, 2]))
    var result = batch.reduce(sparse_plus, items, sparse_zero, axis=1)
    assert_equal(result.shape(), List[Int]([1]))
    assert_equal(sparse.evaluate(result.item(), (2, 3)), 5)
    assert_equal(sparse.evaluate(p, (2, 3)), 8)
    print("polynomial batch axis checks passed")
