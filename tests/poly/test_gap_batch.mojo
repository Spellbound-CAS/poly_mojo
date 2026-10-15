"""Batch shapes, broadcasting and reusable budget-aware functions."""
from std.testing import assert_equal, assert_true
from polynomials import Budget, PolyError
from poly import dense, batch


def plus(a: dense.Array[0], b: dense.Array[0], mut budget: Budget) raises PolyError -> dense.Array[0]:
    return dense.add(a, b, budget)


def slope(a: dense.Array[0], mut budget: Budget) raises PolyError -> dense.Array[0]:
    return dense.derivative(a, budget)


def main() raises:
    var p = dense.asarray([1, 2])
    var q = dense.asarray([3])
    var a = batch.asarray([p.copy(), q.copy()], shape=List[Int]([2, 1]))
    var b = batch.asarray([p.copy(), q.copy(), p.copy()], shape=List[Int]([1, 3]))
    var grid = batch.vmap(plus, a, b)
    assert_equal(grid.shape()[0], 2)
    assert_equal(grid.shape()[1], 3)
    assert_equal(grid[4][0], 6)
    var shared = batch.asarray(grid)
    assert_equal(shared.ndim(), 2)
    assert_equal(shared.shape()[1], 3)
    var d = batch.vmap(slope, grid)
    assert_equal(len(d[4]), 0)
    var sum = batch.reduce(plus, a, dense.asarray(List[Int]()))
    assert_equal(sum[0], 4)
    var scan = batch.accumulate(plus, a, dense.asarray(List[Int]()))
    assert_equal(scan[1][0], 4)
    var product = batch.outer(plus, a, b)
    assert_equal(product.ndim(), 4)
    assert_equal(len(product), 6)
    var scalar = batch.scalar(p)
    assert_equal(len(batch.broadcast_to(scalar, List[Int]([0, 3]))), 0)
    assert_equal(len(batch.vmap(plus, batch.asarray(List[dense.Array[0]]()), scalar)), 0)
    print("array gap batch checks passed")
