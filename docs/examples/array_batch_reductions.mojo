"""Reduce rows of polynomials while keeping their batch dimensions."""
from poly import batch, dense
from polynomials import Budget, PolyError
from std.testing import assert_equal


def add_polynomials(a: dense.Array[0], b: dense.Array[0], mut budget: Budget) raises PolyError -> dense.Array[0]:
    return dense.add(a, b, budget)


def main() raises:
    var grid = batch.asarray([
        dense.asarray([1, 1]), dense.asarray([2]), dense.asarray([3, -1]),
        dense.asarray([4]), dense.asarray([5, 1]), dense.asarray([6]),
    ], shape=List[Int]([2, 3]))
    var zero = dense.asarray(List[Int]())
    var rows = batch.reduce(add_polynomials, grid, zero, axis=1, keepdims=True)
    assert_equal(rows.shape(), List[Int]([2, 1]))
    assert_equal(dense.evaluate(rows[0], 2), 6)
    assert_equal(dense.evaluate(rows[1], 2), 17)
    print("Row sums at x=2:", dense.evaluate(rows[0], 2), dense.evaluate(rows[1], 2))
    var total = batch.reduce(add_polynomials, grid, zero, axis=None)
    assert_equal(total.ndim(), 0)
    assert_equal(dense.evaluate(total.item(), 2), 23)
    print("Total at x=2:", dense.evaluate(total.item(), 2))
