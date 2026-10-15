"""Accumulate each row while preserving the batch shape."""
from poly import batch
from polynomials import Budget, PolyError
from std.testing import assert_equal

def add(a: Int, b: Int, mut budget: Budget) raises PolyError -> Int:
    budget.charge(1)
    return a + b

def main() raises:
    var values = batch.asarray([1, 2, 3, 4, 5, 6], shape=List[Int]([2, 3]))
    var running = batch.accumulate(add, values, 0, axis=-1)
    assert_equal(List(running), List[Int]([1, 3, 6, 4, 9, 15]))
    assert_equal(running.shape(), values.shape())
    print("Row totals:", running[2], running[5])
