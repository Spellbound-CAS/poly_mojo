"""Array batching."""
from poly import dense, batch
from polynomials import Budget, PolyError


def add_polynomials(a: dense.Array[0], b: dense.Array[0], mut budget: Budget) raises PolyError -> dense.Array[0]:
    return dense.add(a, b, budget)


def main() raises:
    var p = dense.asarray([1, 2])
    var q = dense.asarray([3])
    var rows = batch.asarray([p.copy(), q.copy()], shape=List[Int]([2, 1]))
    var columns = batch.asarray([p.copy(), q.copy(), p.copy()], shape=List[Int]([1, 3]))
    var grid = batch.vmap(add_polynomials, rows, columns)
    print("Batch dimensions:", grid.shape()[0], grid.shape()[1])
    print("Middle item in row two:", grid[4][0])
    var total = batch.reduce(add_polynomials, rows, dense.asarray(List[Int]()))
    print("Sum coefficients:", total[0], total[1])
