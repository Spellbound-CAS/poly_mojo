# expect: axis must be None, Int, Optional[Int], or Int axes
from poly import batch
from polynomials import Budget, PolyError

def plus(a: Int, b: Int, mut budget: Budget) raises PolyError -> Int:
    budget.charge(1)
    return a + b

def main() raises:
    _ = batch.reduce(plus, batch.asarray([1, 2]), 0, axis=Float64(0))
