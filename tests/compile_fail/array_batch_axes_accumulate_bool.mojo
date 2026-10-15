# expect: accumulation axis must be None, Int, or Optional[Int]
from poly import batch
from polynomials import Budget, PolyError
def plus(a: Int, b: Int, mut budget: Budget) raises PolyError -> Int:
    return a + b
def main() raises:
    _ = batch.accumulate(plus, batch.asarray([1, 2]), 0, axis=True)
