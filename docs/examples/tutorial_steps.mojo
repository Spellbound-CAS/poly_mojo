"""Count the steps a calculation takes, then bound it with max_steps."""
import poly
from polynomials import Budget, PolyErrorKind


def main() raises:
    var work = Budget()
    _ = poly.dense.power([1, 1], 100, work)
    print("Steps for (1+x)^100:", UInt64.MAX - work.remaining())

    var stopped = False
    try:
        _ = poly.dense.power([1, 1], 100, max_steps=10000)
    except error:
        stopped = error.kind == PolyErrorKind.LIMIT
    print("Stopped by max_steps=10000:", stopped)

    var result = poly.dense.power([1, 1], 100, max_steps=1000000)
    print("Finished within max_steps=1000000, degree", poly.dense.degree(result))
