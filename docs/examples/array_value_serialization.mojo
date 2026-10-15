"""Save complete series, basis representations, fractions, and batches."""
from poly import io, series, basis, rational, batch
from apn_mojo import ArithmeticContext, FloatFormat


def main() raises:
    var s = series.asarray([1, 2], order=6)
    var restored = io.loads[type_of(s)](io.dumps(s))
    print("Known series coefficients:", restored.order())

    var b = basis.asarray([1, 2, 0], basis="newton", nodes=[0, 3])
    var saved_basis = io.loads[type_of(b)](io.dumps(b))
    print("Basis:", saved_basis.basis(), "coefficient count:", len(saved_basis))
    print("Last Newton node:", basis.nodes(saved_basis)[1])

    var f = rational.asarray([-1, 0, 1], [-1, 1])
    var saved_fraction = io.loads[type_of(f)](io.dumps(f))
    print("Restored fraction at three:", rational.evaluate(saved_fraction, 3))
    print("Fraction verified:", rational.check(saved_fraction))

    var grid = batch.broadcast_to(batch.scalar(s), List[Int]([2, 3]))
    var saved_grid = io.loads[type_of(grid)](io.dumps(grid))
    print("Batch shape:", saved_grid.shape()[0], "by", saved_grid.shape()[1])

    var context = ArithmeticContext(format=FloatFormat(precision=96))
    var numerical = series.asarray([1], order=4, context=context)
    var saved_numerical = io.loads[type_of(numerical)](io.dumps(numerical))
    print("Implicit zero precision:", saved_numerical[3].precision())
