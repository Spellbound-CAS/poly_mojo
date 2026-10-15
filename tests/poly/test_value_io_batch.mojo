"""Batch records retain logical shape and serialize broadcasting in row-major order."""
from std.testing import assert_equal
from poly import io, batch, dense, series, basis, rational


def same[T: Copyable & Deinitable](value: T) raises:
    var text = io.dumps(value)
    assert_equal(io.dumps(io.loads[T](text)), text)


def main() raises:
    var p = dense.asarray([1, 2])
    var q = dense.asarray([3])
    var a = batch.asarray([p.copy(), q.copy()], shape=List[Int]([2, 1]))
    var broadcast = batch.broadcast_to(a, List[Int]([2, 3]))
    var copy = io.loads[type_of(broadcast)](io.dumps(broadcast))
    assert_equal(copy.shape(), List[Int]([2, 3]))
    for i in range(6):
        assert_equal(copy[i][0], 1 if i < 3 else 3)
    same(broadcast)
    same(batch.scalar(p))
    same(batch.broadcast_to(batch.scalar(p), List[Int]([0, 3])))
    same(batch.asarray([batch.scalar(p), batch.scalar(q)]))
    same(batch.asarray([series.asarray([1], order=3), series.asarray([2], order=7)]))
    same(batch.asarray([basis.asarray([1, 0], basis="bernstein"), basis.asarray([2], basis="chebyshev")]))
    same(batch.asarray([rational.asarray([1], [1, 1]), rational.asarray([2], [1])]))
    same(batch.asarray([1, -2, 3]))
    same(batch.asarray(List[series.Array[1]](), shape=List[Int]([2, 0, 3])))
    var golden = String("POLYNOMIALS-VALUE 1\nBATCH(Int)\n0\n1\nInt\n7\n")
    assert_equal(io.dumps(batch.scalar(7)), golden)
    assert_equal(io.loads[batch.Array[Int]](golden).ndim(), 0)
    print("value serialization batch checks passed")
