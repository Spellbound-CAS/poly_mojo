"""Use APN's array layout operations around a polynomial batch reduction."""
from apn_mojo import Batch, Integer
from poly import batch
from polynomials import Budget, PolyError, ZZ
from std.testing import assert_equal


def add_integers(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    return ZZ(a).add(ZZ(b), budget).value


def main() raises:
    var source = Batch[Integer]([1, 2, 3, 4, 5, 6]).reshape([2, 3])
    var arranged = source.transpose()
    var prepared = batch.asarray(arranged, shape=arranged.shape())
    assert_equal(prepared.shape(), List[Int]([3, 2]))
    var sums = batch.reduce(add_integers, prepared, Integer(0), axis=1)
    var restored = Batch[Integer](sums).reshape(sums.shape())
    assert_equal(restored.shape(), List[Int]([3]))
    assert_equal(restored[0], Integer(5))
    assert_equal(restored[1], Integer(7))
    assert_equal(restored[2], Integer(9))
    print("Column sums:", restored[0], restored[1], restored[2])
