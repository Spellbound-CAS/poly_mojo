"""Point iteration, shape validation, and runtime field errors at the array boundary."""
import poly
from apn_mojo import Integer, Rational, Batch
from polynomials import Budget, PolyErrorKind, ModInt
from std.iter import Iterable, Iterator, StopIteration
from std.memory import ArcPointer
from std.testing import assert_equal, assert_true


@fieldwise_init
struct Counts(Movable):
    var iterations: Int
    var reads: Int


@fieldwise_init
struct Points(Copyable, Movable, Iterable):
    var counts: ArcPointer[Counts]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = PointIterator

    def __iter__(self) -> PointIterator:
        self.counts[].iterations += 1
        return PointIterator(self.counts, 0)


@fieldwise_init
struct PointIterator(Iterator):
    comptime Element = Int
    var counts: ArcPointer[Counts]
    var index: Int

    def __next__(mut self) raises StopIteration -> Int:
        self.counts[].reads += 1
        if self.index == 3:
            raise StopIteration()
        self.index += 1
        return 2 if self.index == 1 else 0

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        # A hint cannot cause an unchecked upfront allocation.
        return Int.MAX, Optional(Int.MAX)


def test_single_pass() raises:
    var p = poly.dense.asarray([0, 1])
    var sparse = poly.to_sparse(p)
    for scenario in range(3):
        var counts = ArcPointer(Counts(0, 0))
        var points = Points(counts)
        if scenario == 0:
            assert_equal(poly.dense.evaluate_many(p, points).to_list(), [Integer(2), Integer(0), Integer(0)])
        elif scenario == 1:
            assert_equal(poly.sparse.evaluate_many(sparse, points).to_list(), [Integer(2), Integer(0), Integer(0)])
        else:
            assert_equal(poly.dense.from_roots(points).to_list(), poly.dense.asarray([0, 0, -2, 1]).to_list())
        assert_equal(counts[].iterations, 1)
        assert_equal(counts[].reads, 4)
        counts[].reads = 0
        var caught = False
        try:
            if scenario == 0:
                _ = poly.dense.evaluate_many(p, points, max_steps=0)
            elif scenario == 1:
                _ = poly.sparse.evaluate_many(sparse, points, max_steps=0)
            else:
                _ = poly.dense.from_roots(points, max_steps=0)
        except error:
            assert_equal(error.kind, PolyErrorKind.LIMIT)
            caught = True
        assert_true(caught)
        assert_equal(counts[].reads, 0)


def test_fields_and_shapes() raises:
    var grid = Batch[Integer]([1, 2, 3, 4]).reshape([2, 2])
    var mixed_fields = [poly.modint(0, 7), poly.modint(1, 5)]
    for scenario in range(9):
        var caught = False
        try:
            if scenario == 0:
                _ = poly.dense.evaluate_many([0, 1], grid)
            elif scenario == 1:
                _ = poly.sparse.evaluate_many([(1, 1)], grid)
            elif scenario == 2:
                _ = poly.dense.from_roots(grid)
            elif scenario == 3:
                _ = poly.dense.evaluate_many([0, 1], mixed_fields)
            elif scenario == 4:
                _ = poly.sparse.evaluate_many([(1, 1)], mixed_fields)
            elif scenario == 5:
                _ = poly.dense.from_roots(mixed_fields)
            elif scenario == 6:
                _ = poly.dense.evaluate_many([0, 1], List[ModInt]())
            elif scenario == 7:
                _ = poly.sparse.evaluate_many([(1, 1)], List[ModInt]())
            else:
                _ = poly.dense.from_roots(List[ModInt]())
        except error:
            assert_equal(error.kind, PolyErrorKind.DOMAIN)
            caught = True
        assert_true(caught)
    var coordinates = [(1, poly.modint(3, 7)), (0, poly.modint(0, 7))]
    var result = poly.sparse.evaluate_many([((1, 0), 2), ((0, 1), 3)], coordinates)
    assert_equal(result.modulus(), UInt64(7))
    assert_equal(result[0].value, UInt64(4))
    assert_equal(result[1].value, UInt64(0))
    var empty = List[Tuple[Int, Rational]]()
    var values = poly.sparse.evaluate_many([((0, 0), 1)], empty)
    assert_equal(values.domain(), "QQ")
    assert_equal(len(values), 0)


def main() raises:
    test_single_pass()
    test_fields_and_shapes()
    print("array evaluation checks passed")
