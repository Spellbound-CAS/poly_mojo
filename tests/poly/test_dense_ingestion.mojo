"""One-pass dense construction, read limits, and ownership across growth boundaries."""
from apn_mojo import Integer, Rational, ExactComplex
from poly import dense, algebraic, coefficients
from polynomials import Budget, PolyErrorKind
from std.iter import IterableOwned, Iterator, StopIteration, iter
from std.memory import ArcPointer
from std.testing import assert_equal, assert_true


@fieldwise_init
struct Counts(Movable):
    var iterations: Int
    var reads: Int
    var hints: Int
    var destroyed: Int


struct Once[T: Copyable & Deinitable](IterableOwned, Iterator):
    comptime Element = Self.T
    comptime IteratorOwnedType: Iterator = Self
    var values: List[Self.T]
    var position: Int
    var counts: ArcPointer[Counts]

    def __init__(out self, values: List[Self.T], counts: ArcPointer[Counts]):
        self.values = values.copy()
        self.position = 0
        self.counts = counts

    def __deinit__(deinit self):
        self.counts[].destroyed += 1

    def __iter__(var self) -> Self:
        self.counts[].iterations += 1
        return self^

    def __next__(mut self) raises StopIteration -> Self.T:
        self.counts[].reads += 1
        if self.position == len(self.values):
            raise StopIteration()
        var value = self.values[self.position].copy()
        self.position += 1
        return value^

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        self.counts[].hints += 1
        # Deliberately unhelpful: hints must never authorize allocation or reads.
        return Int.MAX, Optional(Int.MAX)


def test_growth_and_zeros() raises:
    var huge = Rational((Integer(1) << 2048) + 1, 3)
    for n in [0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 255, 256, 257]:
        for kept in [0, 1, n // 2, n]:
            var expected = List[Rational]()
            var input = List[Rational]()
            for i in range(n):
                var value = huge.copy() if i < kept else Rational(0)
                input.append(value.copy())
                if i < kept:
                    expected.append(value^)
            var counts = ArcPointer(Counts(0, 0, 0, 0))
            var result = dense.asarray(Once(input, counts))
            assert_equal(result.to_list(), expected)
            assert_equal(counts[].iterations, 1)
            assert_equal(counts[].reads, n + 1)
            assert_equal(counts[].hints, 0)
            assert_equal(counts[].destroyed, 1)
            if n:
                input[0] = Rational(99)
                assert_equal(result.to_list(), expected)
            if len(result) == 0:
                assert_equal(result._core[].coefficients.capacity(), 0)
            else:
                assert_true(result._core[].coefficients.capacity() <= 2 * len(result))


def test_exact_read_limits() raises:
    var values: List[Rational] = [Rational(1, 3), Rational(2, 5), Rational(3, 7)]
    # Hand-accounted trace for three values: read, allocate/move if full,
    # transfer; then EOF, last-coefficient check, and result handle.
    var charges: List[Int] = [1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 1]
    var reads: List[Int] = [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
    for allowance in range(22):
        var remaining = allowance
        var expected_reads = 0
        var success = True
        for i in range(len(charges)):
            if remaining < charges[i]:
                success = False
                break
            remaining -= charges[i]
            expected_reads += reads[i]
        var counts = ArcPointer(Counts(0, 0, 0, 0))
        var budget = Budget(UInt64(allowance))
        var completed = False
        var result = dense.asarray(values)
        try:
            result = dense.asarray(Once(values, counts), budget)
            completed = True
        except error:
            assert_equal(error.kind, PolyErrorKind.LIMIT)
        if completed:
            assert_equal(result.to_list(), values)
        assert_equal(completed, success)
        assert_equal(budget.remaining(), UInt64(remaining))
        assert_equal(counts[].reads, expected_reads)
        assert_equal(counts[].iterations, 1)
        assert_equal(counts[].destroyed, 1)
        assert_equal(counts[].hints, 0)
    # EOF with one unit left must be read even though a transfer would need two.
    var empty_counts = ArcPointer(Counts(0, 0, 0, 0))
    var budget = Budget(1)
    var failed = False
    try:
        _ = dense.asarray(Once(List[Rational](), empty_counts), budget)
    except error:
        failed = error.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    assert_equal(empty_counts[].reads, 1)
    assert_equal(budget.remaining(), UInt64(0))


def check_domain[T: Copyable & Deinitable](values: List[T]) raises:
    var unlimited = Budget()
    var counts = ArcPointer(Counts(0, 0, 0, 0))
    var expected = dense.asarray(Once(values, counts), unlimited)
    var spent = UInt64.MAX - unlimited.remaining()
    for allowance in range(Int(spent) + 1):
        var observed = ArcPointer(Counts(0, 0, 0, 0))
        var budget = Budget(UInt64(allowance))
        var completed = False
        var result = expected.copy()
        try:
            result = dense.asarray(Once(values, observed), budget)
            completed = True
        except error:
            assert_equal(error.kind, PolyErrorKind.LIMIT)
        if completed:
            assert_equal(result.to_list(), expected.to_list())
        assert_equal(completed, UInt64(allowance) == spent)
        assert_true(observed[].reads <= len(values) + 1)
        assert_equal(observed[].destroyed, 1)
        assert_equal(observed[].hints, 0)
        if allowance == 0:
            assert_equal(observed[].reads, 0)
    assert_equal(dense.asarray(iter(values)).to_list(), expected.to_list())


def test_domains() raises:
    var native: List[Int] = [1, -2, 3, 4, 5]
    var unsigned: List[UInt64] = [UInt64.MAX, 0, 1]
    var signed: List[Int64] = [Int64.MIN, 0, 1]
    var integer: List[Integer] = [Integer(1) << 1024, -2, 3]
    var rational: List[Rational] = [Rational(1, 3), Rational(0), Rational(2)]
    var imaginary: List[ExactComplex] = [ExactComplex(0, 1), ExactComplex(0, 0), ExactComplex(1, 2)]
    var real: List[algebraic.Real] = [algebraic.asreal(2), algebraic.asreal(0), algebraic.asreal(-3)]
    var complex_values: List[algebraic.Complex] = [algebraic.ascomplex(0, 1), algebraic.ascomplex(0), algebraic.ascomplex(2, -1)]
    check_domain(native)
    check_domain(unsigned)
    check_domain(signed)
    check_domain(integer)
    check_domain(rational)
    check_domain(imaginary)
    check_domain(real)
    check_domain(complex_values)
    # Mixed-input promotion and modular validation remain compatible.
    var mixed = coefficients((1, Rational(1, 2), ExactComplex(0, 1)))
    assert_equal(dense.asarray(mixed).to_list(), [ExactComplex(1), ExactComplex(Rational(1, 2)), ExactComplex(0, 1)])
    var modular = dense.asarray([1, 2, 3], modulus=7).to_list()
    check_domain(modular)
    var other = dense.asarray([1], modulus=11).to_list()
    modular.append(other[0].copy())
    var failed = False
    try:
        _ = dense.asarray(iter(modular))
    except error:
        failed = error.kind == PolyErrorKind.DOMAIN
    assert_true(failed)


def main() raises:
    test_growth_and_zeros()
    test_exact_read_limits()
    test_domains()
    print("dense ingestion checks passed")
