"""Chunk boundaries, independent core storage, and calculations without flattening."""
from apn_mojo import Integer, Rational, ExactComplex
from poly import dense, series, io, algebraic, modint
from polynomials import Budget, CoefficientBuffer, DensePoly, QQ, NoContext, PolyErrorKind
from std.iter import iter
from std.testing import assert_equal, assert_true, assert_false
from test_dense_ingestion import Counts, Once
from std.memory import ArcPointer


def test_boundaries() raises:
    comptime first = CoefficientBuffer[Int]._initial_size
    comptime prefix = first + CoefficientBuffer[Int]._second_size
    comptime boundary = prefix + CoefficientBuffer[Int]._block_size
    for n in [first - 1, first, first + 1, 1023, 1024, 1025, prefix - 1, prefix, prefix + 1, 4097, boundary - 1, boundary, boundary + 1, boundary + CoefficientBuffer[Int]._block_size + 1]:
        for kept in [0, 1, min(n, first - 1), min(n, first), min(n, first + 1), min(n, prefix - 1), min(n, prefix), min(n, boundary - 1), min(n, boundary), n]:
            var values = List[Rational](capacity=n)
            for i in range(n):
                values.append(Rational(i + 1, 3) if i < kept else Rational(0))
            var counts = ArcPointer(Counts(0, 0, 0, 0))
            var p = dense.asarray(Once(values, counts))
            assert_equal(len(p), kept)
            assert_equal(counts[].iterations, 1)
            assert_equal(counts[].reads, n + 1)
            assert_equal(counts[].hints, 0)
            assert_equal(counts[].destroyed, 1)
            assert_equal(p._core[].coefficients.is_contiguous(), kept <= first)
            if kept <= first:
                assert_equal(p._core[].coefficients._chunks.capacity(), 0)
            assert_true(p._core[].coefficients.capacity() <= 2 * kept)
            for i in range(kept):
                assert_equal(p[i], Rational(i + 1, 3))


def test_storage_and_operations() raises:
    var values = List[Rational](length=2051, fill=Rational(1))
    var p = dense.asarray(iter(values))
    var q = p.copy()
    assert_false(p._core[].coefficients.is_contiguous())
    var work = Budget()
    var copy = p._core[].coefficients.copy()
    copy[1024] = QQ(Rational(7))
    assert_equal(q[1024], Rational(1))
    assert_equal(copy[1024].value, Rational(7))
    var listed = p._core[].coefficients.to_list(work)
    listed[0] = QQ(Rational(9))
    assert_equal(p[0], Rational(1))
    var index = 0
    for coefficient in p._core[].coefficients:
        assert_equal(coefficient.value, Rational(1))
        index += 1
    assert_equal(index, 2051)
    var canonical = DensePoly[QQ].from_coefficients(0, NoContext(), copy^, work)
    assert_equal(canonical.degree(), 2050)
    var added = dense.add(p, p)
    var product = dense.multiply(p, [1, -1])
    assert_equal(len(product), 2052)
    assert_equal(product[0], Rational(1))
    assert_equal(product[2051], Rational(-1))
    for i in range(1, 2051):
        assert_equal(product[i], Rational(0))
        assert_equal(added[i], Rational(2))
    assert_equal(dense.evaluate(p, 1), Rational(2051))
    var derivative = dense.derivative(p)
    assert_equal(derivative[1024], Rational(1025))
    var divided = dense.divmod(product, [1, -1])
    assert_equal(divided[0].to_list(), p.to_list())
    assert_equal(len(divided[1]), 0)
    var formal = series.asarray(p, order=2051)
    assert_equal(series.derivative(formal)[1024], Rational(1025))
    var text = io.dumps(p)
    assert_equal(io.loads[dense.Array[1]](text).to_list(), p.to_list())
    # First and repeated operations read the same chunks; neither may replace
    # the shared input with a deferred contiguous copy.
    assert_equal(dense.evaluate(q, 0), Rational(1))
    assert_false(p._core[].coefficients.is_contiguous())
    assert_false(q._core[].coefficients.is_contiguous())
    var core = DensePoly[QQ].from_coefficients(0, NoContext(), p._core[].coefficients.copy(), work)
    while len(core.coefficients) > 1023:
        _ = core.coefficients.pop()
    core.coefficients.append(QQ(Rational(3)))
    assert_equal(core.coefficients[1023].value, Rational(3))
    assert_equal(p[1023], Rational(1))



def test_bulk_runs() raises:
    # Exercise contiguous-run traversal across two bulk allocations, including
    # the unmatched suffix when the second operand is shorter.
    comptime prefix = CoefficientBuffer[Int]._initial_size + CoefficientBuffer[Int]._second_size
    comptime n = prefix + CoefficientBuffer[Int]._block_size + 3
    var values = List[Rational](length=n, fill=Rational(1))
    var p = dense.asarray(iter(values))
    var doubled = dense.add(p, p)
    var shifted = dense.add(p, [1])
    assert_equal(len(doubled), n)
    assert_equal(len(shifted), n)
    for i in range(n):
        assert_equal(doubled[i], Rational(2))
        assert_equal(shifted[i], Rational(2) if i == 0 else Rational(1))
    assert_false(p._core[].coefficients.is_contiguous())


def check_chunk_domain[T: Copyable & Deinitable](value: T) raises:
    var values = List[T](length=1025, fill=value)
    var p = dense.asarray(iter(values))
    assert_false(p._core[].coefficients.is_contiguous())
    assert_equal(p.to_list(), dense.asarray(values).to_list())


def test_mutable_buffer() raises:
    comptime first = CoefficientBuffer[Int]._initial_size
    comptime second = CoefficientBuffer[Int]._second_size
    comptime bulk = CoefficientBuffer[Int]._block_size
    var chunks = List[List[Int]]()
    chunks.append(List[Int](length=first, fill=7))
    chunks.append(List[Int](length=second, fill=8))
    chunks.append(List[Int](length=bulk, fill=9))
    var buffer = CoefficientBuffer[Int]._from_chunks(chunks^, first + second + bulk)
    buffer.append(10)
    assert_equal(len(buffer), first + second + bulk + 1)
    assert_equal(buffer.pop(), 10)
    for _ in range(bulk):
        assert_equal(buffer.pop(), 9)
    assert_false(buffer.is_contiguous())
    for _ in range(second):
        assert_equal(buffer.pop(), 8)
    assert_true(buffer.is_contiguous())
    assert_equal(buffer._chunks.capacity(), 0)
    assert_equal(len(buffer), first)
    assert_equal(buffer[first - 1], 7)
    buffer.clear()
    assert_equal(len(buffer), 0)
    assert_equal(buffer.capacity(), 0)
    buffer.append(3)
    assert_equal(buffer[0], 3)
    # A compacted tail can grow back to its scheduled length. Its allocation
    # must stop there, so capacity remains accurate after a new chunk is added.
    var short_chunks = List[List[Int]]()
    short_chunks.append(List[Int](length=first, fill=7))
    short_chunks.append(List[Int](length=17, fill=8))
    var growing = CoefficientBuffer[Int]._from_chunks(short_chunks^, first + 17)
    for _ in range(17, second):
        growing.append(8)
    growing.append(9)
    assert_equal(growing._chunks[1].capacity(), second)
    assert_equal(growing.capacity(), first + second + bulk)
    assert_equal(growing[first + second], 9)
    var exact_chunks = List[List[QQ]]()
    exact_chunks.append(List[QQ](length=first, fill=QQ(Rational(1))))
    exact_chunks.append(List[QQ](length=2, fill=QQ(Rational(0))))
    var work = Budget()
    var p = DensePoly[QQ].from_coefficients(0, NoContext(), CoefficientBuffer[QQ]._from_chunks(exact_chunks^, first + 2), work)
    assert_equal(p.degree(), first - 1)
    assert_true(p.coefficients.is_contiguous())
    p.check_invariants(work)


def test_domains() raises:
    check_chunk_domain(Int(-7))
    check_chunk_domain(Integer(1) << 2048)
    check_chunk_domain(Rational(2, 3))
    check_chunk_domain(ExactComplex(1, 2))
    check_chunk_domain(algebraic.asreal(2))
    check_chunk_domain(algebraic.ascomplex(0, 1))
    check_chunk_domain(modint(3, 7))


def test_limits() raises:
    comptime n = CoefficientBuffer[Int]._initial_size + CoefficientBuffer[Int]._second_size + 1
    var values = List[Rational](length=n, fill=Rational(1, 3))
    var full = Budget()
    _ = dense.asarray(iter(values), full)
    var spent = UInt64.MAX - full.remaining()
    # Sample the whole allowance range and cover every finalization charge.
    # The small ingestion regression separately checks every unit of its trace.
    var allowances = List[Int]()
    for allowance in range(0, Int(spent), 31):
        allowances.append(allowance)
    for allowance in range(3000, 3300):
        allowances.append(allowance)
    for allowance in range(Int(spent) - 32, Int(spent) + 1):
        allowances.append(allowance)
    for allowance in allowances:
        var counts = ArcPointer(Counts(0, 0, 0, 0))
        var budget = Budget(UInt64(allowance))
        var success = False
        var observed_length = 0
        try:
            var result = dense.asarray(Once(values, counts), budget)
            success = True
            observed_length = len(result)
        except error:
            assert_equal(error.kind, PolyErrorKind.LIMIT)
        if success:
            assert_equal(observed_length, n)
        assert_equal(success, UInt64(allowance) == spent)
        assert_true(counts[].reads <= n + 1)
        assert_equal(counts[].destroyed, 1)
        assert_equal(counts[].hints, 0)
    var entries: List[Int] = [1, 2, 3]
    var buffer = CoefficientBuffer(entries^)
    var budget = Budget(5)
    var failed = False
    try:
        _ = buffer.to_list(budget)
    except error:
        failed = error.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    assert_equal(buffer[2], 3)


def main() raises:
    test_boundaries()
    test_storage_and_operations()
    test_bulk_runs()
    test_domains()
    test_mutable_buffer()
    test_limits()
    print("chunked dense checks passed")
