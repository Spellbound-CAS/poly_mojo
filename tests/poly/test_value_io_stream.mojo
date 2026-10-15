"""Fragmented streaming, bounded fields, and failed transports preserve input values."""
from std.testing import assert_equal, assert_true
from apn_mojo import Integer
from poly import io, series, basis, batch, dense, rational, algebraic
from polynomials import Budget, PolyError, PolyErrorKind, PolynomialTextSource, PolynomialTextSink


@fieldwise_init
struct Sink(PolynomialTextSink):
    var text: String
    var chunks: Int
    var largest: Int

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(UInt64(chunk.byte_length()) + 1)
        self.largest = max(self.largest, chunk.byte_length())
        self.text.write(chunk)
        self.chunks += 1


@fieldwise_init
struct Source(PolynomialTextSource):
    var text: String
    var position: Int
    var width: Int
    var empty: Bool

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        self.empty = not self.empty
        if self.empty:
            return String()
        if self.position == self.text.byte_length():
            return None
        var end = min(self.text.byte_length(), self.position + self.width)
        var chunk = String(self.text[byte=self.position:end])
        self.position = end
        return chunk^


@fieldwise_init
struct BrokenSource(PolynomialTextSource):
    var first: Bool

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        if self.first:
            self.first = False
            return String("POLYNOMIALS-VALUE 1\nSERIES(QQ)\n4\n1\n")
        raise PolyError(PolyErrorKind.ENCODING, "read failed")


@fieldwise_init
struct BrokenSink(PolynomialTextSink):
    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(1)
        raise PolyError(PolyErrorKind.ENCODING, "write failed")


def fragmented[T: Copyable & Deinitable](value: T) raises:
    var sink = Sink(String(), 0, 0)
    io.dump(value, sink)
    var text = io.dumps(value)
    assert_equal(sink.text, text)
    for width in [1, 2, 7, 4096]:
        var copy = io.load[T](Source(text.copy(), 0, width, False))
        assert_equal(io.dumps(copy), text)
    var work = Budget()
    _ = io.load[T](Source(text.copy(), 0, 1, False), work)
    var spent = UInt64.MAX - work.remaining()
    var exact = Budget(spent)
    _ = io.load[T](Source(text.copy(), 0, 1, False), exact)
    assert_equal(exact.remaining(), UInt64(0))
    var exhausted = False
    try:
        _ = io.load[T](Source(text.copy(), 0, 1, False), max_steps=spent - 1)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    assert_true(exhausted)
    var small = False
    try:
        _ = io.load[T](Source(text.copy(), 0, 1, False), max_field_bytes=2)
    except e:
        small = e.kind == PolyErrorKind.CAPACITY
    assert_true(small)


def main() raises:
    var value = series.asarray([1, 2], order=10)
    fragmented(value)
    fragmented(basis.asarray([1, 2, 0], basis="newton", nodes=[0, 1]))
    fragmented(rational.asarray([1], [1, 1]))
    var root = dense.root_between([-2, 0, 1], 1, 2)
    fragmented(algebraic.ascomplex(root, 1))
    var b = batch.broadcast_to(batch.scalar(value), List[Int]([300]))
    var sink = Sink(String(), 0, 0)
    io.dump(b, sink)
    assert_true(sink.chunks > 1)
    assert_true(sink.largest <= 4096 + 32)
    var copy = io.load[type_of(b)](Source(sink.text.copy(), 0, 17, False))
    assert_equal(len(copy), 300)
    assert_equal(copy[299].order(), 10)
    var failed = False
    try:
        _ = io.load[type_of(value)](BrokenSource(True))
    except e:
        failed = e.kind == PolyErrorKind.ENCODING and e.message == "read failed"
    assert_true(failed)
    var broken = BrokenSink()
    failed = False
    try:
        io.dump(value, broken)
    except e:
        failed = e.kind == PolyErrorKind.ENCODING and e.message == "write failed"
    assert_true(failed)
    var digits = String()
    for _ in range(600):
        digits.write("12345678")
    var large_text = String("POLYNOMIALS-VALUE 1\nZZ\n", digits, "\n")
    var large = io.loads[Integer](large_text, max_field_bytes=4800)
    var large_sink = Sink(String(), 0, 0)
    io.dump(large, large_sink)
    assert_equal(large_sink.text, large_text)
    assert_true(large_sink.largest > 4096)
    failed = False
    try:
        _ = io.load[Integer](Source(large_text.copy(), 0, 17, False), max_field_bytes=4799)
    except e:
        failed = e.kind == PolyErrorKind.CAPACITY
    assert_true(failed)
    assert_equal(value.order(), 10)
    print("value serialization stream checks passed")
