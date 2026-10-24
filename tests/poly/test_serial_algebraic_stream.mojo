"""Fragmented transports, bounded fields and complete-call work limits."""
from poly import dense, sparse, algebraic
from polynomials import Budget, PolyError, PolyErrorKind, PolynomialTextSource, PolynomialTextSink
from std.testing import assert_equal, assert_true


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
    var calls: Int

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        self.calls += 1
        if self.calls == 1:
            return String("POLYNOMIALS-ALGEBRAIC 1\nDENSE\nAA\n0\n1\n")
        raise PolyError(PolyErrorKind.ENCODING, "transport failure")


@fieldwise_init
struct BrokenSink(PolynomialTextSink):
    var calls: Int

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.calls += 1
        raise PolyError(PolyErrorKind.ENCODING, "transport failure")


def main() raises:
    var root = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray([root.copy(), algebraic.asreal(1)])
    var sink = Sink(String(), 0, 0)
    dense.dump(p, sink)
    assert_equal(sink.text, dense.dumps(p))
    for width in [1, 2, 7, 4096]:
        var restored = dense.load[algebraic.Real](Source(sink.text.copy(), 0, width, False))
        assert_equal(restored.to_list(), p.to_list())
    var z = algebraic.ascomplex(root, 1)
    var sp = sparse.asarray([((1, 0), z.copy()), ((0, 2), algebraic.ascomplex(1))])
    var ss = Sink(String(), 0, 0)
    sparse.dump(sp, ss)
    assert_equal(ss.text, sparse.dumps(sp))
    assert_equal(sparse.load[algebraic.Complex, ndim=2](Source(ss.text.copy(), 0, 1, False)).to_list(), sp.to_list())
    var too_large = False
    try:
        _ = dense.load[algebraic.Real](Source(sink.text.copy(), 0, 1, False), max_field_bytes=2)
    except e:
        too_large = e.kind == PolyErrorKind.CAPACITY
    assert_true(too_large)
    var work = Budget()
    _ = dense.dumps(p, work)
    var spent = UInt64.MAX - work.remaining()
    var exact = Budget(spent)
    assert_equal(dense.dumps(p, exact), sink.text)
    assert_equal(exact.remaining(), UInt64(0))
    var limited = False
    try:
        _ = dense.dumps(p, max_steps=spent - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    work = Budget()
    _ = dense.load[algebraic.Real](Source(sink.text.copy(), 0, 1, False), work)
    spent = UInt64.MAX - work.remaining()
    exact = Budget(spent)
    _ = dense.load[algebraic.Real](Source(sink.text.copy(), 0, 1, False), exact)
    assert_equal(exact.remaining(), UInt64(0))
    limited = False
    try:
        _ = dense.load[algebraic.Real](Source(sink.text.copy(), 0, 1, False), max_steps=spent - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(p[0] == root)
    var many = List[algebraic.Real](length=512, fill=algebraic.asreal(1))
    var large = Sink(String(), 0, 0)
    dense.dump(many, large)
    assert_true(large.chunks > 1)
    assert_true(large.largest <= 4096 + 24)
    assert_equal(dense.load[algebraic.Real](Source(large.text.copy(), 0, 17, False)).to_list(), many)
    var failed = False
    try:
        _ = dense.load[algebraic.Real](BrokenSource(0))
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    var broken = BrokenSink(0)
    failed = False
    try:
        dense.dump(p, broken)
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    assert_equal(broken.calls, 1)
    failed = False
    try:
        _ = dense.load[algebraic.Real](Source(sink.text.copy(), 0, 1, False), max_field_bytes=0)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    print("array algebraic serialization stream checks passed")
