"""Incremental evaluation, numeric roots, fragmented transport and corruption."""
from std.testing import assert_equal, assert_true, assert_false
from apn_mojo import Rational, ExactComplex, Float, Ball, ComplexBall, BallContext, ArithmeticContext, FloatFormat, contains
from polynomials import Budget, PolyError, PolyErrorKind, NumericDomain, NumericOptions, NumericBox, QQ, EnclosureStream, EnclosureStreamChecker, PolynomialTextSource, PolynomialTextSink
from poly import dense, sparse, roots


@fieldwise_init
struct Sink(PolynomialTextSink):
    var text: String
    var chunks: Int

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(UInt64(chunk.byte_length()) + 1)
        self.text.write(chunk)
        self.chunks += 1


@fieldwise_init
struct Source(PolynomialTextSource):
    var text: String
    var position: Int

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        if self.position == self.text.byte_length():
            return None
        var chunk = String(self.text[byte=self.position:self.position + 1])
        self.position += 1
        return chunk^


def consume[D: NumericDomain](var stream: EnclosureStream[D], var checker: EnclosureStreamChecker[D], mut budget: Budget) raises -> D.Scalar:
    var premature = False
    try:
        _ = stream.finish(budget)
    except e:
        premature = e.kind == PolyErrorKind.DOMAIN
    assert_true(premature)
    var count = 0
    while True:
        var step = stream.next(budget)
        if not step:
            break
        assert_true(checker.accept(step.take(), budget))
        count += 1
    var receipt = stream.finish(budget)
    assert_equal(count, receipt.steps)
    assert_true(checker.finish(receipt, budget))
    receipt.steps += 1
    assert_false(checker.finish(receipt, budget))
    return receipt.value(budget)


def main() raises:
    var budget = Budget()
    var bc = BallContext(precision=80)
    var fc = ArithmeticContext(format=FloatFormat(precision=80))
    var value = consume(dense.enclose_stream([1, 2, 3], Rational(1, 3), context=bc), dense.enclose_checker([1, 2, 3], Rational(1, 3), context=bc), budget)
    assert_true(contains(value, Rational(2)))
    var complex = [ExactComplex(1, 0), ExactComplex(0, 1)]
    var cv = consume(dense.enclose_stream(complex, 2, context=bc), dense.enclose_checker(complex, 2, context=bc), budget)
    assert_true(contains(cv.real(), Rational(1)))
    assert_true(contains(cv.imag(), Rational(2)))
    var p = [((1, 0), 2), ((0, 1), 3)]
    var sv = consume(sparse.enclose_stream(p, (1, 2), context=bc), sparse.enclose_checker(p, (1, 2), context=bc), budget)
    assert_true(contains(sv, Rational(8)))
    var numeric = dense.asarray([1, 2, 3], context=fc)
    var nv = consume(dense.enclose_stream(numeric, 2, context=bc), dense.enclose_checker(numeric, 2, context=bc), budget)
    assert_true(contains(nv, Rational(17)))
    var sn = sparse.asarray(p, context=bc)
    var nsv = consume(sparse.enclose_stream(sn, (1, 2), context=bc), sparse.enclose_checker(sn, (1, 2), context=bc), budget)
    assert_true(contains(nsv, Rational(8)))
    var single = sparse.asarray([(1, 2)], context=fc)
    var usv = consume(sparse.enclose_stream(single, 3, context=bc), sparse.enclose_checker(single, 3, context=bc), budget)
    assert_true(contains(usv, Rational(6)))
    var exact_single = consume(sparse.enclose_stream([(1, 2)], 3, context=bc), sparse.enclose_checker([(1, 2)], 3, context=bc), budget)
    assert_true(contains(exact_single, Rational(6)))
    var bad = dense.enclose_checker([1, 2], 3, context=bc)
    assert_false(bad.accept(NumericBox.point(QQ(Rational(999)), QQ(Rational(0)), budget), budget))
    var source = dense.enclose_stream([1, 2], 3, context=bc)
    var next_box = source.next(budget)
    assert_false(bad.accept(next_box.take(), budget))
    var truncated = dense.enclose_checker([1, 2], 3, context=bc)
    while source.next(budget):
        pass
    assert_false(truncated.finish(source.finish(budget), budget))
    var options = NumericOptions(absolute_bits=2, max_subdivisions=256)
    var real = dense.asarray([-2, 0, 1], context=fc)
    var rc = dense.real_root_checker(real, options=options)
    var rr = dense.real_root_cover_to(real, rc, options=options, complete_only=True)
    assert_true(rc.finish(rr, budget))
    var cp = dense.asarray([1, 0, 1], context=bc)
    var cc = dense.complex_root_checker(cp, options=options)
    var cr = dense.complex_root_cover_to(cp, cc, options=options, complete_only=True)
    assert_true(cc.finish(cr, budget))
    var sr = sparse.asarray([((0,), -2), ((2,), 1)], context=fc)
    var src = sparse.real_root_checker(sr, context=fc, options=options)
    var srr = sparse.real_root_cover_to(sr, src, context=fc, options=options, complete_only=True)
    assert_true(src.finish(srr, budget))
    assert_true(roots.check(sparse.real_root_cover(sr, context=fc, options=options)))
    var sc = sparse.asarray([(0, 1), (2, 1)], context=bc)
    var scc = sparse.complex_root_checker(sc, context=bc, options=options)
    var scr = sparse.complex_root_cover_to(sc, scc, context=bc, options=options, complete_only=True)
    assert_true(scc.finish(scr, budget))
    assert_true(roots.check(sparse.complex_root_cover(sc, context=bc, options=options)))
    var sink = Sink("", 0)
    dense.dump(numeric, sink, chunk_bytes=1)
    assert_true(sink.chunks > 1)
    assert_equal(sink.text, dense.dumps(numeric))
    var restored = dense.load_numeric[Float](Source(sink.text.copy(), 0))
    assert_equal(dense.dumps(restored), sink.text)
    var ssink = Sink("", 0)
    sparse.dump(sn, ssink, chunk_bytes=2)
    var restored_sparse = sparse.load_numeric[Ball, ndim=2](Source(ssink.text.copy(), 0))
    assert_equal(sparse.dumps(restored_sparse), sparse.dumps(sn))
    var failed = False
    try:
        _ = dense.load_numeric[Float](Source(sink.text + "garbage", 0))
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    failed = False
    try:
        _ = dense.load_numeric[Float](Source("", 0))
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    failed = False
    try:
        dense.dump(numeric, sink, chunk_bytes=0)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    print("array parity stream checks passed")
