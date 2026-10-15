"""Canonical bytes and adversarial chunk boundaries for interchange transport."""
from apn_mojo import Rational, Integer, ExactComplex, ArithmeticContext, ComplexContext, BallContext
from polynomials import *
from support.check import Checker


@fieldwise_init
struct Source(PolynomialTextSource):
    var text: String
    var position: Int
    var width: Int

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        if self.position == self.text.byte_length():
            return None
        var end = self.position + min(self.width, self.text.byte_length() - self.position)
        var chunk = String(self.text[byte=self.position:end])
        self.position = end
        return chunk^


@fieldwise_init
struct Sink(PolynomialTextSink):
    var text: String
    var calls: Int

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.text.write(chunk)
        self.calls += 1


def exact[C: CoefficientRing](value: C, mut c: Checker, mut work: Budget) raises:
    var p = DensePoly[C].from_coefficients(-7, value.context(), [value.copy(), value.copy()], work)
    var sparse = p.to_sparse(MonomialOrder(), work)
    var sink = Sink(String(), 0)
    serialize_dense_to(p, sink, work, chunk_bytes=1)
    c.check(sink.text == serialize_dense(p, work) and sink.calls >= 3, "streamed dense bytes identical")
    for width in [1, 2, 7, 10000]:
        c.check(deserialize_dense_from[C](Source(sink.text.copy(), 0, width), value.context(), work) == p, "dense arbitrary chunk boundary")
    sink = Sink(String(), 0)
    serialize_sparse_to(sparse, sink, work, chunk_bytes=3)
    c.check(sink.text == serialize_sparse(sparse, work), "sparse bytes identical")
    for width in [1, 2, 7, 10000]:
        c.check(deserialize_sparse_from[C](Source(sink.text.copy(), 0, width), value.context(), work) == sparse, "sparse arbitrary chunk boundary")
    for text in [String(sink.text, "x"), String(sink.text[byte=0:sink.text.byte_length() - 1]), String("POLYNOMIALS 1\nSPARSE\nQQ\n0\n0\n9223372036854775807\n")]:
        try:
            _ = deserialize_sparse_from[C](Source(text, 0, 1), value.context(), work)
            c.check(False, "malformed chunk stream accepted")
        except e:
            c.check(e.kind == PolyErrorKind.ENCODING, "malformed/truncated stream rejected")
    try:
        _ = deserialize_sparse_from[C](Source(sink.text.copy(), 0, 1), value.context(), work, max_field_bytes=2)
        c.check(False, "field limit ignored")
    except e:
        c.check(e.kind == PolyErrorKind.CAPACITY, "bounded field parser")


def numerical[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    var value = numeric_value[D](NumericBox.point(QQ(Rational(3, 2)), QQ(Rational(0)), work), context, work)
    var p = NumericDensePoly[D].from_coefficients(-7, [value.copy(), value.copy()], work)
    var sink = Sink(String(), 0)
    serialize_numeric_dense_to(p, sink, work, chunk_bytes=8)
    c.check(sink.text == serialize_numeric_dense(p, work), "numeric stream preserves scalar JSON")
    var parsed = deserialize_numeric_dense_from[D](Source(sink.text.copy(), 0, 1), work)
    c.check(serialize_numeric_dense(parsed, work) == sink.text, "numeric scalar JSON split at every byte")
    var sparse = NumericSparsePoly[D].from_terms([NumericTerm[D](Monomial(), value.copy())], MonomialOrder(), context, work)
    sink = Sink(String(), 0)
    serialize_numeric_sparse_to(sparse, sink, work, chunk_bytes=8)
    var parsed_sparse = deserialize_numeric_sparse_from[D](Source(sink.text.copy(), 0, 1), work)
    c.check(serialize_numeric_sparse(parsed_sparse, work) == sink.text, "numeric sparse stream")


def main() raises:
    var c = Checker("text stream")
    var work = Budget()
    exact(ZZ(Integer(17)), c, work)
    exact(QQ(Rational(-5, 7)), c, work)
    exact(QQi(ExactComplex(Rational(2, 3), Rational(-4, 5))), c, work)
    exact(ModInt(7, Modulus(101, work)), c, work)
    numerical[FloatDomain](ArithmeticContext(), c, work)
    numerical[ComplexDomain](ComplexContext(ArithmeticContext()), c, work)
    numerical[BallDomain](BallContext(128), c, work)
    numerical[ComplexBallDomain](BallContext(128), c, work)
    c.finish()
