"""Consume polynomial text and proofs without retaining completed payloads."""
from apn_mojo import Rational
from polynomials import *


@fieldwise_init
struct ByteCounter(PolynomialTextSink):
    var count: Int

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(UInt64(chunk.byte_length()))
        self.count += chunk.byte_length()


def main() raises:
    var work = Budget()
    var polynomial = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], work)
    var text_sink = ByteCounter(0)
    serialize_dense_to(polynomial, text_sink, work, chunk_bytes=16)
    if text_sink.count != serialize_dense(polynomial, work).byte_length():
        raise Error("text byte count mismatch")
    print("Canonical text consumed in chunks")

    var options = NumericOptions(absolute_bits=2, max_subdivisions=1024)
    var real_checker = DescartesStreamChecker.from_dense(polynomial, options, work)
    var real_receipt = real_roots_exact_to(polynomial, options, real_checker, work)
    if not real_checker.finish(real_receipt, work):
        raise Error("incomplete real root witness")
    var complex_checker = ComplexRootStreamChecker.from_dense(polynomial, options, work)
    var complex_receipt = complex_roots_exact_to(polynomial, options, complex_checker, work)
    if not complex_checker.finish(complex_receipt, work):
        raise Error("incomplete complex root witness")
    print("Real and complex root streams checked")

    var inputs: List[SparsePoly[QQ]] = [SparsePoly[QQ].variable(0, MonomialOrder(), NoContext(), work)]
    var cad_checker = CADStreamChecker.prepare(inputs, work)
    var cad_receipt = cad_to(inputs, [Int64(0)], cad_checker, work)
    if not cad_checker.finish(cad_receipt, work):
        raise Error("incomplete CAD witness")
    print("Cylindrical coverage checked incrementally")
