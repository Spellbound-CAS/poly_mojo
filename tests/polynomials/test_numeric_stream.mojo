"""Incremental witnesses, independent consumers and fail-closed completion."""
from apn_mojo import Rational, BallContext
from polynomials import *
from polynomials.numeric_evaluate import _enclosure_context
from support.check import Checker


def consume[D: NumericDomain](var stream: EnclosureStream[D], var verifier: EnclosureStreamChecker[D], expected: EnclosureEvaluation[D], mut check: Checker, mut work: Budget) raises:
    check.check(stream.step_count() == len(expected.steps) and verifier.step_count() == stream.step_count(), "stream and checker retain the original witness length")
    try:
        _ = stream.finish(work)
        check.check(False, "incomplete producer returned a completed receipt")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "incomplete producer has no receipt")
    var receipt = StreamedEnclosure[D](expected.precision, len(expected.steps), [expected.value(work)])
    check.check(not verifier.finish(receipt, work), "empty checker rejects a complete-looking receipt")
    for allowance in [UInt64(0), UInt64(1)]:
        var tiny = Budget(allowance)
        try:
            _ = stream.next(tiny)
            check.check(False, "producer ignored budget")
        except e:
            check.check(e.kind == PolyErrorKind.LIMIT, "producer preserves LIMIT before and inside a step, and can retry")
    var i = 0
    while True:
        var next = stream.next(work)
        if not next:
            break
        var box = next.take()
        check.check(box == expected.steps[i], "streamed step equals original materialized witness")
        if i == 0:
            for allowance in [UInt64(0), UInt64(1)]:
                var tiny = Budget(allowance)
                try:
                    _ = verifier.accept(box.copy(), tiny)
                    check.check(False, "checker ignored budget")
                except e:
                    check.check(e.kind == PolyErrorKind.LIMIT, "checker can retry budget failure before and inside a step")
        check.check(verifier.accept(box^, work), "independent exact checker accepts streamed bound")
        i += 1
    receipt = stream.finish(work)
    check.check(verifier.finish(receipt, work), "complete stream and receipt accepted")
    check.check(numeric_bounds[D](receipt.value(work), work) == numeric_bounds[D](expected.value(work), work), "streamed final scalar equals materialized value")
    check.check(not stream.next(work), "producer EOF is stable")
    var wrong = receipt.copy()
    wrong.steps -= 1
    check.check(not verifier.finish(wrong, work), "wrong receipt length rejected")
    wrong = receipt.copy()
    wrong.precision += 1
    check.check(not verifier.finish(wrong, work), "wrong receipt precision rejected")
    wrong = receipt.copy()
    wrong._value = List[D.Scalar]()
    check.check(not verifier.finish(wrong, work), "malformed value storage rejected")
    wrong = receipt.copy()
    wrong._value = [numeric_value[D](NumericBox.point(QQ(Rational(999)), QQ(Rational(0)), work), _enclosure_context[D](receipt.precision, 8192), work)]
    check.check(not verifier.finish(wrong, work), "wrong reported scalar rejected independently")
    check.check(not verifier.accept(expected.steps[0].copy(), work) and not verifier.finish(receipt, work), "extra bound permanently invalidates completion")


def domain[D: NumericDomain](mut check: Checker, mut work: Budget) raises:
    var context = _enclosure_context[D](32, 8192)
    var x = numeric_value[D](NumericBox.real_interval(QQ(Rational(1, 2)), QQ(Rational(3, 4)), work), context, work)
    for n in [0, 1, 20, 21, 65, 257]:
        var coefficients = List[QQ]()
        for i in range(n):
            coefficients.append(QQ(Rational(1 if i % 2 else -1, 3)))
        var poly = DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work)
        var expected = enclose_dense[D](poly, x, context, work)
        consume[D](EnclosureStream[D].from_dense(poly, x, context, work), EnclosureStreamChecker[D].from_dense(poly, x, work), expected, check, work)
        var numerical = numeric_from_dense[D, QQ](poly, context, work)
        var numerical_proof = enclose_numeric[D](numerical, x, context, work)
        consume[D](EnclosureStream[D].from_numeric(numerical, x, context, work), EnclosureStreamChecker[D].from_numeric(numerical, x, work), numerical_proof, check, work)
        var sparse = poly.to_sparse(MonomialOrder(MonomialOrder.LEX, [Int64(0)], work), work)
        var assignments: List[Tuple[Int64, D.Scalar]] = [(Int64(0), x.copy())]
        var sparse_proof = enclose_sparse[D](sparse, assignments, context, work)
        consume[D](EnclosureStream[D].from_sparse(sparse, assignments, context, work), EnclosureStreamChecker[D].from_sparse(sparse, assignments, work), sparse_proof, check, work)
        var numerical_sparse = numeric_from_sparse[D, QQ](sparse, context, work)
        var numerical_sparse_proof = enclose_numeric_sparse[D](numerical_sparse, assignments, context, work)
        consume[D](EnclosureStream[D].from_numeric_sparse(numerical_sparse, assignments, context, work), EnclosureStreamChecker[D].from_numeric_sparse(numerical_sparse, assignments, work), numerical_sparse_proof, check, work)
        var corrupt = EnclosureStreamChecker[D].from_dense(poly, x, work)
        var good = EnclosureStream[D].from_dense(poly, x, context, work)
        var bad = NumericBox.point(QQ(Rational(999)), QQ(Rational(0)), work)
        check.check(not corrupt.accept(bad^, work), "forged first step rejected")
        var first = good.next(work)
        check.check(not corrupt.accept(first.take(), work), "valid later bound cannot repair an invalid stream")
        corrupt = EnclosureStreamChecker[D].from_dense(poly, x, work)
        bad = NumericBox(QQ(Rational(2)), QQ(Rational(1)), QQ(Rational(0)), QQ(Rational(0)))
        check.check(not corrupt.accept(bad^, work), "reversed endpoint interval rejected")


def main() raises:
    var check = Checker("streamed numerical witnesses")
    var work = Budget()
    domain[BallDomain](check, work)
    domain[ComplexBallDomain](check, work)
    check.finish()
