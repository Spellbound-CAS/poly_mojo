"""Algebraic bounds remain checked through streaming, overflow and work limits."""
from poly import dense, sparse, algebraic, evaluation
from apn_mojo import Float, BallContext, ArithmeticContext, FloatFormat, Rational, contains
from polynomials import Budget, PolyErrorKind, NumericDomain, EnclosureStream, EnclosureStreamChecker, QQ
from std.testing import assert_true, assert_false, assert_equal


def consume[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](var stream: evaluation.AlgebraicEnclosureStream[C, D, numerical, sparse, False], var checker: evaluation.AlgebraicEnclosureStream[C, D, numerical, sparse, True], mut budget: Budget) raises -> D.Scalar:
    while True:
        var step = stream.next(budget)
        if not step:
            break
        assert_true(checker.accept(step.take(), budget))
    var receipt = stream.finish(budget)
    assert_true(checker.finish(receipt, budget))
    receipt.steps += 1
    assert_false(checker.finish(receipt, budget))
    return receipt.value(budget)


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var bc = BallContext(precision=40)
    var p = dense.asarray([a.copy(), algebraic.asreal(1)])
    var budget = Budget()
    var value = consume(dense.enclose_stream(p, 1, context=bc), dense.enclose_checker(p, 1, context=bc), budget)
    assert_false(contains(value, Rational(5, 2)))
    var z = algebraic.ascomplex(a, 1)
    var sp = sparse.asarray([((1, 0), z.copy())])
    var sv = consume(sparse.enclose_stream(sp, (1, 0), context=bc), sparse.enclose_checker(sp, (1, 0), context=bc), budget)
    assert_true(contains(sv.imag(), Rational(1)))
    var warm = Budget()
    _ = dense.enclose(p, a, warm, context=bc)
    var spent = UInt64.MAX - warm.remaining()
    var exact = Budget(spent)
    var result = dense.enclose(p, a, exact, context=bc)
    assert_equal(exact.remaining(), UInt64(0))
    var rejected = False
    try:
        _ = dense.enclose(p, a, context=bc, max_steps=spent - 1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    warm = Budget()
    assert_true(evaluation.check(result, warm))
    spent = UInt64.MAX - warm.remaining()
    exact = Budget(spent)
    assert_true(evaluation.check(result, exact))
    assert_equal(exact.remaining(), UInt64(0))
    rejected = False
    try:
        _ = evaluation.check(result, max_steps=spent - 1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    # Same polynomial and root counts, but a different isolating interval.
    result._data[].algebraic.value().witnesses[0].real.low = QQ(Rational(3))
    assert_false(evaluation.check(result))
    var overflow = ArithmeticContext(format=FloatFormat(precision=8, emin=-2, emax=2))
    var large = algebraic.multiply(a, 100)
    rejected = False
    try:
        _ = dense.asarray([large.copy()], context=overflow)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    assert_true(algebraic.check(a))
    assert_equal(p[0], a)
    print("array algebraic numerical limit checks passed")
