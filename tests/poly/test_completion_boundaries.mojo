"""Degenerate determinants, invalid fields, quotient poles, and budgets."""
from poly import dense, exact, modint, rational, io
from apn_mojo import ExactComplex, BallContext, Rational
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true, assert_false

def main() raises:
    var p = [ExactComplex(1), ExactComplex(0), ExactComplex(1)]
    var q = [ExactComplex(0), ExactComplex(1)]
    var proof = dense.resultant_certified(p, q)
    assert_true(exact.check(proof))
    assert_equal(exact.value(proof)[0], ExactComplex(1))
    var singular = dense.resultant_certified(p, p)
    assert_true(exact.check(singular))
    assert_equal(len(exact.value(singular)), 0)
    var constant = dense.resultant_certified([ExactComplex(2)], p)
    assert_true(exact.check(constant))
    assert_equal(exact.value(constant)[0], ExactComplex(4))
    assert_true(exact.check(dense.resultant_certified(List[ExactComplex](), p)))
    proof._data[].proof.steps[0].pivot_row = 100
    assert_false(exact.check(proof))
    var failed = False
    try:
        _ = dense.resultant_certified(p, q, max_steps=0)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var f = rational.asarray([modint(1, 5)], [modint(1, 5), modint(1, 5)])
    assert_true(rational.check(rational.compose(f, f)))
    assert_equal(rational.residue(f, modint(4, 5)).value, UInt64(1))
    var parts = rational.partial_fractions(f)
    parts._factors[].certificates[0].prime = 7
    assert_false(rational.check(parts))
    var context = BallContext(precision=64)
    var real = rational.asarray([1], [1, 1], context=context)
    for mode in range(3):
        failed = False
        try:
            if mode == 0: _ = rational.enclose(real, -1, context=context)
            elif mode == 1: _ = rational.power(real, 2, context=context, max_steps=0)
            else: _ = rational.asarray([1], [0], context=context)
        except e:
            failed = e.kind == (PolyErrorKind.LIMIT if mode == 1 else PolyErrorKind.DOMAIN)
        assert_true(failed)
    var enclosed = rational.enclose(real, 0, context=context)
    enclosed._value[] = rational.enclose(real, 1, context=context).value()
    assert_false(rational.check(enclosed))
    print("completion boundary checks passed")
