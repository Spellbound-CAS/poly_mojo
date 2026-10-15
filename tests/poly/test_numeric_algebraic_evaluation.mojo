"""Exact coefficient and point witnesses accompany independent numerical traces."""
from poly import dense, sparse, algebraic, evaluation, coefficients
from apn_mojo import Ball, BallContext, Rational, contains
from polynomials import NumericOptions, NumericOutcome
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var bc = BallContext(precision=80)
    var p = dense.asarray(coefficients((a.copy(), 1)))
    var e = dense.enclose(p, 2, context=bc)
    assert_true(evaluation.check(e))
    var square = dense.enclose([algebraic.multiply(a, a)], 0, context=bc)
    assert_true(evaluation.check(square))
    assert_true(contains(square.value(), Rational(2)))
    var options = NumericOptions(initial_precision=8, max_precision=128, absolute_bits=60)
    var adaptive = dense.adaptive_evaluate(p, 2, options=options)
    assert_true(evaluation.check(adaptive))
    assert_true(adaptive.accurate())
    assert_true(adaptive.precision() > 8)
    var ceiling = dense.adaptive_evaluate(p, 2, options=NumericOptions(initial_precision=8, max_precision=8, absolute_bits=60))
    assert_true(evaluation.check(ceiling))
    assert_equal(ceiling.outcome().value(), NumericOutcome.PRECISION_CEILING)
    var at_root = dense.enclose([0, 0, 1], a, context=bc)
    assert_true(evaluation.check(at_root))
    assert_true(contains(at_root.value(), Rational(2)))
    var numerical = dense.asarray([0, 0, 1], context=bc)
    var numeric_at_root = dense.enclose(numerical, a, context=bc)
    assert_true(evaluation.check(numeric_at_root))
    assert_true(contains(numeric_at_root.value(), Rational(2)))
    var z = algebraic.ascomplex(a, 1)
    var cp = dense.asarray([z.copy(), algebraic.ascomplex(1)])
    var ce = dense.enclose(cp, 2, context=bc)
    assert_true(evaluation.check(ce))
    assert_true(contains(ce.value().imag(), Rational(1)))
    var sp = sparse.asarray((((1, 0), a.copy()), ((0, 1), 1)))
    var se = sparse.enclose(sp, (a.copy(), 1), context=bc)
    assert_true(evaluation.check(se))
    assert_true(contains(se.value(), Rational(3)))
    assert_true(evaluation.check(sparse.adaptive_evaluate(sp, (a.copy(), 1), options=options)))
    var zero_axes = sparse.asarray([((), a.copy())])
    assert_true(evaluation.check(sparse.enclose(zero_axes, (), context=bc)))
    # Replacing a valid refinement by the opposite root must fail before the
    # arithmetic trace can establish a claim about the wrong coefficient.
    var corrupted = dense.enclose(p, 2, context=bc)
    corrupted._data[].algebraic.value().witnesses[0].real = dense.root_between([-2, 0, 1], -2, -1).witness().copy()
    assert_false(evaluation.check(corrupted))
    corrupted = dense.enclose(p, 2, context=bc)
    corrupted._data[].algebraic = None
    assert_false(evaluation.check(corrupted))
    var wrong_point = dense.enclose([0, 0, 1], a, context=bc)
    wrong_point._data[].algebraic_points[0].source = algebraic.ascomplex(3)
    assert_false(evaluation.check(wrong_point))
    print("array algebraic numerical evaluation checks passed")
