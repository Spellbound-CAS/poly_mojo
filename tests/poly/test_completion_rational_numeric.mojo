"""Numerical fractions and certified quotient enclosures use explicit contexts."""
from poly import rational, algebraic, sparse, io
from apn_mojo import BallContext, Rational, ArithmeticContext, FloatFormat
from polynomials import PolyError, PolyErrorKind, Budget, BallDomain, numeric_bounds
from std.testing import assert_true, assert_false

def main() raises:
    var context = BallContext(precision=80)
    var exact = rational.asarray([1], [1, 1])
    var f = rational.asarray(exact, context=context)
    var result = rational.enclose(f, 1, context=context)
    assert_true(rational.check(result))
    assert_true(rational.check(rational.enclose(exact, 1, context=context)))
    assert_true(rational.check(rational.enclose(exact, (1,), context=context)))
    var restored = io.loads[type_of(f)](io.dumps(f))
    assert_true(rational.check(rational.enclose(restored, 1, context=context)))
    _ = rational.derivative(f, context=context)
    _ = rational.power(f, -2, context=context)
    _ = rational.substitute(f, [(0, sparse.asarray([(0, 1), (1, 1)], context=context))], context=context)
    _ = rational.compose(f, f, context=context)
    var core = rational.to_core(f)
    _ = rational.from_core[1](core)
    _ = rational.asarray(f, context=context)
    var saved_result = io.loads[type_of(result)](io.dumps(result))
    assert_true(rational.check(saved_result))
    var wrong = rational.enclose(f, 2, context=context)
    saved_result._denominator = wrong._denominator.copy()
    assert_false(rational.check(saved_result))
    var sum = rational.add(f, f, context=context)
    var work = Budget()
    var bounds = numeric_bounds[BallDomain](rational.evaluate(sum, 1, context=context), work)
    assert_true(bounds.real_low.value <= Rational(1) and bounds.real_high.value >= Rational(1))
    var failed = False
    try:
        _ = rational.enclose(f, -1, context=context)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var numerical = rational.asarray([1, 2], [1], context=ArithmeticContext(format=FloatFormat(precision=80)))
    _ = rational.evaluate(numerical, 1, context=ArithmeticContext(format=FloatFormat(precision=80)))
    print("completion numerical rational checks passed")
