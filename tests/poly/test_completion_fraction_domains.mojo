"""Certified fractions bind complex, algebraic, and multivariate inputs."""
from poly import dense, sparse, rational, algebraic, io, coefficients, evaluation
from apn_mojo import Rational, ExactComplex, BallContext, contains
from std.testing import assert_true, assert_equal
from polynomials import Budget, ComplexBallDomain, numeric_bounds

def main() raises:
    var context = BallContext(precision=80)
    var f = rational.asarray([ExactComplex(1)], [ExactComplex(0, 1)])
    var result = rational.enclose(f, 0, context=context)
    assert_true(rational.check(result))
    var work = Budget()
    var box = numeric_bounds[ComplexBallDomain](result.value(), work)
    assert_equal(box.real_low.value, Rational(0))
    assert_equal(box.real_high.value, Rational(0))
    assert_equal(box.imag_low.value, Rational(-1))
    assert_equal(box.imag_high.value, Rational(-1))
    assert_true(rational.check(io.loads[type_of(result)](io.dumps(result))))
    var root = dense.root_between([-2, 0, 1], 1, 2)
    var aa = rational.asarray([root.copy()], [1])
    var algebraic_result = rational.enclose(aa, root, context=context)
    assert_true(rational.check(algebraic_result))
    assert_true(rational.check(io.loads[type_of(algebraic_result)](io.dumps(algebraic_result))))
    var polynomial_result = dense.enclose(coefficients((root.copy(), 1)), root, context=context)
    assert_true(evaluation.check(io.loads[type_of(polynomial_result)](io.dumps(polynomial_result))))
    var p = rational.from_terms([((1, 0), 1)], [((0, 0), 1), ((0, 1), 1)], context=context)
    assert_true(contains(rational.evaluate(p, (3, 2), context=context), Rational(1)))
    var bounded = rational.enclose(p, (3, 2), context=context)
    assert_true(rational.check(bounded))
    var restored = io.loads[type_of(p)](io.dumps(p))
    assert_equal(restored.ndim(), 2)
    print("completion fraction domain checks passed")
