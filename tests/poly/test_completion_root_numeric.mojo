"""Algebraic interval queries and checked numerical diagnostics."""
from poly import dense, sparse, algebraic, roots, coefficients
from apn_mojo import BallContext, Rational, ComplexContext, ArithmeticContext, FloatFormat
from polynomials import Budget, PolyError, PolyErrorKind, NumericOptions
from std.testing import assert_equal, assert_true, assert_false

def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray(coefficients((algebraic.negative(a), 1)))
    var info = roots.diagnostics(p, Rational(3, 2), context=BallContext(precision=80))
    assert_true(roots.check(info))
    info._algebraic = None
    assert_false(roots.check(info))
    var context = ComplexContext(ArithmeticContext(format=FloatFormat(precision=80)))
    var polished = roots.polish(p, Rational(3, 2), context=context, radius=Rational(1, 16), options=NumericOptions(absolute_bits=2))
    assert_true(roots.check(polished))
    var cp = dense.asarray(coefficients((algebraic.ascomplex(algebraic.negative(a), -1), 1)))
    assert_true(roots.check(roots.diagnostics(cp, 1, context=BallContext(precision=80))))
    print("completion algebraic numerical root checks passed")
