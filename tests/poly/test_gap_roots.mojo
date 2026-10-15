"""Multiplicity, diagnostics, and validated polishing have distinct claims."""
from std.testing import assert_equal, assert_true
from apn_mojo import Rational, ArithmeticContext, FloatFormat, ComplexContext
from polynomials import NumericOptions
from poly import dense, roots


def main() raises:
    var options = NumericOptions(absolute_bits=2, max_subdivisions=256)
    var cover = dense.real_root_cover([1, -2, 1], options=options, complete_only=True)
    var repeated = roots.multiplicities(cover)
    assert_equal(len(repeated.counts()), 1)
    assert_equal(repeated.counts()[0], 2)
    assert_true(roots.check(repeated))
    var info = roots.diagnostics([-2, 0, 1], Rational(3, 2))
    assert_equal(info.evidence().residual.real_low.value, Rational(1, 4))
    assert_equal(info.evidence().derivative.real_low.value, Rational(3))
    assert_true(roots.check(info))
    var stationary = roots.diagnostics([-2, 0, 1], 0)
    assert_true(not stationary.evidence().newton_correction_upper)
    var context = ComplexContext(ArithmeticContext(format=FloatFormat(precision=80)))
    var polished = roots.polish([-2, 0, 1], Rational(3, 2), context=context, radius=Rational(1, 32), options=options)
    assert_true(roots.check(polished))
    var numeric = dense.asarray([-2, 0, 1], context=ArithmeticContext(format=FloatFormat(precision=80)))
    assert_true(roots.check(roots.diagnostics(numeric, Rational(3, 2))))
    assert_true(roots.check(roots.polish(numeric, Rational(3, 2), context=context, radius=Rational(1, 32), options=options)))
    print("array gap roots checks passed")
