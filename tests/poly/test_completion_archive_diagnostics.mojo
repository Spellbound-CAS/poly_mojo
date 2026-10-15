"""Root covers, multiplicities, diagnostics, and polishing retain evidence."""
from poly import dense, roots, io
from poly._archive_meta import _Archive
from polynomials import Budget, NumericOptions
from apn_mojo import Rational, ComplexContext, ArithmeticContext, FloatFormat
from std.testing import assert_true, assert_equal

def roundtrip[T: _Archive](value: T) raises:
    var text = io.dumps(value)
    var restored = io.loads[T](text)
    var work = Budget()
    assert_true(restored._archive_check(work))
    assert_equal(io.dumps(restored), text)

def main() raises:
    var p = [-1, 1]
    var cover = dense.real_root_cover(p, options=NumericOptions(absolute_bits=2))
    roundtrip(cover)
    roundtrip(roots.multiplicities(cover))
    roundtrip(dense.complex_root_cover(p, options=NumericOptions(absolute_bits=2)))
    var context = ArithmeticContext(format=FloatFormat(precision=80))
    roundtrip(dense.real_root_cover(dense.asarray(p, context=context), options=NumericOptions(absolute_bits=2)))
    roundtrip(roots.diagnostics(p, 0))
    roundtrip(roots.polish(p, 1, context=ComplexContext(context), radius=Rational(1, 16), options=NumericOptions(absolute_bits=2)))
    print("completion diagnostics archive checks passed")
