"""Root covers, multiplicities, diagnostics, and polishing retain evidence."""
from poly import dense, roots, io
from poly._archive_meta import _Archive
from polynomials import Budget, NumericOptions, PolyError, PolyErrorKind
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
    # Archives written before real covers carried Descartes evidence hold a
    # Sturm cell cover; they are rejected with a typed error, never misread.
    var old_cover = open("tests/reference/real_cover_sturm_archive.txt", "r").read()
    var kind = PolyErrorKind.LIMIT
    try:
        _ = io.loads[type_of(cover)](old_cover)
    except e:
        kind = e.kind
    assert_true(kind == PolyErrorKind.ENCODING)
    var old_multiplicities = open("tests/reference/multiplicities_sturm_archive.txt", "r").read()
    kind = PolyErrorKind.LIMIT
    try:
        _ = io.loads[type_of(roots.multiplicities(cover))](old_multiplicities)
    except e:
        kind = e.kind
    assert_true(kind == PolyErrorKind.ENCODING)
    roundtrip(dense.complex_root_cover(p, options=NumericOptions(absolute_bits=2)))
    var context = ArithmeticContext(format=FloatFormat(precision=80))
    roundtrip(dense.real_root_cover(dense.asarray(p, context=context), options=NumericOptions(absolute_bits=2)))
    roundtrip(roots.diagnostics(p, 0))
    roundtrip(roots.polish(p, 1, context=ComplexContext(context), radius=Rational(1, 16), options=NumericOptions(absolute_bits=2)))
    print("completion diagnostics archive checks passed")
