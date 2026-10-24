"""Checked archives round-trip evidence rather than regenerating results."""
from poly import dense, sparse, exact, io, algebraic, coefficients, roots, rational, cad, modint
from poly._archive_meta import _Archive
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_true, assert_equal

def roundtrip[T: _Archive](value: T) raises:
    var text = io.dumps(value)
    var restored = io.loads[T](text)
    var budget = Budget()
    assert_true(restored._archive_check(budget))
    assert_equal(io.dumps(restored), text)

def main() raises:
    roundtrip(dense.real_roots([-2, 0, 1]))
    roundtrip(dense.count_roots([-2, 0, 1], -2, 2))
    roundtrip(dense.count_roots([algebraic.asreal(-2), algebraic.asreal(0), algebraic.asreal(1)], -2, 2))
    roundtrip(dense.complex_roots([1, 0, 1]))
    # Revision 11: counts carry Descartes evidence under the CHECKED-2 tag, and
    # a revision-10 archive of the same count, with Sturm evidence, fails cleanly.
    var text = io.dumps(dense.count_roots([-2, 0, 1], -2, 2))
    assert_true(text.find("\nCHECKED-2:poly.roots.RootCount:") >= 0)
    var sturm_era = "POLYNOMIALS-VALUE 1\nCHECKED-1:poly.roots.RootCount:axes=-1\n0\n3\n-2\n0\n1\n3\n0\n3\n-2\n0\n1\n0\n2\n0\n2\n0\n1\n2\n2\n0\n2\n0\n1/2\n0\n2\n0\n1\n-2\n2\nInt\n2\n"
    var kind = PolyErrorKind.LIMIT
    try:
        _ = io.loads[roots.RootCount](sturm_era)
    except e:
        kind = e.kind
    assert_true(kind == PolyErrorKind.ENCODING)
    print("completion roots archive checks passed")
