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
    print("completion roots archive checks passed")
