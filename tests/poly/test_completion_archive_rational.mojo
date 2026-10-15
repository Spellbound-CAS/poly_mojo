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
    var f = rational.asarray([1, 1], [0, 0, 1])
    roundtrip(rational.partial_fractions(f))
    roundtrip(rational.laurent(f, 0, terms=3))
    var a = algebraic.asreal(2)
    roundtrip(rational.partial_fractions(rational.asarray([a.copy()], [1, 0, 1])))
    print("completion rational archive checks passed")
