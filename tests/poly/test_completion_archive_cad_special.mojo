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
    var p = [(0, -1), (2, 1)]
    roundtrip(cad.cell_at([p.copy()], [0]))
    roundtrip(cad.cover([cad.negative(p), cad.nonnegative(p)]))
    roundtrip(cad.lazard_evaluate(p, List[Tuple[Int, Int]]()))
    print("completion cad_special archive checks passed")
