"""Checked archives round-trip evidence rather than regenerating results."""
from poly import dense, sparse, exact, io, algebraic, coefficients, roots, rational, cad, modint
from poly._archive_meta import _Archive
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_true, assert_equal
from apn_mojo import BallContext

def roundtrip[T: _Archive](value: T) raises:
    var text = io.dumps(value)
    var restored = io.loads[T](text)
    var budget = Budget()
    assert_true(restored._archive_check(budget))
    assert_equal(io.dumps(restored), text)

def main() raises:
    roundtrip(dense.factor_checked([algebraic.asreal(1), algebraic.asreal(0), algebraic.asreal(1)]))
    roundtrip(dense.resultant_certified([algebraic.ascomplex(1), algebraic.ascomplex(0), algebraic.ascomplex(1)], [0, 1]))
    roundtrip(roots.diagnostics([algebraic.asreal(-1), algebraic.asreal(1)], 0, context=BallContext(precision=64)))
    print("completion algebraic archive checks passed")
