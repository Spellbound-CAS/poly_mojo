"""Refine a checked real-root result to a tighter absolute width."""
from std.testing import assert_equal, assert_true
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, NumericOptions, real_roots_exact, refine_real_roots_exact, check_real_roots_exact, descartes_intervals


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var polynomial = DensePoly[QQ].from_sparse(ring.parse("x^2 - 2", work), work, variable=Int64(0))
    var coarse = NumericOptions(absolute_bits=2, max_subdivisions=64)
    var fine = NumericOptions(absolute_bits=8, max_subdivisions=256)
    var initial = real_roots_exact(polynomial, coarse, work, complete_only=True)
    var refined = refine_real_roots_exact(polynomial, initial, coarse, fine, work, complete_only=True)
    assert_true(check_real_roots_exact(polynomial, refined, fine, work))
    assert_true(refined.complete)
    assert_true(refined.accurate)
    assert_true(refined.isolated)
    var count = len(descartes_intervals(refined.certificate, work))
    assert_equal(count, 2)
    print("Distinct real roots:", count)
    print("Complete cover:", refined.complete)
    print("Requested width reached:", refined.accurate)
