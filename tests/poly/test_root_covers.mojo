"""Root cover claims, functional refinement, and complex multiplicities."""
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, NumericOptions
from poly import dense, sparse, roots, coefficients
from std.testing import assert_equal, assert_true


def main() raises:
    var coarse = NumericOptions(absolute_bits=2, max_subdivisions=256)
    var fine = NumericOptions(absolute_bits=6, max_subdivisions=1024)
    var cover = dense.real_root_cover([-2, 0, 1], options=coarse, complete_only=True)
    assert_true(cover.complete() and cover.accurate() and cover.isolated())
    assert_true(roots.check(cover))
    assert_equal(len(roots.intervals(cover)), 2)
    var refined = roots.refine(cover, fine, complete_only=True)
    assert_true(roots.check(refined))
    assert_equal(cover.options().absolute_bits, 2)
    assert_equal(refined.options().absolute_bits, 6)
    for interval in roots.intervals(refined):
        assert_true(interval[1] - interval[0] <= Rational(1, 64))
        assert_equal(interval[2], 1)
    var complex_cover = dense.complex_root_cover([1, 0, 1], options=coarse, complete_only=True)
    assert_true(roots.check(complex_cover))
    assert_equal(len(roots.rectangles(complex_cover)), 2)
    var complex_refined = roots.refine(complex_cover, fine, complete_only=True)
    assert_true(roots.check(complex_refined))
    var cluster_options = NumericOptions(absolute_bits=2, max_subdivisions=256, require_isolation=False)
    var repeated = dense.complex_root_cover([1, -2, 1], options=cluster_options, complete_only=True)
    assert_true(roots.check(repeated))
    assert_equal(roots.rectangles(repeated)[0][4], 2)
    assert_true(not repeated.isolated())
    var linear = dense.complex_root_cover(coefficients((ExactComplex(Rational(1, 3), Rational(-1, 7)), 1)), options=coarse, complete_only=True)
    assert_true(roots.check(linear))
    var sp = [(0, 1), (2, 1)]
    assert_true(roots.check(sparse.complex_root_cover(sp, options=coarse, complete_only=True)))
    assert_true(roots.check(sparse.real_root_cover([(0, -2), (2, 1)], options=coarse, complete_only=True)))
    var limited = dense.real_root_cover([-2, 0, 1], options=NumericOptions(absolute_bits=64, max_subdivisions=0))
    assert_true(not limited.accurate())
    print("array root cover checks passed")
