"""Find exact roots and refine a certified numerical cover."""
from apn_mojo import Rational
from polynomials import NumericOptions
from poly import dense, roots, algebraic
from std.testing import assert_equal, assert_true


def main() raises:
    var p = [-2, 0, 1]
    var exact = dense.real_roots(p)
    assert_true(roots.check(exact))
    var positive = exact[1]
    var square = algebraic.multiply(positive, positive)
    assert_true(algebraic.check(square))
    assert_equal(algebraic.compare(square, 2), 0)
    print("Distinct real roots:", len(exact))
    print("Positive root squared:", algebraic.bounds(square)[0])

    var coarse = dense.real_root_cover(p, options=NumericOptions(absolute_bits=2), complete_only=True)
    var fine = roots.refine(coarse, NumericOptions(absolute_bits=8), complete_only=True)
    assert_true(roots.check(fine))
    for interval in roots.intervals(fine):
        assert_true(interval[1] - interval[0] <= Rational(1, 256))
    print("Refined intervals:", len(roots.intervals(fine)))
    print("Complete and accurate:", fine.complete() and fine.accurate())
