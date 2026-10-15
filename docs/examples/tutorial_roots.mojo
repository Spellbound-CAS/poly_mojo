"""Count real roots and distinguish rational roots from isolating intervals."""
from std.testing import assert_equal, assert_true
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, PolynomialRing, DensePoly, count_roots, check_root_count, isolate_real_roots


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var polynomial = DensePoly[QQ].from_sparse(ring.parse("(x^2 - 2)*(3*x - 1)", work), work, variable=Int64(0))
    var count = count_roots(polynomial, QQ(Rational(-2)), QQ(Rational(2)), work)
    assert_equal(count.count, 3)
    assert_true(check_root_count(polynomial, count, work))
    var roots = isolate_real_roots(polynomial, work)
    assert_equal(len(roots.exact_roots), 1)
    assert_equal(roots.exact_roots[0].value, Rational(1, 3))
    assert_equal(len(roots.intervals), 2)
    print("Distinct real roots:", count.count)
    print("Rational root:", roots.exact_roots[0].value)
    print("Irrational-root intervals:", len(roots.intervals))
