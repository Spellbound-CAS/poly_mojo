"""Gaussian rational inputs stay exact through enclosure and root certification."""
from apn_mojo import Rational, ExactComplex, BallContext
from polynomials import Budget, PolyErrorKind, QQ, QQi, NoContext, DensePoly, MonomialOrder, NumericBox, ComplexBallDomain, BallDomain, numeric_from_dense, numeric_from_sparse, numeric_bounds, NumericOptions, complex_roots_exact, check_complex_roots_exact, refine_complex_roots_exact, ComplexRootStreamChecker, complex_roots_exact_to
from std.testing import assert_true


def main() raises:
    var work = Budget()
    var c = ExactComplex(Rational(1, 3), Rational(-1, 7))
    var p = DensePoly[QQi].from_coefficients(0, NoContext(), [QQi(c.copy()), QQi(ExactComplex(1))], work)
    var exact_coefficient = NumericBox.point(QQ(Rational(1, 3)), QQ(Rational(-1, 7)), work)
    var enclosure = numeric_from_dense[ComplexBallDomain](p, BallContext(16), work)
    assert_true(numeric_bounds[ComplexBallDomain](enclosure.coefficients[0], work).contains(exact_coefficient, work))
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0)], work)
    var sparse = p.to_sparse(order, work)
    var enclosed_sparse = numeric_from_sparse[ComplexBallDomain](sparse, BallContext(16), work)
    assert_true(numeric_bounds[ComplexBallDomain](enclosed_sparse.terms[0].coefficient, work).contains(exact_coefficient, work))
    var context = BallContext(16)
    var failed = False
    try:
        _ = numeric_from_dense[BallDomain](p, context, work)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var options = NumericOptions(absolute_bits=2, max_subdivisions=256)
    var cover = complex_roots_exact(p, options, work, complete_only=True)
    assert_true(check_complex_roots_exact(p, cover, options, work))
    var exact_root = NumericBox.point(QQ(Rational(-1, 3)), QQ(Rational(1, 7)), work)
    assert_true(cover.regions[0].rectangle.contains(exact_root, work))
    var tighter = NumericOptions(absolute_bits=5, max_subdivisions=1024)
    var refined = refine_complex_roots_exact(p, cover, options, tighter, work, complete_only=True)
    assert_true(check_complex_roots_exact(p, refined, tighter, work))
    assert_true(refined.regions[0].rectangle.contains(exact_root, work))
    var checker = ComplexRootStreamChecker.from_dense(p, options, work)
    var receipt = complex_roots_exact_to(p, options, checker, work, complete_only=True)
    assert_true(checker.finish(receipt, work))
    print("Gaussian rational numerical checks passed")
