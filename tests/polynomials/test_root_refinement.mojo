"""Root refinement reuses certified cells and rejects stale or forged results."""
from apn_mojo import Rational, Ball
from polynomials import Budget, QQ, NoContext, DensePoly, NumericOptions, PolyErrorKind, real_roots_exact, complex_roots_exact, refine_real_roots_exact, refine_complex_roots_exact, check_real_roots_exact, check_complex_roots_exact, NumericDensePoly, BallDomain, real_roots_numeric, refine_real_roots_numeric, check_real_roots_numeric
from support.check import Checker


def main() raises:
    var c = Checker("root_refinement")
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], work)
    var coarse = NumericOptions(absolute_bits=1, max_subdivisions=2)
    var fine = NumericOptions(absolute_bits=7, max_subdivisions=256)
    var real = real_roots_exact(f, coarse, work)
    var refined = refine_real_roots_exact(f, real, coarse, fine, work, complete_only=True)
    c.check(refined.accurate and refined.isolated and check_real_roots_exact(f, refined, fine, work), "resume complete coarse Sturm cover")
    var retained = 0
    for cell in real.cells:
        if cell.count == 0:
            for after in refined.cells:
                if after.interval == cell.interval:
                    retained += 1
    c.check(retained > 0, "excluded cells reused")
    var coarse_complex = NumericOptions(absolute_bits=1, max_subdivisions=64, require_isolation=False)
    var initial = complex_roots_exact(f, coarse_complex, work)
    var fine_complex = NumericOptions(absolute_bits=4, max_subdivisions=1024)
    var refined_complex = refine_complex_roots_exact(f, initial, coarse_complex, fine_complex, work, complete_only=True)
    c.check(check_complex_roots_exact(f, refined_complex, fine_complex, work), "resume complex contour regions")
    var family = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(0, Rational(9, 10)), Ball(1)], work)
    var no_splits = NumericOptions(absolute_bits=1, max_subdivisions=0)
    var family_before = real_roots_numeric(family, no_splits, work)
    var family_after = refine_real_roots_numeric(family, family_before, no_splits, fine, work)
    c.check(family_after.complete and not family_after.accurate and check_real_roots_numeric(family, family_after, fine, work), "failed family cut retains certified parent")
    var uncertain = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(-2, Rational(1, 4096)), Ball(0), Ball(1)], work)
    var unresolved = real_roots_numeric(uncertain, no_splits, work)
    var resolved = refine_real_roots_numeric(uncertain, unresolved, no_splits, fine, work, complete_only=True)
    c.check(check_real_roots_numeric(uncertain, resolved, fine, work), "resume unresolved interval-family cover")
    var forged = real.copy()
    forged.cells[0].count += 1
    try:
        _ = refine_real_roots_exact(f, forged, coarse, fine, work)
        c.check(False, "forged previous root cover accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "previous root cover validated before reuse")
    var tiny = Budget(0)
    try:
        _ = refine_complex_roots_exact(f, initial, coarse_complex, fine_complex, tiny)
        c.check(False, "refinement ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "refinement budget preserves original")
    c.finish()
