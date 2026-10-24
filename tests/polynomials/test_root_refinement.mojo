"""Root refinement keeps certified leaves, agrees with Sturm counts, and rejects stale or forged results."""
from apn_mojo import Rational, Ball
from polynomials import Budget, QQ, NoContext, DensePoly, NumericOptions, PolyErrorKind, sturm_sequence, sign_changes, real_roots_exact, complex_roots_exact, refine_real_roots_exact, refine_complex_roots_exact, check_real_roots_exact, check_complex_roots_exact, NumericDensePoly, BallDomain, real_roots_numeric, refine_real_roots_numeric, check_real_roots_numeric
from polynomials import DescartesRootResult, descartes_intervals, evaluate_dense
from support.check import Checker


def sturm_agrees(f: DensePoly[QQ], result: DescartesRootResult, mut work: Budget) raises -> Bool:
    # Sturm counts, independent of the Descartes path: each leaf holds as many
    # roots as it has brackets (plus an exact root at its upper end, which
    # (low, high] counts), each open bracket holds one root, and each exact
    # root is a root. The chain is that of the squarefree part, so endpoints
    # at repeated roots of f are simple roots of it.
    var squarefree = List[QQ]()
    for value in result.squarefree.coefficients:
        squarefree.append(QQ(Rational(value.value)))
    var sequence = sturm_sequence(DensePoly[QQ].from_coefficients(0, NoContext(), squarefree^, work), work)
    ref certificate = result.certificate
    var next_root = 0
    var next = 0
    for leaf in certificate.leaves:
        var expected = leaf.variations if leaf.resolved else -1
        var at_high = next_root < len(certificate.exact_roots) and certificate.exact_roots[next_root] == leaf.high
        if at_high:
            next_root += 1
        var count = sign_changes(sequence, leaf.low, work) - sign_changes(sequence, leaf.high, work)
        if expected >= 0 and count != expected + (1 if at_high else 0):
            return False
        if leaf.resolved:
            next += leaf.variations
    for root in descartes_intervals(certificate, work):
        if root.low == root.high:
            if not evaluate_dense(f, root.low, work).is_zero():
                return False
        elif sign_changes(sequence, root.low, work) - sign_changes(sequence, root.high, work) != 1:
            return False
    return next == len(certificate.brackets)


def as_list(f: DensePoly[QQ]) -> List[QQ]:
    var values = List[QQ]()
    for v in f.coefficients:
        values.append(v.copy())
    return values^


def deep_refinement(mut c: Checker) raises:
    # Every leaf and bracket count must be the Sturm count, including repeated
    # roots, dyadic rational roots on the bisection grid, and nearly
    # coincident roots.
    var work = Budget()
    var cases = List[List[QQ]]()
    # (x - 1/2)^3 (x^2 - 2) (x + 3)
    cases.append([QQ(Rational(3, 4)), QQ(Rational(-17, 4)), QQ(Rational(57, 8)), QQ(Rational(-7, 8)), QQ(Rational(-23, 4)), QQ(Rational(3, 2)), QQ(Rational(1))])
    # (x - 1)(x - 1 - 2^-60) (x + 5/3)
    var close = Rational(1) + Rational(1, 1 << 60)
    var p0 = Rational(5, 3) * close
    var p1 = Rational(5, 3) * (Rational(-1) - close) + close
    var p2 = Rational(5, 3) - Rational(1) - close
    cases.append([QQ(p0^), QQ(p1^), QQ(p2^), QQ(Rational(1))])
    # prod_{i=1..12} (x - i): integer roots on the dyadic bisection grid.
    var wilkinson = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1))], work)
    for i in range(1, 13):
        wilkinson = wilkinson.mul(DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-i)), QQ(Rational(1))], work), work)
    cases.append(as_list(wilkinson))
    # T_9^2: nine double roots clustered towards +-1, and a root at 0.
    var t9 = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(0)), QQ(Rational(9)), QQ(Rational(0)), QQ(Rational(-120)), QQ(Rational(0)), QQ(Rational(432)), QQ(Rational(0)), QQ(Rational(-576)), QQ(Rational(0)), QQ(Rational(256))], work)
    cases.append(as_list(t9.mul(t9, work)))
    # x^30 - 2 (100 x - 1)^2: a pair of roots about 10^-32 apart near 1/100.
    var mignotte = List[QQ]()
    for i in range(31):
        mignotte.append(QQ(Rational(-2 if i == 0 else (400 if i == 1 else (-20000 if i == 2 else (1 if i == 30 else 0))))))
    cases.append(mignotte^)
    var distinct: List[Int] = [4, 3, 12, 9, 4]
    for k in range(len(cases)):
        var f = DensePoly[QQ].from_coefficients(0, NoContext(), cases[k].copy(), work)
        var coarse = NumericOptions(absolute_bits=1, max_subdivisions=1 << 20)
        var deep = NumericOptions(absolute_bits=200, max_subdivisions=1 << 20)
        var cover = real_roots_exact(f, coarse, work)
        var refined = refine_real_roots_exact(f, cover, coarse, deep, work, complete_only=True)
        c.check(check_real_roots_exact(f, refined, deep, work), String("deep refinement checks ", k))
        c.check(sturm_agrees(f, refined, work), String("deep refinement counts are Sturm counts ", k))
        c.check(len(descartes_intervals(refined.certificate, work)) == distinct[k], String("every distinct root covered ", k))
        c.check(sturm_agrees(f, cover, work), String("coarse counts are Sturm counts ", k))


def main() raises:
    var c = Checker("root_refinement")
    deep_refinement(c)
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], work)
    var coarse = NumericOptions(absolute_bits=1, max_subdivisions=2)
    var fine = NumericOptions(absolute_bits=7, max_subdivisions=256)
    var real = real_roots_exact(f, coarse, work)
    var refined = refine_real_roots_exact(f, real, coarse, fine, work, complete_only=True)
    c.check(refined.accurate and refined.isolated and check_real_roots_exact(f, refined, fine, work), "resume a complete coarse result")
    var retained = len(real.certificate.leaves) == len(refined.certificate.leaves)
    for i in range(len(real.certificate.leaves)):
        retained = retained and real.certificate.leaves[i].low == refined.certificate.leaves[i].low and real.certificate.leaves[i].high == refined.certificate.leaves[i].high
    c.check(retained, "leaves reused, only brackets tightened")
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
    forged.certificate.leaves[0].variations += 1
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
