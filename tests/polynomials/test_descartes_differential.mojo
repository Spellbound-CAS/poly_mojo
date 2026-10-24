"""Descartes real roots agree with independent Sturm counts on the corpus and generated families.

For every univariate exact polynomial in the reference corpus, and for
families with clustered, repeated, rational and many real roots, both
Descartes outputs, real_roots_exact and isolate_real_roots, must agree with
the Sturm path: as many distinct real roots as the Sturm count over a Cauchy
bound, in increasing disjoint intervals that each hold exactly one root by a
Sturm count, exact roots that are roots, the same roots in the same order in
both outputs, and isolate_real_roots' rational roots found exactly.
"""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, NumericOptions, RootInterval, SturmSequence, real_roots_exact, check_real_roots_exact, isolate_real_roots, sturm_sequence, sign_changes, descartes_intervals, evaluate_dense, refine
from support.corpus import read_cases, parse_polynomial
from support.rng import SplitMix64
from support.check import Checker


def meets(a: RootInterval, b: RootInterval) raises -> Bool:
    # Closed intervals intersect.
    return a.low.value <= b.high.value and b.low.value <= a.high.value


def cauchy(f: DensePoly[QQ]) raises -> QQ:
    # Every root lies below 1 + max |a_i / a_n| in absolute value.
    var n = f.degree()
    var top = Rational(0)
    for i in range(n):
        var ratio = f.coefficients[i].value / f.coefficients[n].value
        if ratio < Rational(0):
            ratio = -ratio
        if ratio > top:
            top = ratio^
    return QQ(top + Rational(1))


def sturm_holds(f: DensePoly[QQ], sequence: SturmSequence, roots: List[RootInterval], mut work: Budget) raises -> Bool:
    # Increasing disjoint roots, each an exact root or an open interval whose
    # (low, high] Sturm count, with nonroot ends, is one.
    for i in range(len(roots)):
        ref root = roots[i]
        if i and not roots[i - 1].high.value <= root.low.value:
            return False
        if root.low == root.high:
            if not evaluate_dense(f, root.low, work).is_zero():
                return False
        elif not root.low.value < root.high.value or sign_changes(sequence, root.low, work) - sign_changes(sequence, root.high, work) != 1:
            return False
    return True


def agree(mut c: Checker, f: DensePoly[QQ], name: String, mut work: Budget, rationals: Int = -1) raises:
    if f.degree() < 1:
        return
    var sequence = sturm_sequence(f, work)
    var bound = cauchy(f)
    var total = sign_changes(sequence, bound.neg(work), work) - sign_changes(sequence, bound, work)
    var options = NumericOptions(absolute_bits=30, max_subdivisions=1 << 30)
    var result = real_roots_exact(f, options, work, complete_only=True)
    c.check(check_real_roots_exact(f, result, options, work), name + ": certificate")
    var descartes = descartes_intervals(result.certificate, work)
    c.check(len(descartes) == total, name + ": real_roots_exact finds the Sturm number of roots")
    c.check(sturm_holds(f, sequence, descartes, work), name + ": real_roots_exact intervals hold one root each")
    var isolation = isolate_real_roots(f, work)
    var merged = List[RootInterval]()
    var e = 0
    var k = 0
    var narrow = True
    while e < len(isolation.exact_roots) or k < len(isolation.intervals):
        if k >= len(isolation.intervals) or (e < len(isolation.exact_roots) and isolation.exact_roots[e].value < isolation.intervals[k].low.value):
            merged.append(RootInterval(isolation.exact_roots[e].copy(), isolation.exact_roots[e].copy()))
            e += 1
        else:
            narrow = narrow and (isolation.intervals[k].high.value - isolation.intervals[k].low.value) < Rational(1, 1 << 40)
            merged.append(isolation.intervals[k].copy())
            k += 1
    c.check(len(merged) == total, name + ": isolate_real_roots finds the Sturm number of roots")
    c.check(sturm_holds(f, sequence, merged, work) and narrow, name + ": isolate_real_roots intervals hold one root each, narrower than 2^-40")
    var refined = List[RootInterval]()
    var tight = True
    var target = QQ(Rational(1, 1 << 60) * Rational(1, 1 << 60))
    for root in merged:
        if root.low == root.high:
            refined.append(root.copy())
        else:
            var r = refine(f, root, target, work)
            tight = tight and r.high.value - r.low.value < target.value and meets(r, root)
            refined.append(r^)
    c.check(tight and sturm_holds(f, sequence, refined, work), name + ": refine narrows every interval below 2^-120 around its root")
    if rationals >= 0:
        c.check(len(isolation.exact_roots) == rationals, name + ": every rational root exact")
    if len(merged) == len(descartes):
        var same = True
        for i in range(len(merged)):
            same = same and meets(merged[i], descartes[i])
        c.check(same, name + ": both outputs list the same roots in order")


def poly(values: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var c = List[QQ]()
    for v in values:
        c.append(QQ(Rational(v)))
    return DensePoly[QQ].from_coefficients(0, NoContext(), c^, work)


def main() raises:
    var c = Checker("descartes differential")
    var work = Budget()
    var compared = 0
    for record in read_cases():
        var exact = True
        for line in record.lines:
            if line.startswith("ring ") and line != "ring QQ" and line != "ring ZZ":
                exact = False
        if not exact:
            continue
        for line in record.lines:
            if not line.startswith("in "):
                continue
            var text = String(line[byte=3:])
            var space = text.find(" ")
            if space < 0:
                continue
            var body = String(text[byte=space + 1:])
            var f: DensePoly[QQ]
            try:
                f = DensePoly[QQ].from_sparse(parse_polynomial(body, MonomialOrder(), work), work)
            except:
                continue
            if f.degree() < 1 or f.degree() > 200:
                continue
            agree(c, f, record.name, work)
            compared += 1
    c.check(compared >= 20, "the corpus supplies univariate polynomials")
    var rng = SplitMix64(20261020)
    for trial in range(40):
        var n = 3 + rng.below(60)
        var bits = 1 + rng.below(40)
        var values = List[Int]()
        for _ in range(n + 1):
            values.append(rng.below(2 * (1 << bits) + 1) - (1 << bits))
        if values[n] == 0:
            values[n] = 1
        agree(c, poly(values, work), "random " + String(trial), work)
    # Monomials: the squarefree part is x, whose only root is zero.
    agree(c, poly([0, 1], work), "x", work, rationals=1)
    agree(c, poly([0, 0, 0, 7], work), "7x^3", work, rationals=1)
    for n in [10, 25, 40]:
        var m = poly([1, -100], work).pow(2, work).mul(poly([-2], work), work)
        var power = List[Int](length=n + 1, fill=0)
        power[n] = 1
        agree(c, poly(power, work).add(m, work), "Mignotte " + String(n), work)
    var w = poly([1], work)
    for i in range(1, 31):
        w = w.mul(poly([-i, 1], work), work)
    agree(c, w, "Wilkinson 30", work, rationals=30)
    var t0 = poly([1], work)
    var t1 = poly([0, 1], work)
    for _ in range(2, 61):
        var t2 = poly([0, 2], work).mul(t1, work).sub(t0, work)
        t0 = t1^
        t1 = t2^
    agree(c, t1, "Chebyshev 60", work, rationals=0)
    agree(c, t1.mul(t1, work).mul(poly([-1, 0, 4], work), work), "Chebyshev 60 squared times 4x^2 - 1", work, rationals=2)
    # Roots i / 3 with multiplicities, many of them off the dyadic grid.
    var repeated = poly([1], work)
    for i in range(-6, 7):
        repeated = repeated.mul(poly([-i, 3], work).pow(UInt64(abs(i) % 3 + 1), work), work)
    agree(c, repeated, "repeated rational roots", work, rationals=13)
    # (5x - 2)(7x + 3)(x^2 - 2)(x^2 + 1): two non-dyadic rational roots
    # beside two irrational ones.
    var mixed = poly([-2, 5], work).mul(poly([3, 7], work), work).mul(poly([-2, 0, 1], work), work).mul(poly([1, 0, 1], work), work)
    agree(c, mixed, "rational beside irrational roots", work, rationals=2)
    # One root in (0, 3) beside a complex pair near 2 +- i/10: the Descartes
    # bound on the whole interval overcounts, so refine has to subdivide.
    var near = poly([-1, 1], work).mul(poly([401, -400, 100], work), work)
    var wide = RootInterval(QQ(Rational(0)), QQ(Rational(3)))
    var r = refine(near, wide, QQ(Rational(1, 1000)), work)
    c.check(r.low.value < Rational(1) and Rational(1) < r.high.value and r.high.value - r.low.value < Rational(1, 1000), "refine subdivides an overcounting interval")
    var two = RootInterval(QQ(Rational(0)), QQ(Rational(5, 2)))
    var refused = False
    try:
        _ = refine(poly([-2, 3, -1], work), two, QQ(Rational(1, 1000)), work)
    except:
        refused = True
    c.check(refused, "refine rejects an interval with two roots")
    c.finish()
