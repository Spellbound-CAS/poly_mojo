"""Descartes isolation: counts against Sturm, refined brackets, and certificate rejection."""
from apn_mojo import Integer, Rational
from polynomials import Budget, RootInterval, QQ, NoContext, DensePoly, NumericOptions, NumericOutcome, count_roots
from polynomials.descartes import _descartes_roots, _descartes_refine, _check_descartes_roots, descartes_intervals, DescartesRootResult, DescartesLeaf
from polynomials._root_ops import _fixed
from support.rng import SplitMix64
from support.check import Checker


def poly(values: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var c = List[QQ]()
    for v in values:
        c.append(QQ(Rational(v)))
    return DensePoly[QQ].from_coefficients(0, NoContext(), c^, work)


def linear(root: Rational, mut work: Budget) raises -> DensePoly[QQ]:
    return DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(-root), QQ.one(NoContext())], work)


def sturm_total(f: DensePoly[QQ], result: DescartesRootResult, mut work: Budget) raises -> Int:
    return count_roots(f, result.certificate.low, result.certificate.high, work).count


def agree(mut c: Checker, f: DensePoly[QQ], name: String, mut work: Budget) raises:
    # Root count against Sturm, one Sturm root per open bracket, ordered and
    # disjoint answers, a checked certificate, and refinement inside it.
    var options = NumericOptions(absolute_bits=1, max_subdivisions=1 << 30)
    var result = _descartes_roots(f, options, work)
    var roots = descartes_intervals(result.certificate, work)
    c.check(result.outcome == NumericOutcome.SUCCESS and result.complete and result.accurate, name + ": success")
    c.check(len(roots) == sturm_total(f, result, work), name + ": count")
    for i in range(len(roots)):
        if roots[i].low != roots[i].high:
            c.check(count_roots(f, roots[i].low, roots[i].high, work).count == 1, name + ": one root per bracket")
        if i:
            c.check(roots[i - 1].high.value <= roots[i].low.value and (roots[i - 1].low != roots[i].low or roots[i - 1].high != roots[i].high), name + ": ordered")
    c.check(_check_descartes_roots(f, result, options, work), name + ": certificate")
    var fine = NumericOptions(absolute_bits=60, max_subdivisions=1 << 30)
    var refined = _descartes_refine(result, fine, work)
    c.check(_check_descartes_roots(f, refined, fine, work), name + ": refined certificate")
    var tight = descartes_intervals(refined.certificate, work)
    var bound = Rational(Integer(1), Integer(1) << 60)
    for i in range(len(tight)):
        c.check(tight[i].high.value - tight[i].low.value <= bound, name + ": refined width")
        c.check(roots[i].low.value <= tight[i].low.value and tight[i].high.value <= roots[i].high.value, name + ": refined inside")


def rejections(mut c: Checker, mut work: Budget) raises:
    # (x - 1/2)(x^2 - 2)(x + 3)^2: an exact dyadic root, two irrational ones and a repeated one.
    var f = linear(Rational(1, 2), work).mul(poly([-2, 0, 1], work), work).mul(poly([3, 1], work).pow(2, work), work)
    var options = NumericOptions(absolute_bits=8, max_subdivisions=1 << 30)
    var good = _descartes_roots(f, options, work)
    c.check(_check_descartes_roots(f, good, options, work), "valid certificate")
    c.check(len(descartes_intervals(good.certificate, work)) == 4, "four distinct roots")
    var bad = good.copy()
    bad.certificate.leaves[0].variations += 1
    c.check(not _check_descartes_roots(f, bad, options, work), "wrong leaf bound")
    bad = good.copy()
    _ = bad.certificate.leaves.pop(1)
    c.check(not _check_descartes_roots(f, bad, options, work), "gap in the leaves")
    bad = good.copy()
    if len(bad.certificate.exact_roots):
        _ = bad.certificate.exact_roots.pop()
        c.check(not _check_descartes_roots(f, bad, options, work), "missing exact root")
    bad = good.copy()
    bad.certificate.exact_roots.append(QQ(Rational(7)))
    c.check(not _check_descartes_roots(f, bad, options, work), "extra exact root")
    bad = good.copy()
    bad.bound_exponent += 1
    c.check(not _check_descartes_roots(f, bad, options, work), "bound does not match the leaves")
    bad = good.copy()
    bad.squarefree = bad.squarefree.mul(bad.cofactor, work)
    c.check(not _check_descartes_roots(f, bad, options, work), "wrong squarefree part")
    bad = good.copy()
    var open_bracket = 0
    while bad.certificate.brackets[open_bracket].low == bad.certificate.brackets[open_bracket].high:
        open_bracket += 1
    bad.certificate.brackets[open_bracket].high = bad.certificate.brackets[open_bracket].low.copy()
    c.check(not _check_descartes_roots(f, bad, options, work), "bracket without its root")
    bad = good.copy()
    bad.certificate.brackets[0] = RootInterval(bad.certificate.leaves[0].low.copy(), bad.certificate.leaves[0].high.copy())
    c.check(not _check_descartes_roots(f, bad, options, work), "bracket outside its leaf")
    bad = good.copy()
    _ = bad.certificate.brackets.pop()
    c.check(not _check_descartes_roots(f, bad, options, work), "missing bracket")
    bad = good.copy()
    bad.accurate = False
    c.check(not _check_descartes_roots(f, bad, options, work), "false accuracy flag")
    var other = f.mul(poly([-5, 1], work), work)
    c.check(not _check_descartes_roots(other, good, options, work), "another polynomial")
    # Splitting a leaf in two at a nonroot keeps a valid certificate only with recomputed bounds.
    bad = good.copy()
    var split = False
    for i in range(len(bad.certificate.leaves)):
        if bad.certificate.leaves[i].variations == 0:
            split = True
            var leaf = bad.certificate.leaves[i].copy()
            var middle = QQ((leaf.low.value + leaf.high.value) / Rational(2))
            bad.certificate.leaves[i] = DescartesLeaf(leaf.low.copy(), middle.copy(), 1, True)
            bad.certificate.leaves.insert(i + 1, DescartesLeaf(middle^, leaf.high.copy(), 0, True))
            break
    c.check(not split or not _check_descartes_roots(f, bad, options, work), "invented root in a split leaf")
    bad = good.copy()
    for i in range(len(bad.certificate.leaves)):
        if bad.certificate.leaves[i].variations >= 1:
            bad.certificate.leaves[i].resolved = False
            break
    c.check(not _check_descartes_roots(f, bad, options, work), "unresolved leaf beside brackets")


def pairs(mut c: Checker, mut work: Budget) raises:
    # Mignotte's close pair resolves at a leaf with bound two and a sign
    # change between the roots; that middle point is what the checker trusts.
    var n = 30
    var power = List[Int](length=n + 1, fill=0)
    power[n] = 1
    var f = poly(power, work).add(poly([1, -100], work).pow(2, work).mul(poly([-2], work), work), work)
    var options = NumericOptions(absolute_bits=8, max_subdivisions=1 << 30)
    var good = _descartes_roots(f, options, work)
    c.check(_check_descartes_roots(f, good, options, work), "pair certificate")
    var leaf = -1
    var first = 0
    var seen = 0
    for i in range(len(good.certificate.leaves)):
        ref l = good.certificate.leaves[i]
        if l.variations == 2 and l.resolved and leaf < 0:
            leaf = i
            first = seen
        if l.resolved:
            seen += l.variations
    c.check(leaf >= 0, "a leaf with two roots")
    if leaf < 0:
        return
    var bad = good.copy()
    bad.certificate.leaves[leaf].resolved = False
    c.check(not _check_descartes_roots(f, bad, options, work), "unresolved pair with its brackets")
    bad = good.copy()
    _ = bad.certificate.brackets.pop(first)
    bad.certificate.leaves[leaf].variations = 1
    c.check(not _check_descartes_roots(f, bad, options, work), "pair counted as one root")
    bad = good.copy()
    ref l = bad.certificate.leaves[leaf]
    var outside = QQ((l.low.value * Rational(3) + l.high.value) / Rational(4))
    bad.certificate.brackets[first].high = outside.copy()
    bad.certificate.brackets[first + 1].low = outside^
    c.check(not _check_descartes_roots(f, bad, options, work), "middle point without a sign change")
    bad = good.copy()
    bad.certificate.brackets[first + 1] = bad.certificate.brackets[first].copy()
    c.check(not _check_descartes_roots(f, bad, options, work), "the same bracket twice")


def grids(mut c: Checker, mut work: Budget) raises:
    # Wilkinson's well separated roots show all their sign changes on a grid,
    # so one leaf with a large bound resolves many roots.
    var f = poly([1], work)
    for i in range(1, 25):
        f = f.mul(poly([-2 * i - 1, 2], work), work)
    var options = NumericOptions(absolute_bits=4, max_subdivisions=1 << 30)
    var good = _descartes_roots(f, options, work)
    c.check(_check_descartes_roots(f, good, options, work), "grid certificate")
    var leaf = -1
    var first = 0
    var seen = 0
    for i in range(len(good.certificate.leaves)):
        ref l = good.certificate.leaves[i]
        if l.variations >= 3 and l.resolved and leaf < 0:
            leaf = i
            first = seen
        if l.resolved:
            seen += l.variations
    c.check(leaf >= 0, "a leaf with at least three roots")
    if leaf < 0:
        return
    var bad = good.copy()
    # Swallow the next bracket's root: the two brackets overlap.
    bad.certificate.brackets[first].high = bad.certificate.brackets[first + 1].high.copy()
    c.check(not _check_descartes_roots(f, bad, options, work), "overlapping brackets")
    bad = good.copy()
    bad.certificate.leaves[leaf].variations += 2
    c.check(not _check_descartes_roots(f, bad, options, work), "grid leaf with a larger bound")


def reference_fixed(c: List[Integer], p: Integer, shift: Int, scale: Int) raises -> Integer:
    # The Integer floor-rounded Horner that _fixed's limb path replaces.
    var n = len(c) - 1
    var t = c[n] << scale
    for i in range(n - 1, -1, -1):
        t = ((t * p) >> shift) + (c[i] << scale)
    return t^


def fixed_horner(mut c: Checker) raises:
    # The limb path returns the same integer as the Integer path: negative
    # points and values, shifts across limb boundaries, wide coefficients.
    var rng = SplitMix64(20261021)
    var same = True
    for trial in range(400):
        var n = 1 + rng.below(30)
        var bits = 1 + rng.below(300)
        var coefficients = List[Integer]()
        for _ in range(n + 1):
            var x = Integer(Int(rng.next() >> 2))
            x = x << rng.below(bits)
            coefficients.append(-x if rng.below(2) == 1 else x^)
        var shift = [0, 1, 31, 63, 64, 65, 127, 128, 200][rng.below(9)]
        var p = Integer(Int(rng.next() >> 2)) << rng.below(66)
        if rng.below(2) == 1:
            p = -p
        var scale = rng.below(200)
        var got = _fixed(coefficients, p, Integer(1) << shift, shift, scale)
        same = same and got == reference_fixed(coefficients, p, shift, scale)
    c.check(same, "limb fixed-point Horner equals the Integer one")


def main() raises:
    var work = Budget()
    var c = Checker("descartes")
    fixed_horner(c)
    agree(c, poly([5], work), "constant", work)
    agree(c, poly([-3, 1], work), "integer root", work)
    agree(c, poly([0, -2, 0, 1], work), "root at zero", work)
    agree(c, poly([-2, 0, 1], work), "two irrational roots", work)
    agree(c, poly([1, 0, 1], work), "no real roots", work)
    var f = poly([1], work)
    for i in range(-4, 5):
        f = f.mul(linear(Rational(i, 4), work).pow(UInt64(abs(i) % 3 + 1), work), work)
    agree(c, f, "repeated dyadic roots", work)
    var g = poly([1], work)
    for i in range(1, 21):
        g = g.mul(linear(Rational(i), work), work)
    agree(c, g, "Wilkinson 20", work)
    # Mignotte: x^n - 2 (a x - 1)^2 has two roots within about a^(-n/2) of 1/a.
    for n in [8, 15, 24, 60]:
        var m = poly([1, -100], work).pow(2, work).mul(poly([-2], work), work)
        var power = List[Int](length=n + 1, fill=0)
        power[n] = 1
        agree(c, poly(power, work).add(m, work), "Mignotte " + String(n), work)
    # Chebyshev: only real roots, crowding at the ends of (-1, 1).
    var t0 = poly([1], work)
    var t1 = poly([0, 1], work)
    var x2 = poly([0, 2], work)
    for n in range(2, 81):
        var t2 = x2.mul(t1, work).sub(t0, work)
        t0 = t1^
        t1 = t2^
        if n == 40 or n == 80:
            agree(c, t1, "Chebyshev " + String(n), work)
    var rng = SplitMix64(20261020)
    for trial in range(12):
        var n = 5 + rng.below(40)
        var values = List[Int]()
        for _ in range(n + 1):
            values.append(rng.below(2001) - 1000)
        if values[n] == 0:
            values[n] = 1
        agree(c, poly(values, work), "random " + String(trial), work)
    var r = poly([3, -7, 2], work).mul(poly([-1, 0, 0, 1], work), work)
    agree(c, r.mul(r, work).mul(poly([-2, 0, 1], work), work), "repeated random factors", work)
    # A tiny root needs depth below the unit interval.
    agree(c, poly([-1, 1 << 40], work).mul(poly([1, 0, -3], work), work), "tiny root", work)
    rejections(c, work)
    pairs(c, work)
    grids(c, work)
    # A subdivision allowance of zero leaves the whole interval unresolved, and says so.
    var few = NumericOptions(absolute_bits=1, max_subdivisions=0)
    var partial = _descartes_roots(g, few, work)
    c.check(partial.outcome == NumericOutcome.SUBDIVISION_LIMIT and not partial.complete, "subdivision limit")
    c.check(_check_descartes_roots(g, partial, few, work), "partial certificate")
    c.finish()
