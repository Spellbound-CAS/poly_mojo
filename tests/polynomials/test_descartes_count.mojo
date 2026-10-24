"""Descartes root counts on intervals agree with independent Sturm counts.

For corpus polynomials and generated families with clustered, repeated,
rational and many real roots, the count from a complete Descartes isolation
on an open interval must equal the Sturm count there, including intervals
whose ends cut an isolating bracket, and the checker must accept the count
and reject a wrong one, a wrong polynomial and endpoints that are roots.
"""
from apn_mojo import Rational, Integer
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, PolyError, DescartesCountCertificate, descartes_count, check_descartes_count
from polynomials.roots import count_roots
from polynomials._root_ops import _sign_at
from support.corpus import read_cases, parse_polynomial
from support.rng import SplitMix64
from support.check import Checker


def poly(values: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var c = List[QQ]()
    for v in values:
        c.append(QQ(Rational(v)))
    return DensePoly[QQ].from_coefficients(0, NoContext(), c^, work)


def linear(numerator: Int, denominator: Int, mut work: Budget) raises -> DensePoly[QQ]:
    # denominator * x - numerator, the root numerator / denominator.
    return poly([-numerator, denominator], work)


def endpoints(f: DensePoly[QQ], mut work: Budget) raises -> List[QQ]:
    # Nonroot points around and between the roots: a coarse odd grid plus a
    # point a third of the way into every bracket and exact root neighbour.
    var probe = descartes_count(f, QQ(Rational(-1, 1000003)), QQ(Rational(1, 1000033)), work)
    ref isolation = probe.isolation
    var points = List[QQ]()
    var bound = Rational(Integer(1) << (isolation.bound_exponent + 1))
    points.append(QQ(-bound))
    points.append(QQ(bound.copy()))
    for i in range(-9, 10):
        points.append(QQ(Rational(Integer(2 * i + 1), Integer(7))))
    for bracket in isolation.certificate.brackets:
        if bracket.low != bracket.high:
            var third = (bracket.high.value - bracket.low.value) / Rational(3)
            points.append(QQ(bracket.low.value + third))
            points.append(QQ(bracket.high.value - third))
    for r in isolation.certificate.exact_roots:
        points.append(QQ(r.value - Rational(1, 1000)))
        points.append(QQ(r.value + Rational(1, 1000)))
    var nonroots = List[QQ]()
    for p in points:
        if _sign_at(f, p, work) != 0:
            nonroots.append(p.copy())
    return nonroots^


def agree(mut c: Checker, f: DensePoly[QQ], name: String, mut work: Budget, mut rng: SplitMix64) raises:
    var points = endpoints(f, work)
    var tried = 0
    var same = True
    var checked = True
    var forged = False
    for _ in range(24):
        var i = rng.below(len(points))
        var j = rng.below(len(points))
        if points[i].value >= points[j].value:
            continue
        var result = descartes_count(f, points[i], points[j], work)
        var sturm = count_roots(f, points[i], points[j], work)
        same = same and result.count == sturm.count
        checked = checked and check_descartes_count(f, result, work)
        var wrong = result.copy()
        wrong.count += 1
        forged = forged or check_descartes_count(f, wrong, work)
        tried += 1
    c.check(tried > 0, name + ": intervals sampled")
    c.check(same, name + ": Descartes counts equal Sturm counts")
    c.check(checked, name + ": the checker accepts every count")
    c.check(not forged, name + ": the checker rejects a count one too high")


def main() raises:
    var c = Checker("descartes count")
    var work = Budget()
    var rng = SplitMix64(20261022)
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
            if f.degree() < 1 or f.degree() > 60:
                continue
            agree(c, f, record.name, work, rng)
            compared += 1
    c.check(compared >= 20, "the corpus supplies univariate polynomials")
    for trial in range(30):
        var n = 3 + rng.below(40)
        var bits = 1 + rng.below(30)
        var values = List[Int]()
        for _ in range(n + 1):
            values.append(rng.below(2 * (1 << bits) + 1) - (1 << bits))
        if values[n] == 0:
            values[n] = 1
        agree(c, poly(values, work), "random " + String(trial), work, rng)
    # Clustered rational roots, repeated factors and a Mignotte close pair.
    var cluster = poly([1], work)
    for i in range(1, 9):
        cluster = cluster.mul(linear(1000 * i + 1, 999 * i, work), work)
    agree(c, cluster, "rational cluster", work, rng)
    var repeated = poly([-2, 0, 1], work).mul(poly([-2, 0, 1], work), work).mul(poly([-1, 1], work), work).mul(poly([3, 0, 1], work), work)
    agree(c, repeated, "repeated", work, rng)
    var mignotte = poly([-2, 400, -20000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], work)
    agree(c, mignotte, "mignotte 20", work, rng)
    # The checker binds the count to its polynomial and nonroot endpoints.
    var f = poly([-6, 11, -6, 1], work)
    var a = QQ(Rational(1, 2))
    var b = QQ(Rational(5, 2))
    var result = descartes_count(f, a, b, work)
    c.check(result.count == 2, "x^3 - 6x^2 + 11x - 6 has two roots in (1/2, 5/2)")
    c.check(not check_descartes_count(poly([-6, 11, -6, 2], work), result, work), "the checker rejects another polynomial")
    var reversed = DescartesCountCertificate(result.isolation.copy(), b.copy(), a.copy(), 2)
    c.check(not check_descartes_count(f, reversed, work), "the checker rejects reversed endpoints")
    var at_root = DescartesCountCertificate(result.isolation.copy(), QQ(Rational(1)), b.copy(), 1)
    c.check(not check_descartes_count(f, at_root, work), "the checker rejects a root endpoint")
    var refused = False
    try:
        _ = descartes_count(f, QQ(Rational(2)), b, work)
    except e:
        refused = True
    c.check(refused, "a root endpoint raises")
    c.finish()
