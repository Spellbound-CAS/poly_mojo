"""Canonical fractions, cancellation witnesses, and field identities."""
from polynomials import Budget, PolyErrorKind, QQ, NoContext, SparsePoly, MonomialOrder
from polynomials.fractions import RationalFunction
from polynomials.certificates import check_cancellation
from support.check import Checker
from support.corpus import read_cases, parse_polynomial, parse_order
from support.gen import polynomial
from support.rng import SplitMix64

def main() raises:
    var c = Checker("fraction")
    var work = Budget()
    var count = 0
    for record in read_cases():
        if record.field("op ") != "normalize_fraction":
            continue
        var order = MonomialOrder()
        for line in record.lines:
            if line.startswith("order "):
                order = parse_order(String(line[byte=6:]), work)
        var pt = record.field("in p ")
        var qt = record.field("in q ")
        var nt = record.field("out numerator ")
        var dt = record.field("out denominator ")
        var p = parse_polynomial(pt, order, work)
        var q = parse_polynomial(qt, order, work)
        var n = parse_polynomial(nt, order, work)
        var d = parse_polynomial(dt, order, work)
        var result = RationalFunction.normalize(p, q, work)
        c.check(result[0].numerator == n and result[0].denominator == d, record.name)
        c.check(check_cancellation(p, q, result[1], work), "cancellation identity")
        var corrupt = result[1].copy()
        corrupt.common = SparsePoly[QQ].zero(order, NoContext(), work)
        c.check(not check_cancellation(p, q, corrupt, work), "corrupt cancellation rejected")
        count += 1
    c.check(count == 7, "all seven fraction corpus cases")
    var order = MonomialOrder(MonomialOrder.LEX, [0, 1], work)
    var zero = SparsePoly[QQ].zero(order, NoContext(), work)
    var one = SparsePoly[QQ].constant(QQ.one(NoContext()), order, work)
    var rfzero = RationalFunction.normalize(zero, one, work)[0].copy()
    var rfone = RationalFunction.normalize(one, one, work)[0].copy()
    var rng = SplitMix64(1991)
    for _ in range(8):
        var p = polynomial[QQ](NoContext(), order, rng, work, terms=3, variables=2, degree=1)
        var q = polynomial[QQ](NoContext(), order, rng, work, terms=3, variables=2, degree=1)
        var r = polynomial[QQ](NoContext(), order, rng, work, terms=3, variables=2, degree=1)
        if p.is_zero() or q.is_zero() or r.is_zero():
            continue
        var a = RationalFunction.normalize(p, q, work)[0].copy()
        var b = RationalFunction.normalize(q, r, work)[0].copy()
        var d = RationalFunction.normalize(r, p, work)[0].copy()
        c.check(a == RationalFunction.normalize(p.mul(r, work), q.mul(r, work), work)[0], "equal fractions unique")
        c.check(a.add(rfzero, work) == a and a.mul(rfone, work) == a, "field identities")
        c.check(a.sub(a, work) == rfzero and a.add(a.neg(work), work) == rfzero, "additive inverse")
        c.check(a.div(a, work) == rfone, "multiplicative inverse")
        c.check(a.add(b, work) == b.add(a, work), "commutative sum")
        c.check(a.mul(b, work) == b.mul(a, work), "commutative product")
        c.check(a.mul(b, work).mul(d, work) == a.mul(b.mul(d, work), work), "associative product")
        c.check(a.add(b, work).add(d, work) == a.add(b.add(d, work), work), "associative sum")
        c.check(a.mul(b.add(d, work), work) == a.mul(b, work).add(a.mul(d, work), work), "distributivity")
    try:
        _ = RationalFunction.normalize(one, zero, work)
        c.check(False, "zero denominator accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero denominator DOMAIN")
    try:
        _ = rfone.div(rfzero, work)
        c.check(False, "division by zero fraction accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero fraction divisor DOMAIN")
    c.finish()
