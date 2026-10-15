"""Independent gcd corpus, including integer contents and an unlucky prime."""
from apn_mojo import Integer
from polynomials import Budget, QQ, ZZ, ModInt, Modulus, NoContext, CoefficientRing, SparsePoly, DensePoly, MonomialOrder, prime, PolyErrorKind
from polynomials._prs import _prs
from support.gen import coefficient
from support.rng import SplitMix64
from polynomials._gcd import _gcd_value, _modular_primitive, _primitive_zz
from support.corpus import read_cases, parse_order, parse_polynomial, as_zz
from support.check import Checker
from polynomials.gcd import gcd
from polynomials.certificates import check_common_divisor, check_bezout, check_coprimality

def generated[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    var work = Budget()
    var rng = SplitMix64(941)
    var order = MonomialOrder(MonomialOrder.GREVLEX, [0, 1, 2], work)
    var one = SparsePoly[C].constant(C.one(ctx), order, work)
    var x = SparsePoly[C].variable(0, order, ctx, work)
    var y = SparsePoly[C].variable(1, order, ctx, work)
    var z = SparsePoly[C].variable(2, order, ctx, work)
    var a = x.scale(coefficient[C](ctx, 3, work), work)
    var b = x.add(one, work).scale(coefficient[C](ctx, 4, work), work)
    for i in range(16):
        var expected = one.copy()
        for _ in range(1 + i % 3):
            var shift = SparsePoly[C].constant(coefficient[C](ctx, rng.below(13) - 6, work), order, work)
            var linear = x.add(shift, work)
            if i % 3 > 0:
                linear = linear.add(y, work)
            if i % 3 > 1:
                linear = linear.add(z, work)
            expected = expected.mul(linear, work)
        comptime if C == ZZ:
            expected = expected.scale(coefficient[C](ctx, 2, work), work)
        var f = expected.mul(a, work)
        var h = expected.mul(b, work)
        var got = gcd(f, h, work, certify_maximal=i < 4)
        c.check(got.gcd == expected and check_common_divisor(f, h, got, work), "generated planted gcd and cofactors")
        if got.bezout:
            c.check(check_bezout(f, h, got, work), "generated univariate maximality")
        elif got.coprimality:
            c.check(check_coprimality(got.cofactors[0], got.cofactors[1], got.coprimality.value(), work), "generated content maximality")
    var tiny = Budget(0)
    try:
        _ = gcd(x, one, tiny)
        c.check(False, "gcd ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "gcd budget exhaustion")


def modular_cases(mut c: Checker) raises:
    var work = Budget()
    var order = MonomialOrder()
    var x = SparsePoly[ZZ].variable(0, order, NoContext(), work)
    var one = SparsePoly[ZZ].constant(ZZ.one(NoContext()), order, work)
    var p0 = Integer(Int(prime(0, work)))
    var p1 = Integer(Int(prime(1, work)))
    var height = p0 + Integer(123)
    var factor = x.scale(ZZ(height), work).add(one, work)
    var a = x.copy()
    var b = x.add(SparsePoly[ZZ].constant(ZZ(p1), order, work), work)
    # The first prime has the correct image degree but cannot reconstruct
    # height. At the second prime a and b coincide, so its larger-degree
    # gcd image must be discarded. A later lucky image finishes the CRT.
    c.check(_modular_primitive(factor.mul(a, work), factor.mul(b, work), 0, work) == factor, "discard later unlucky image")
    var drop = x.scale(ZZ(p0), work).add(one, work)
    c.check(_modular_primitive(factor.mul(drop, work), factor.mul(x.add(one, work), work), 0, work) == factor, "skip prime dropping one input degree")
    var high = x.pow(12, work).add(SparsePoly[ZZ].constant(ZZ(height), order, work), work)
    c.check(gcd(high.mul(x, work), high.mul(x.add(one, work), work), work).gcd == high, "public modular dispatch and multiple-prime reconstruction")
    var rng = SplitMix64(913)
    for degree in range(2, 14):
        var aa = List[ZZ]()
        var bb = List[ZZ]()
        for _ in range(degree):
            aa.append(ZZ(Integer(rng.below(11) - 5)))
            bb.append(ZZ(Integer(rng.below(11) - 5)))
        aa.append(ZZ.one(NoContext()))
        bb.append(ZZ.one(NoContext()))
        var f = DensePoly[ZZ].from_coefficients(0, NoContext(), aa^, work).to_sparse(order, work)
        var h = DensePoly[ZZ].from_coefficients(0, NoContext(), bb^, work).to_sparse(order, work)
        var sequence = _prs(f, h, 0, False, work)
        var expected = _primitive_zz(sequence[len(sequence) - 1].polynomial, work)
        c.check(_modular_primitive(f, h, 0, work) == expected, "modular and PRS gcd agree")


def main() raises:
    var c = Checker("gcd")
    var work = Budget()
    var count = 0
    for record in read_cases():
        if record.field("op ") != "gcd":
            continue
        var order = MonomialOrder()
        for line in record.lines:
            if line.startswith("order "):
                order = parse_order(String(line[byte=6:]), work)
        var ft = record.field("in f ")
        var ht = record.field("in h ")
        var et = record.field("out gcd ")
        var f = parse_polynomial(ft, order, work)
        var h = parse_polynomial(ht, order, work)
        var expected = parse_polynomial(et, order, work)
        if record.field("ring ") == "ZZ":
            var fz = as_zz(f, work)
            var hz = as_zz(h, work)
            var ez = as_zz(expected, work)
            var got = gcd(fz, hz, work, certify_maximal=True)
            c.check(got.gcd == ez and check_common_divisor(fz, hz, got, work), record.name)
            if got.coprimality:
                c.check(check_coprimality(got.cofactors[0], got.cofactors[1], got.coprimality.value(), work), "integer maximality")
            if record.name == "gcd.uni.unlucky_first_prime":
                var fp = _primitive_zz(fz, work)
                var hp = _primitive_zz(hz, work)
                var ep = _primitive_zz(ez, work)
                c.check(_modular_primitive(fp, hp, 0, work) == ep, "unlucky first prime restart")
        else:
            var got = gcd(f, h, work, certify_maximal=True)
            c.check(got.gcd == expected and check_common_divisor(f, h, got, work), record.name)
            if got.bezout:
                c.check(check_bezout(f, h, got, work), "field Bezout maximality")
            elif got.coprimality:
                c.check(check_coprimality(got.cofactors[0], got.cofactors[1], got.coprimality.value(), work), "field recursive maximality")
        count += 1
    c.check(count == 21, "all 21 gcd corpus cases")
    generated[ZZ](NoContext(), c)
    generated[QQ](NoContext(), c)
    generated[ModInt](Modulus(101, work), c)
    modular_cases(c)
    c.finish()
