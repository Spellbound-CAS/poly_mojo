"""Membership witnesses, nonmembers, empty ideals, and exact boundary errors."""
from polynomials import Budget, PolyErrorKind, QQ, ModInt, Modulus, CoefficientField, NoContext, SparsePoly, MonomialOrder, groebner, reduce, membership, MembershipCertificate, check_membership
from support.check import Checker
from support.corpus import read_cases, parse_polynomial, parse_order
from support.groebner import verify_basis
from support.gen import polynomial
from support.rng import SplitMix64


def generated[C: CoefficientField](ctx: C.Context, mut c: Checker) raises:
    var work = Budget()
    var rng = SplitMix64(50317)
    for kind in range(1, 4):
        var order = MonomialOrder(UInt8(kind), [0, 1], work)
        var x = SparsePoly[C].variable(0, order, ctx, work)
        var y = SparsePoly[C].variable(1, order, ctx, work)
        var one = SparsePoly[C].constant(C.one(ctx), order, work)
        var zero = SparsePoly[C].zero(order, ctx, work)
        var generators: List[SparsePoly[C]] = [x.pow(2, work).sub(y, work), x.mul(y, work).sub(one, work), zero.copy()]
        var basis = groebner(generators, order, work, track_cofactors=True)
        verify_basis(basis, generators, c, work)
        for _ in range(8):
            var a = polynomial[C](ctx, order, rng, work, variables=2, terms=3, degree=2)
            var b = polynomial[C](ctx, order, rng, work, variables=2, terms=3, degree=2)
            var f = a.mul(generators[0], work).add(b.mul(generators[1], work), work)
            var proof = membership(f, generators, order, work)
            c.check(Bool(proof) and check_membership(f, generators, proof.value(), work), "planted ideal member")
            var corrupt = proof.take()
            corrupt.cofactors[0] = corrupt.cofactors[0].add(one, work)
            c.check(not check_membership(f, generators, corrupt, work), "corrupt multiplier")
            _ = corrupt.cofactors.pop()
            c.check(not check_membership(f, generators, corrupt, work), "wrong multiplier count")
        c.check(not membership(one, generators, order, work), "one outside proper ideal")
        var empty = List[SparsePoly[C]]()
        var zero_proof = membership(zero, empty, order, work)
        c.check(Bool(zero_proof) and check_membership(zero, empty, zero_proof.value(), work), "zero in empty ideal")
        c.check(not membership(one, empty, order, work), "one outside empty ideal")
        var original = generators.copy()
        var tiny = Budget(100)
        try:
            _ = groebner(generators, order, tiny, track_cofactors=True)
            c.check(False, "Groebner budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT and generators == original, "basis budget preserves inputs")
        var kernel = MonomialOrder()
        try:
            _ = groebner(generators, kernel, work)
            c.check(False, "KERNEL basis accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "KERNEL basis rejected")


def main() raises:
    var work = Budget()
    var c = Checker("membership")
    var count = 0
    for record in read_cases():
        if record.field("op ") != "membership":
            continue
        var order_text = record.field("order ")
        var order = parse_order(order_text, work)
        var generators = List[SparsePoly[QQ]]()
        for line in record.lines:
            if line.startswith("in generator "):
                var text = String(line[byte=13:])
                generators.append(parse_polynomial(text, order, work))
        var text = record.field("in f ")
        var f = parse_polynomial(text, order, work)
        var expected = record.field("out member ") == "true"
        var proof = membership(f, generators, order, work)
        c.check(Bool(proof) == expected, record.name)
        if proof:
            c.check(check_membership(f, generators, proof.value(), work), "corpus membership identity")
        var basis = groebner(generators, order, work)
        var rem_text = record.field("out remainder ")
        c.check(reduce(f, basis, work).remainder == parse_polynomial(rem_text, order, work), "corpus normal form")
        count += 1
    c.check(count == 6, "all six membership corpus cases")
    generated[QQ](NoContext(), c)
    generated[ModInt](Modulus(7, work), c)
    generated[ModInt](Modulus(101, work), c)
    c.finish()
