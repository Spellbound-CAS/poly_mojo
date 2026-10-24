"""All Gröbner corpus bases, pair criteria, and generator-order independence."""
from polynomials import Budget, QQ, ModInt, Modulus, CoefficientField, NoContext, SparsePoly, MonomialOrder, groebner, MembershipCertificate, check_membership
from support.check import Checker
from support.corpus import read_cases, parse_polynomial, parse_order, as_mod
from support.groebner import verify_basis


def compare_case[C: CoefficientField](generators: List[SparsePoly[C]], expected: List[SparsePoly[C]], order: MonomialOrder, name: String, mut c: Checker, mut work: Budget) raises:
    var got = groebner(generators, order, work)
    c.check(got.basis == expected, name)
    c.check(not got.cofactors, "cofactor tracking disabled")
    verify_basis(got, generators, c, work)
    var reversed_generators = List[SparsePoly[C]]()
    for i in range(len(generators) - 1, -1, -1):
        reversed_generators.append(generators[i].copy())
    var reverse = groebner(reversed_generators, order, work)
    c.check(reverse.basis == expected, "generator order independence")
    # The trace: the same basis, a replay that reproduces it, and expanded
    # rows that check as membership certificates; a changed basis fails.
    var traced = groebner(generators, order, work, trace=True)
    c.check(traced.basis == expected and Bool(traced.trace) and not traced.cofactors, "trace preserves basis")
    ref trace = traced.trace.value()
    c.check(trace.replay(generators, traced.basis, work), "trace replays to the basis")
    if len(expected) <= 3:
        for i in range(len(traced.basis)):
            var proof = MembershipCertificate[C](trace.cofactors(i, generators, work))
            c.check(check_membership(traced.basis[i], generators, proof, work), "expanded trace row identity")
    if len(traced.basis) > 0 and len(traced.basis[0].terms) > 0:
        var altered = traced.basis.copy()
        altered[0] = altered[0].add(altered[0], work)
        c.check(not trace.replay(generators, altered, work), "trace rejects a changed basis element")
    if len(expected) <= 3:
        var tracked = groebner(generators, order, work, track_cofactors=True)
        c.check(tracked.basis == expected and Bool(tracked.cofactors), "tracking preserves basis")
        for i in range(len(tracked.basis)):
            var proof = MembershipCertificate[C](tracked.cofactors.value()[i].copy())
            c.check(check_membership(tracked.basis[i], generators, proof, work), "basis cofactor identity")


def main() raises:
    var work = Budget()
    var c = Checker("groebner")
    var count = 0
    for record in read_cases():
        if record.field("op ") != "groebner":
            continue
        var order_text = record.field("order ")
        var order = parse_order(order_text, work)
        var generators = List[SparsePoly[QQ]]()
        var expected = List[SparsePoly[QQ]]()
        for line in record.lines:
            if line.startswith("in generator "):
                var text = String(line[byte=13:])
                generators.append(parse_polynomial(text, order, work))
            elif line.startswith("out basis "):
                var text = String(line[byte=10:])
                expected.append(parse_polynomial(text, order, work))
        var ring = record.field("ring ")
        if ring == "QQ":
            compare_case(generators, expected, order, record.name, c, work)
        else:
            var ctx = Modulus(UInt64(Int(ring[byte=3:])), work)
            var gm = List[SparsePoly[ModInt]]()
            var em = List[SparsePoly[ModInt]]()
            for p in generators:
                gm.append(as_mod(p, ctx, work))
            for p in expected:
                em.append(as_mod(p, ctx, work))
            compare_case(gm, em, order, record.name, c, work)
        count += 1
    c.check(count == 24, "all twenty-four Groebner corpus cases")
    # The reference's old inter-reduction discarded both copies of LM x.
    for kind in range(1, 4):
        var order = MonomialOrder(UInt8(kind), [0, 1], work)
        var x = SparsePoly[QQ].variable(0, order, NoContext(), work)
        var y = SparsePoly[QQ].variable(1, order, NoContext(), work)
        var inputs: List[SparsePoly[QQ]] = [x.copy(), x.add(y, work), y.copy(), x.copy()]
        var expected: List[SparsePoly[QQ]] = [x.copy(), y.copy()]
        compare_case(inputs, expected, order, "equal-leading-monomial regression", c, work)
    c.finish()
