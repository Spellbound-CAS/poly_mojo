"""Cross-phase validation, charged proof checking, and reduction regressions."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, ModInt, Modulus, NoContext, DensePoly, SparsePoly, MonomialOrder, MembershipCertificate, groebner, membership, check_membership, factor, check_factorization, squarefree, check_squarefree_product, count_roots, check_root_count, sturm_sequence, SturmSequence, RootCountCertificate, sign_changes, isolate_real_roots, refine, RootInterval
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var work = Budget()
    var c = Checker("p3_p5_boundaries")
    var order = MonomialOrder(MonomialOrder.LEX, [0, 1], work)
    var p7 = Modulus(7, work)
    var p11 = Modulus(11, work)
    var z7 = SparsePoly[ModInt].zero(order, p7, work)
    var z11 = SparsePoly[ModInt].zero(order, p11, work)
    var mixed: List[SparsePoly[ModInt]] = [z7.copy(), z11.copy()]
    try:
        _ = groebner(mixed, order, work)
        c.check(False, "mixed zero generators accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "mixed zero generators rejected")
    var zeros: List[SparsePoly[ModInt]] = [z7.copy()]
    try:
        _ = membership(z11, zeros, order, work)
        c.check(False, "mixed-field zero member accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "mixed-field zero member rejected")
    var x = SparsePoly[QQ].variable(0, order, NoContext(), work)
    var short_order = MonomialOrder(MonomialOrder.LEX, [1], work)
    var generators: List[SparsePoly[QQ]] = [x.copy()]
    try:
        _ = groebner(generators, short_order, work)
        c.check(False, "unknown order variable accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "unknown order variable rejected")
    # Proofs remain identities when the computational order differs from storage.
    var stored = x.reorder(MonomialOrder(), work)
    var stored_generators: List[SparsePoly[QQ]] = [stored.copy()]
    var member = stored.pow(3, work)
    var proof = membership(member, stored_generators, order, work)
    c.check(Bool(proof) and check_membership(member, stored_generators, proof.value(), work), "membership across explicit orders")
    var text = String("x0^4 - 2*x0^2 + 1")
    var f = DensePoly[QQ].from_sparse(parse_polynomial(text, order, work), work)
    var factors = factor(f, work)
    var blocks = squarefree(f, work)
    var a = QQ(Rational(-2))
    var b = QQ(Rational(2))
    var roots = count_roots(f, a, b, work)
    for operation in range(7):
        var tiny = Budget(0)
        try:
            if operation == 0:
                _ = check_factorization(f, factors, tiny)
            elif operation == 1:
                _ = check_squarefree_product(f, blocks, tiny)
            elif operation == 2:
                _ = check_membership(member, stored_generators, proof.value(), tiny)
            elif operation == 3:
                _ = check_root_count(f, roots, tiny)
            elif operation == 4:
                _ = sign_changes(roots.sequence, a, tiny)
            elif operation == 5:
                _ = membership(member, stored_generators, order, tiny)
            else:
                _ = refine(f, RootInterval(QQ(Rational(0)), b.copy()), QQ.one(NoContext()), tiny)
            c.check(False, "operation ignored empty budget")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "all new checkers and operations charge")
    var corrupt = roots.copy()
    corrupt.sequence.polynomials.clear()
    c.check(not check_root_count(f, corrupt, work), "empty Sturm sequence refused")
    corrupt = roots.copy()
    corrupt.sequence.quotients[len(corrupt.sequence.quotients) - 1] = DensePoly[QQ](0, NoContext())
    c.check(not check_root_count(f, corrupt, work), "corrupt terminal Sturm quotient refused")
    # T-C1 relies on the identities, not strictly decreasing degrees.
    var xd = DensePoly[QQ].from_sparse(x, work)
    var one = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], work)
    var alternate = SturmSequence([xd.copy(), one.copy(), one.copy()], [xd.add(one, work), one.copy()])
    var alternate_proof = RootCountCertificate(alternate^, a.copy(), b.copy(), 1)
    c.check(check_root_count(xd, alternate_proof, work), "valid generalized Sturm identity chain")
    # One-term multiplication must retain the KERNEL fallback: this order is
    # not invariant under multiplication, unlike all three admissible orders.
    var input_text = String("x0 + x1^2 + 3*x0*x1")
    var multiplier_text = String("2*x0*x1")
    var expected_text = String("2*x0^2*x1 + 2*x0*x1^3 + 6*x0^2*x1^2")
    for kind in range(4):
        var selected = MonomialOrder()
        if kind:
            selected = MonomialOrder(UInt8(kind), [0, 1], work)
        var left = parse_polynomial(input_text, selected, work)
        var right = parse_polynomial(multiplier_text, selected, work)
        var expected = parse_polynomial(expected_text, selected, work)
        var got = left.mul(right, work)
        got.check_invariants(work)
        c.check(got == expected and right.mul(left, work) == expected, "single-term multiplication in every order")
    c.finish()
