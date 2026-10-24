"""Sparse ring laws and every independent arithmetic vector."""
from apn_mojo import Integer, Rational
from polynomials import Budget, PolyErrorKind, CoefficientRing, QQ, ZZ, ModInt, Modulus, NoContext, Monomial, MonomialOrder, Term, SparsePoly
from support.check import Checker
from support.rng import SplitMix64
from support.gen import polynomial
from support.kernel_vectors import read_lines, parse_poly
from polynomials.sparse import _sorted_indices, _ExactTermOrder
from polynomials._ordering import _stable_indices

def laws[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    var rng = SplitMix64(101)
    var work = Budget()
    for kind in range(4):
        var order = MonomialOrder()
        if kind:
            order = MonomialOrder(UInt8(kind), [2, 0, 1], work)
        for _ in range(12):
            var a = polynomial[C](ctx, order, rng, work)
            var b = polynomial[C](ctx, order, rng, work)
            var d = polynomial[C](ctx, order, rng, work)
            var zero = SparsePoly[C].zero(order, ctx, work)
            var one = SparsePoly[C].constant(C.one(ctx), order, work)
            var sum = a.add(b, work)
            var product = a.mul(b, work)
            var negative = a.neg(work)
            sum.check_invariants(work)
            product.check_invariants(work)
            negative.check_invariants(work)
            c.check(a.add(zero, work) == a and a.mul(one, work) == a, "identities")
            c.check(a.sub(a, work) == zero and a.add(negative, work) == zero, "inverse")
            c.check(sum == b.add(a, work), "sum commutes")
            c.check(product == b.mul(a, work), "product commutes")
            c.check(sum.add(d, work) == a.add(b.add(d, work), work), "sum associates")
            c.check(product.mul(d, work) == a.mul(b.mul(d, work), work), "product associates")
            c.check(a.mul(b.add(d, work), work) == product.add(a.mul(d, work), work), "distributivity")
            c.check(a.pow(3, work) == a.mul(a, work).mul(a, work), "repeated power")
            c.check(zero.pow(0, work) == one, "zero to zero keeps ring")
            c.check(a.scale(C.zero(ctx), work) == zero, "zero scale")
            c.check(a.reorder(MonomialOrder(), work).reorder(order, work) == a, "order round trip")
            var reversed = List[Term[C]]()
            for i in range(len(a.terms)):
                reversed.append(a.terms[len(a.terms) - 1 - i].copy())
            c.check(SparsePoly[C].from_terms(reversed^, order, ctx, work) == a, "deterministic normalization")

def kernel_keys(mut c: Checker) raises:
    """Packed KERNEL sort keys order terms exactly as the pairwise comparison does."""
    var work = Budget()
    var order = MonomialOrder()
    var rng = SplitMix64(733)
    var same = True
    for trial in range(60):
        var terms = List[Term[ZZ]]()
        var n = 1 + rng.below(40)
        for _ in range(n):
            var pairs = List[Tuple[Int64, UInt32]]()
            # Cases 0..39 fit keys; 40..49 have a fifth pair; 50..59 a wide exponent.
            var count = rng.below(5) if trial < 40 else (5 if trial < 50 else 1 + rng.below(3))
            var variable = Int64(0)
            for i in range(count):
                variable += Int64(1 + rng.below(3))
                var e = UInt32(1 + rng.below(4))
                if trial >= 50 and i == 0:
                    e = UInt32(0x10000 + rng.below(3))
                pairs.append((variable, e))
            # Repeat monomials often, so equal keys test stability.
            if rng.below(3) == 0 and len(terms) > 0:
                terms.append(Term[ZZ](terms[rng.below(len(terms))].monomial.copy(), ZZ(Integer(rng.below(9) - 4))))
            else:
                terms.append(Term[ZZ](Monomial(pairs^, work), ZZ(Integer(rng.below(9) - 4))))
        var got = _sorted_indices(terms, order, NoContext(), work)
        var expected = _stable_indices[_ExactTermOrder[ZZ]](terms, order, work)
        same = same and got == expected
    c.check(same, "packed KERNEL keys sort stably as the pairwise comparison")
    var large = List[Term[ZZ]]()
    large.append(Term[ZZ](Monomial([(Int64(0xFFFE), UInt32(1))], work), ZZ(Integer(1))))
    large.append(Term[ZZ](Monomial([(Int64(0xFFFF), UInt32(1))], work), ZZ(Integer(1))))
    large.append(Term[ZZ](Monomial([(Int64(3), UInt32(0xFFFF))], work), ZZ(Integer(1))))
    c.check(_sorted_indices(large, order, NoContext(), work) == _stable_indices[_ExactTermOrder[ZZ]](large, order, work), "variables and exponents at the key limits")


def main() raises:
    var c = Checker("sparse")
    var work = Budget()
    laws[ZZ](NoContext(), c)
    laws[QQ](NoContext(), c)
    laws[ModInt](Modulus(101, work), c)
    kernel_keys(c)
    var count = 0
    var counts = List[Int](length=4, fill=0)
    for line in read_lines("tests/reference/kernel/poly_ops.txt"):
        var f = line.split(" ")
        var op = String(f[0])
        if op != "add" and op != "mul" and op != "neg" and op != "pow":
            continue
        var a = parse_poly(String(f[1]), work)
        var result: SparsePoly[QQ]
        var expected: SparsePoly[QQ]
        if op == "neg":
            result = a.neg(work)
            expected = parse_poly(String(f[2]), work)
            counts[2] += 1
        elif op == "pow":
            result = a.pow(UInt64(Int(f[2])), work)
            expected = parse_poly(String(f[3]), work)
            counts[3] += 1
        else:
            var b = parse_poly(String(f[2]), work)
            result = a.add(b, work) if op == "add" else a.mul(b, work)
            expected = parse_poly(String(f[3]), work)
            counts[0 if op == "add" else 1] += 1
        result.check_invariants(work)
        c.check(result == expected, String("kernel ", op, " vector ", count))
        count += 1
    print("kernel arithmetic vectors:", count, counts)
    c.check(count == 240, "all arithmetic vectors executed")
    var order = MonomialOrder()
    var x = SparsePoly[QQ].variable(0, order, NoContext(), work)
    var y = SparsePoly[QQ].variable(1, order, NoContext(), work)
    var p = x.add(y, work)
    var repeated: List[Term[QQ]] = [Term[QQ](x.terms[0].monomial.copy(), QQ(Rational(2))), Term[QQ](Monomial(), QQ(Rational(3))), Term[QQ](x.terms[0].monomial.copy(), QQ(Rational(-2))), Term[QQ](x.terms[0].monomial.copy(), QQ(Rational(5))), Term[QQ](y.terms[0].monomial.copy(), QQ.zero(NoContext()))]
    var saved_terms = repeated.copy()
    var normalized = SparsePoly[QQ].from_terms(repeated.copy(), order, NoContext(), work)
    var expected = x.scale(QQ(Rational(5)), work).add(SparsePoly[QQ].constant(QQ(Rational(3)), order, work), work)
    c.check(normalized == expected, "duplicates cancel and reappear during owned normalization")
    normalized.check_invariants(work)
    normalized.terms[1].monomial.pairs[0] = (Int64(99), UInt32(1))
    c.check(repeated == saved_terms, "moving an owned copy does not alias the caller's terms")
    for steps in range(50):
        var limited = Budget(UInt64(steps))
        try:
            _ = SparsePoly[QQ].from_terms(repeated.copy(), order, NoContext(), limited)
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "normalization budget failure")
        c.check(repeated == saved_terms, "failed normalization preserves caller's copy")
    c.check(p.total_degree(work) == 1 and p.degree_in(0, work) == 1, "degrees")
    c.check(len(p.variables(work)) == 2, "variables")
    var original = p.copy()
    for steps in range(25):
        var small = Budget(UInt64(steps))
        try:
            _ = p.mul(p, small)
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "partial multiplication exhaustion")
        c.check(p == original, "input unchanged on failure")
    try:
        _ = p.add(p.reorder(MonomialOrder(MonomialOrder.LEX, [0, 1], work), work), work)
        c.check(False, "different orders accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "order DOMAIN")
    try:
        _ = SparsePoly[QQ].zero(order, NoContext(), work).leading(work)
        c.check(False, "zero leading accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero leading DOMAIN")
    c.finish()
