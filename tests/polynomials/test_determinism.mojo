"""Replay every algorithm and its complete certificate after term permutations.

Writable's derived representation includes every stored field, recursively;
these transcripts are compared within a run, never used as golden answers.
Correctness and corrupted certificates are checked by the phase tests.
"""
from apn_mojo import Rational
from polynomials import *
from support.check import Checker
from support.gen import polynomial
from support.rng import SplitMix64


def permute[C: CoefficientRing](f: SparsePoly[C], mut work: Budget) raises -> SparsePoly[C]:
    var terms = List[Term[C]]()
    for i in range(len(f.terms) - 1, -1, -1):
        # Also split each coefficient across duplicate monomials and add zero.
        ref term = f.terms[i]
        terms.append(Term[C](term.monomial.copy(), term.coefficient.add(C.one(f.ring), work)))
        terms.append(Term[C](term.monomial.copy(), C.one(f.ring).neg(work)))
        terms.append(Term[C](term.monomial.copy(), C.zero(f.ring)))
    return SparsePoly[C].from_terms(terms^, f.order, f.ring, work)


def compare(a: List[String], b: List[String], mut c: Checker):
    c.check(len(a) == len(b), "replay transcript length")
    for i in range(min(len(a), len(b))):
        c.check(a[i] == b[i], String("replay result and all certificate fields ", i))


def sparse_run[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C]) raises -> List[String]:
    var work = Budget()
    var out = List[String]()
    out.append(String(f.add(g, work)))
    out.append(String(f.sub(g, work)))
    out.append(String(f.mul(g, work)))
    out.append(String(f.neg(work)))
    out.append(String(f.pow(3, work)))
    out.append(String(f.scale(C.one(f.ring), work)))
    out.append(String(derivative(f, 0, work)))
    out.append(String(substitute(f, 0, g, work)))
    var point: List[Tuple[Int64, C]] = [(Int64(0), C.one(f.ring)), (Int64(1), C.one(f.ring))]
    out.append(String(evaluate(f, point, work)))
    out.append(String(content(f, 0, work)))
    out.append(String(primitive_part(f, 0, work)))
    var product = f.mul(g, work)
    out.append(String(divide_exact(product, f, work)))
    out.append(String(gcd(product, f, work, certify_maximal=True)))
    out.append(String(coprimality_certificate(f, g, work)))
    out.append(String(resultant(f, g, 0, work)))
    out.append(String(subresultants(f, g, 0, work)))
    comptime if conforms_to(C, CoefficientField):
        if f.order.kind != MonomialOrder.KERNEL:
            var generators: List[SparsePoly[C]] = [f.copy(), g.copy()]
            out.append(String(divide_by_list(product, generators, work)))
            var basis = groebner(generators, f.order, work, track_cofactors=True)
            out.append(String(basis))
            out.append(String(reduce(product, basis, work)))
            out.append(String(membership(product, generators, f.order, work)))
    comptime if C == QQ:
        ref fq = rebind[SparsePoly[QQ]](f)
        ref gq = rebind[SparsePoly[QQ]](g)
        var normalized = RationalFunction.normalize(fq, gq, work)
        out.append(String(normalized))
        var a = normalized[0].copy()
        var b = RationalFunction.normalize(gq, fq, work)[0].copy()
        out.append(String(a.add(b, work)))
        out.append(String(a.sub(b, work)))
        out.append(String(a.mul(b, work)))
        out.append(String(a.div(b, work)))
        out.append(String(a.neg(work)))
        var encoded = to_kernel_bytes(fq, work)
        out.append(String(encoded))
        out.append(String(from_kernel_bytes(encoded, work)))
    return out^


def sparse_replays[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    var work = Budget()
    var rng = SplitMix64(6102026)
    for kind in range(4):
        var order = MonomialOrder()
        if kind:
            order = MonomialOrder(UInt8(kind), [0, 1], work)
        var x = SparsePoly[C].variable(0, order, ctx, work)
        var y = SparsePoly[C].variable(1, order, ctx, work)
        var one = SparsePoly[C].constant(C.one(ctx), order, work)
        # Shared coefficient content exercises all fields of the coprimality proof arena.
        var f = x.mul(y, work)
        var g = x.add(one, work).mul(y.add(one, work), work)
        var baseline = sparse_run(f, g)
        compare(baseline, sparse_run(f, g), c)
        compare(baseline, sparse_run(permute(f, work), permute(g, work)), c)
        for _ in range(3):
            var a = polynomial[C](ctx, order, rng, work, terms=4, variables=2, degree=2)
            var b = polynomial[C](ctx, order, rng, work, terms=4, variables=2, degree=2)
            var pa = permute(a, work)
            var pb = permute(b, work)
            c.check(pa == a and pb == b, "permuted duplicates normalize identically")
            c.check(a.mul(b, work) == pa.mul(pb, work), "generated arithmetic permutation")
            c.check(a.mul(b, work) == a.mul(b, work), "generated arithmetic repeat")


def dense_run[C: CoefficientRing](f: DensePoly[C], g: DensePoly[C]) raises -> List[String]:
    var work = Budget()
    var out = List[String]()
    out.append(String(f.add(g, work)))
    out.append(String(f.sub(g, work)))
    out.append(String(f.mul(g, work)))
    out.append(String(f.mul_schoolbook(g, work)))
    out.append(String(f.neg(work)))
    out.append(String(f.pow(3, work)))
    out.append(String(f.scale(C.one(f.ring), work)))
    out.append(String(evaluate_dense(f, C.one(f.ring), work)))
    out.append(String(pseudo_divmod(f, g, work)))
    comptime if conforms_to(C, CoefficientField):
        out.append(String(divmod(f, g, work)))
        out.append(String(extended_gcd(f, g, work)))
    comptime if C == ZZ:
        out.append(String(factor_integer(rebind[DensePoly[ZZ]](f), work)))
    comptime if C == QQ:
        ref fq = rebind[DensePoly[QQ]](f)
        out.append(String(squarefree(fq, work)))
        out.append(String(factor(fq, work)))
        out.append(String(sturm_sequence(fq, work)))
        var low = QQ(Rational(-10))
        var high = QQ(Rational(10))
        out.append(String(count_roots(fq, low, high, work)))
        out.append(String(sign_changes(sturm_sequence(fq, work), high, work)))
        out.append(String(sign_on_interval(fq, low, high, work)))
        var roots = isolate_real_roots(fq, work)
        out.append(String(roots))
        for interval in roots.intervals:
            out.append(String(refine(fq, interval, QQ(Rational(1, 100000000000000)), work)))
    comptime if C == ModInt:
        out.append(String(factor_mod_p(rebind[DensePoly[ModInt]](f), work, 123)))
    return out^


def dense_replays[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    var work = Budget()
    var order = MonomialOrder()
    var x = SparsePoly[C].variable(0, order, ctx, work)
    var one = SparsePoly[C].constant(C.one(ctx), order, work)
    var two = one.add(one, work)
    var f = x.sub(one, work).mul(x.pow(2, work).sub(two, work), work)
    var g = x.pow(2, work).add(one, work)
    var df = DensePoly[C].from_sparse(f, work)
    var dg = DensePoly[C].from_sparse(g, work)
    var baseline = dense_run(df, dg)
    compare(baseline, dense_run(df, dg), c)
    compare(baseline, dense_run(DensePoly[C].from_sparse(permute(f, work), work), DensePoly[C].from_sparse(permute(g, work), work)), c)
    comptime if C == ZZ or C == QQ:
        var a = df.pow(4, work).to_sparse(order, work)
        var b = dg.pow(4, work).to_sparse(order, work)
        var proof = String(gcd(a, b, work, certify_maximal=True))
        c.check(proof == String(gcd(a, b, work, certify_maximal=True)), "modular gcd and certificate repeat")
        c.check(proof == String(gcd(permute(a, work), permute(b, work), work, certify_maximal=True)), "modular gcd and certificate permutation")
    # Cross the dispatch threshold to exercise iterative Karatsuba as well.
    var large = df.pow(24, work)
    var rebuilt = DensePoly[C].from_sparse(permute(large.to_sparse(order, work), work), work)
    c.check(large.mul(large, work) == rebuilt.mul(rebuilt, work), "Karatsuba permutation")
    c.check(large.mul(large, work) == large.mul(large, work), "Karatsuba repeat")


def main() raises:
    var c = Checker("determinism")
    var work = Budget()
    var ctx = Modulus(101, work)
    sparse_replays[ZZ](NoContext(), c)
    sparse_replays[QQ](NoContext(), c)
    sparse_replays[ModInt](ctx, c)
    dense_replays[ZZ](NoContext(), c)
    dense_replays[QQ](NoContext(), c)
    dense_replays[ModInt](ctx, c)
    for i in range(12):
        c.check(prime(i, work) == prime(i, work), "cold/warm prime equality")
    c.check(crt([UInt64(1), 2, 4], [UInt64(2), 3, 5], work) == crt([UInt64(4), 1, 2], [UInt64(5), 2, 3], work), "CRT permutation")
    c.finish()
