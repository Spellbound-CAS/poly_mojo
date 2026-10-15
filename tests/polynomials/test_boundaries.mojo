"""Domain, capacity, normalization, and shared-budget regression tests."""
from apn_mojo import Integer, Rational
from polynomials import Budget, PolyErrorKind, ZZ, QQ, ModInt, Modulus, NoContext, Monomial, MonomialOrder, Term, SparsePoly, DensePoly, evaluate, evaluate_dense, derivative, substitute, to_kernel_bytes, from_kernel_bytes, MAX_KERNEL_BYTES
from support.check import Checker
from support.corpus import parse_polynomial
from support.gen import polynomial, planted
from support.rng import SplitMix64


def main() raises:
    var c = Checker("boundaries")
    var work = Budget()
    var kernel = MonomialOrder()
    var ctx = Modulus(101, work)
    var other_ctx = Modulus(103, work)
    var x = SparsePoly[ModInt].variable(0, kernel, ctx, work)
    var y = SparsePoly[ModInt].variable(0, kernel, other_ctx, work)
    var dx = DensePoly[ModInt].from_sparse(x, work)
    var dy = DensePoly[ModInt].from_sparse(y, work)
    for scenario in range(14):
        try:
            if scenario == 0:
                _ = x.add(y, work)
            elif scenario == 1:
                _ = x.sub(y, work)
            elif scenario == 2:
                _ = x.mul(y, work)
            elif scenario == 3:
                _ = x.scale(ModInt.one(other_ctx), work)
            elif scenario == 4:
                var terms = List[Term[ModInt]]()
                terms.append(Term[ModInt](Monomial(), ModInt.one(other_ctx)))
                _ = SparsePoly[ModInt].from_terms(terms^, kernel, ctx, work)
            elif scenario == 5:
                _ = dx.add(dy, work)
            elif scenario == 6:
                _ = dx.sub(dy, work)
            elif scenario == 7:
                _ = dx.mul(dy, work)
            elif scenario == 8:
                _ = dx.scale(ModInt.one(other_ctx), work)
            elif scenario == 9:
                _ = DensePoly[ModInt].from_coefficients(0, ctx, [ModInt.one(other_ctx)], work)
            elif scenario == 10:
                _ = substitute(x, 0, y, work)
            elif scenario == 11:
                _ = evaluate_dense(dx, ModInt.one(other_ctx), work)
            elif scenario == 12:
                var point: List[Tuple[Int64, ModInt]] = [(Int64(0), ModInt.one(other_ctx))]
                _ = evaluate(x, point, work)
            else:
                _ = DensePoly[ModInt].from_sparse(x, work, variable=Int64(2))
            c.check(False, String("accepted domain error ", scenario))
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, String("DOMAIN scenario ", scenario))
    var corrupted = x.copy()
    corrupted.terms[0].coefficient = ModInt.zero(ctx)
    try:
        corrupted.check_invariants(work)
        c.check(False, "zero sparse coefficient accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "invariant zero coefficient")
    corrupted = x.copy()
    var duplicate = corrupted.terms[0].copy()
    corrupted.terms.append(duplicate^)
    try:
        corrupted.check_invariants(work)
        c.check(False, "duplicate term accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "invariant duplicate term")
    var corrupt_dense = dx.copy()
    corrupt_dense.coefficients.append(ModInt.zero(ctx))
    try:
        corrupt_dense.check_invariants(work)
        c.check(False, "zero dense leading coefficient accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "invariant dense leading zero")
    var p = parse_polynomial("2*x0^2 - 3/7*x1 + 5", kernel, work)
    var point: List[Tuple[Int64, QQ]] = [(Int64(0), QQ(Rational(2))), (Int64(1), QQ(Rational(7)))]
    c.check(evaluate(p, point, work) == QQ(Rational(10)), "corpus parser independent evaluation")
    var zero = parse_polynomial("x0 - x0", kernel, work)
    c.check(zero.is_zero(), "normalization cancellation")
    var original = p.copy()
    var fraction = QQ(Rational(3, 7))
    for scenario in range(13):
        var empty = Budget(0)
        try:
            if scenario == 0:
                _ = p.neg(empty)
            elif scenario == 1:
                _ = p.leading(empty)
            elif scenario == 2:
                _ = p.total_degree(empty)
            elif scenario == 3:
                _ = p.degree_in(0, empty)
            elif scenario == 4:
                _ = p.variables(empty)
            elif scenario == 5:
                p.check_invariants(empty)
            elif scenario == 6:
                _ = p.reorder(kernel, empty)
            elif scenario == 7:
                _ = evaluate(p, point, empty)
            elif scenario == 8:
                _ = substitute(p, 0, p, empty)
            elif scenario == 9:
                _ = derivative(p, 0, empty)
            elif scenario == 10:
                _ = fraction.inverse(empty)
            elif scenario == 11:
                _ = p.scale(QQ(Rational(2)), empty)
            else:
                _ = to_kernel_bytes(p, empty)
            c.check(False, String("uncharged operation ", scenario))
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, String("LIMIT scenario ", scenario))
        c.check(p == original, "exhaustion preserves operands")
    var near_limit = Monomial()
    near_limit.pairs.append((Int64(0), UInt32.MAX))
    var huge = SparsePoly[ZZ].from_terms([Term[ZZ](near_limit.copy(), ZZ(Integer(1)))], kernel, NoContext(), work)
    var small = Budget(100)
    try:
        _ = DensePoly[ZZ].from_sparse(huge, small)
        c.check(False, "unbounded dense allocation")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "preflight huge dense allocation")
    var huge_coefficient = QQ(Rational(Integer(1) << (MAX_KERNEL_BYTES * 8)))
    var oversized = SparsePoly[QQ].constant(huge_coefficient, kernel, work)
    try:
        _ = to_kernel_bytes(oversized, work)
        c.check(False, "oversized encode accepted")
    except e:
        c.check(e.kind == PolyErrorKind.CAPACITY, "encoder ceiling")
    var rng = SplitMix64(7)
    var large = polynomial[ZZ](NoContext(), kernel, rng, work, terms=3, coefficient_bits=256)
    large.check_invariants(work)
    c.check(large.terms[0].coefficient.size() >= 32, "coefficient size generator")
    var factors = planted[ZZ](NoContext(), kernel, rng, work)
    factors[0].check_invariants(work)
    factors[1].check_invariants(work)
    factors[2].check_invariants(work)
    c.finish()
