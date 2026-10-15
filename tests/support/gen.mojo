"""Deterministic polynomial generators with bounded terms and exponents."""
from polynomials import Budget, CoefficientRing, Monomial, MonomialOrder, Term, SparsePoly
from .rng import SplitMix64

def coefficient[C: CoefficientRing](ctx: C.Context, value: Int, mut work: Budget) raises -> C:
    var result = C.zero(ctx)
    var one = C.one(ctx)
    for _ in range(abs(value)):
        result = result.add(one, work)
    if value < 0:
        result = result.neg(work)
    return result^

def polynomial[C: CoefficientRing](ctx: C.Context, order: MonomialOrder, mut rng: SplitMix64, mut work: Budget, terms: Int = 5, variables: Int = 3, degree: Int = 3, coefficient_bound: Int = 5, coefficient_bits: Int = 0) raises -> SparsePoly[C]:
    var result = List[Term[C]]()
    for _ in range(terms):
        var pairs = List[Tuple[Int64, UInt32]]()
        for i in range(variables):
            var exponent = UInt32(rng.below(degree + 1))
            if exponent:
                pairs.append((Int64(i), exponent))
        var c = coefficient[C](ctx, rng.below(2 * coefficient_bound + 1) - coefficient_bound, work)
        if coefficient_bits > 0:
            c = C.one(ctx)
            for _ in range(coefficient_bits - 1):
                c = c.add(c, work)
                if rng.below(2):
                    c = c.add(C.one(ctx), work)
            if rng.below(2):
                c = c.neg(work)
        result.append(Term[C](Monomial(pairs^, work), c^))
    return SparsePoly[C].from_terms(result^, order, ctx, work)


def planted[C: CoefficientRing](ctx: C.Context, order: MonomialOrder, mut rng: SplitMix64, mut work: Budget) raises -> Tuple[SparsePoly[C], SparsePoly[C], SparsePoly[C]]:
    """Generate a known factor and two multiples for later gcd properties."""
    var factor = polynomial[C](ctx, order, rng, work, terms=3)
    var a = polynomial[C](ctx, order, rng, work, terms=4)
    var b = polynomial[C](ctx, order, rng, work, terms=4)
    return (factor.copy(), factor.mul(a, work), factor.mul(b, work))
