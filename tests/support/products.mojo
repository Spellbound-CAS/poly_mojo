"""Geometric products with coefficients derived by counting exponent pairs."""
from polynomials import Budget, CoefficientRing, SparsePoly, Term, Monomial, MonomialOrder


def geometric[C: CoefficientRing](n: Int, variable: Int64, alternating: Bool, base: C, order: MonomialOrder, mut budget: Budget) raises -> SparsePoly[C]:
    var terms = List[Term[C]]()
    for i in range(n):
        var monomial = Monomial([(variable, UInt32(i))], budget) if i else Monomial()
        terms.append(Term[C](monomial^, base.neg(budget) if alternating and i % 2 else base.copy()))
    return SparsePoly[C].from_terms(terms^, order, base.context(), budget)


def expected_product[C: CoefficientRing](n: Int, distinct: Bool, alternating: Bool, base: C, order: MonomialOrder, mut budget: Budget) raises -> SparsePoly[C]:
    var terms = List[Term[C]]()
    var square = base.mul(base, budget)
    if distinct:
        for i in range(n):
            for j in range(n):
                var pairs = List[Tuple[Int64, UInt32]]()
                if i:
                    pairs.append((Int64(0), UInt32(i)))
                if j:
                    pairs.append((Int64(1), UInt32(j)))
                terms.append(Term[C](Monomial(pairs^, budget), square.copy()))
    else:
        for k in range(2 * n - 1):
            # For even n, (1+...+x^(n-1))(1-x+...-x^(n-1))
            # has +1 at even k<n, -1 at even k>=n, and zero at odd k.
            var coefficient = C.zero(base.context())
            if alternating:
                if k % 2:
                    continue
                coefficient = square.copy() if k < n else square.neg(budget)
            else:
                var count = min(k + 1, 2 * n - 1 - k)
                for _ in range(count):
                    coefficient = coefficient.add(square, budget)
            if not coefficient.is_zero():
                terms.append(Term[C](Monomial([(Int64(0), UInt32(k))], budget) if k else Monomial(), coefficient^))
    return SparsePoly[C].from_terms(terms^, order, base.context(), budget)
