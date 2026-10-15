"""Independent definition checks for a claimed reduced Gröbner basis."""
from polynomials import Budget, CoefficientField, SparsePoly, Term, GroebnerBasis, divide_by_list, reduce
from .check import Checker


def verify_basis[C: CoefficientField](basis: GroebnerBasis[C], generators: List[SparsePoly[C]], mut c: Checker, mut work: Budget) raises:
    for i in range(len(basis.basis)):
        ref f = basis.basis[i]
        f.check_invariants(work)
        var lt = f.leading(work)
        c.check(lt.coefficient == C.one(f.ring), "monic basis element")
        if i:
            c.check(basis.order.compare(basis.basis[i - 1].leading(work).monomial, lt.monomial, work) > 0, "decreasing leading monomials")
        for j in range(len(basis.basis)):
            if i == j:
                continue
            var other = basis.basis[j].leading(work)
            for term in f.terms:
                c.check(not other.monomial.divides(term.monomial, work), "fully reduced basis")
        for j in range(i):
            ref g = basis.basis[j]
            var other = g.leading(work)
            var lcm = lt.monomial.lcm(other.monomial, work)
            var a = SparsePoly[C].from_terms([Term[C](lcm.div(lt.monomial, work), lt.coefficient.inverse(work))], basis.order, f.ring, work)
            var b = SparsePoly[C].from_terms([Term[C](lcm.div(other.monomial, work), other.coefficient.inverse(work))], basis.order, f.ring, work)
            var s = a.mul(f, work).sub(b.mul(g, work), work)
            c.check(divide_by_list(s, basis.basis, work).remainder.is_zero(), "Buchberger criterion for every pair")
    for g in generators:
        c.check(reduce(g, basis, work).remainder.is_zero(), "original generator in returned ideal")
