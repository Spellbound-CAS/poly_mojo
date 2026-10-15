"""Factorization, ideal membership, and certified real-root counts."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, SparsePoly, MonomialOrder, factor, squarefree, check_factorization, check_squarefree_product, groebner, membership, check_membership, count_roots, check_root_count, isolate_real_roots


def main() raises:
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-1)), QQ.zero(NoContext()), QQ.one(NoContext())], work)
    var factors = factor(f, work)
    var blocks = squarefree(f, work)
    if not check_factorization(f, factors, work) or not check_squarefree_product(f, blocks, work):
        raise Error("factor product did not check")
    print("Factors:", len(factors.factors))

    var order = MonomialOrder(MonomialOrder.LEX, [0, 1], work)
    var x = SparsePoly[QQ].variable(0, order, NoContext(), work)
    var y = SparsePoly[QQ].variable(1, order, NoContext(), work)
    var one = SparsePoly[QQ].constant(QQ.one(NoContext()), order, work)
    var generators: List[SparsePoly[QQ]] = [x.sub(y, work), x.pow(2, work).sub(one, work)]
    var basis = groebner(generators, order, work)
    var member = y.pow(2, work).sub(one, work)
    var proof = membership(member, generators, order, work)
    if not proof or not check_membership(member, generators, proof.value(), work):
        raise Error("membership identity did not check")
    print("Basis elements:", len(basis.basis))

    var low = QQ(Rational(-2))
    var high = QQ(Rational(2))
    var counted = count_roots(f, low, high, work)
    var roots = isolate_real_roots(f, work)
    if counted.count != 2 or not check_root_count(f, counted, work):
        raise Error("root count did not check")
    print("Exact real roots:", len(roots.exact_roots))
    print("All certificates checked")
