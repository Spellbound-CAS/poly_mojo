"""Exact dependent-coordinate specialization and root completeness."""
from apn_mojo import Rational
from polynomials import Budget, QQ, DensePoly, MonomialOrder, RealAlgebraic, real_root, real_rational, algebraic_roots_at, check_algebraic_roots_at, algebraic_evaluate_sparse, check_algebraic_evaluate_sparse, algebraic_specialize, check_algebraic_specialization, PolyErrorKind
from support.corpus import parse_polynomial
from support.check import Checker


def main() raises:
    var work = Budget()
    var c = Checker("algebraic specialization")
    var order = MonomialOrder()
    var definition = DensePoly[QQ].from_sparse(parse_polynomial("x0^2 - 2", order, work), work, variable=Int64(0))
    var sqrt2 = real_root(definition, QQ(Rational(1)), QQ(Rational(2)), work)
    var point: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), sqrt2.copy()), (Int64(1), sqrt2.copy())]
    var dependent = parse_polynomial("x0*x2^2 - x1*x2^2 + x2 - x0", order, work)
    var proof = algebraic_roots_at(dependent, 2, point, work)
    c.check(check_algebraic_roots_at(dependent, 2, point, proof, work), "dependent coordinates and vanished leading coefficient")
    c.check(len(proof.isolation.roots) == 1 and proof.isolation.roots[0].polynomial == definition, "unique root is sqrt(2)")
    var bad = proof.copy()
    bad.specialization.monomials[0] = bad.specialization.monomials[1].copy()
    c.check(not check_algebraic_roots_at(dependent, 2, point, bad, work), "wrong coefficient binding rejected")
    var relation = parse_polynomial("x0*x1 - 2", order, work)
    var evaluated = algebraic_evaluate_sparse(relation, point, work)
    c.check(check_algebraic_evaluate_sparse(relation, point, evaluated, work) and evaluated.value.rational and evaluated.value.low.is_zero(), "dependent irrational coordinates evaluate exactly")
    var wrong = evaluated.copy()
    wrong.steps[0].operation = 4
    c.check(not check_algebraic_evaluate_sparse(relation, point, wrong, work), "wrong arithmetic schedule rejected")
    var partial = algebraic_specialize(dependent, point, work)
    c.check(check_algebraic_specialization(dependent, point, partial, work), "partial specialization includes zero coefficients")
    var zero = parse_polynomial("x0*x2 - x1*x2", order, work)
    var rejected = False
    try:
        _ = algebraic_roots_at(zero, 2, point, work)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    c.check(rejected, "zero specialization has no finite root list")
    c.finish()
