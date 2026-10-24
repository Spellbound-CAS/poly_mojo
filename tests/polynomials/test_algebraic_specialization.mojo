"""Exact dependent-coordinate specialization and root completeness."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, RealAlgebraic, real_root, real_rational, algebraic_roots_at, check_algebraic_roots_at, algebraic_evaluate_sparse, check_algebraic_evaluate_sparse, algebraic_specialize, check_algebraic_specialization, PolyErrorKind
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
    # A valid reducible witness: x^3 - x^2 + x - 1 = (x - 1)(x^2 + 1) on (0, 2)
    # denotes 1, so x0 - 1 vanishes there although z(x - 1) is not 0 modulo it.
    var reducible = DensePoly[QQ].from_sparse(parse_polynomial("x0^3 - x0^2 + x0 - 1", order, work), work, variable=Int64(0))
    var one: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), RealAlgebraic(reducible^, QQ(Rational(0)), QQ(Rational(2)), False))]
    var shifted = parse_polynomial("x0 - 1", order, work)
    var at_one = algebraic_evaluate_sparse(shifted, one, work)
    c.check(at_one.value.rational and at_one.value.low.value == Rational(0) and check_algebraic_evaluate_sparse(shifted, one, at_one, work), "evaluation at a reducible witness of a rational value")
    var forged = at_one.copy()
    forged.value = real_rational(QQ(Rational(1, 2)), work)
    c.check(not check_algebraic_evaluate_sparse(shifted, one, forged, work), "a value that is not f(alpha) is rejected at a reducible witness")
    # From off-center-ellipse: the value's normalized interval is narrower than
    # the image of alpha's interval, so the producer must refine alpha until the
    # image fits inside it, as the checker requires.
    var octic_coefficients: List[QQ] = [QQ(Rational(99, 8960)), QQ(Rational(0)), QQ(Rational(-837, 4480)), QQ(Rational(0)), QQ(Rational(10179, 8960)), QQ(Rational(0)), QQ(Rational(-2207, 840)), QQ(Rational(0)), QQ(Rational(1))]
    var octic = DensePoly[QQ].from_coefficients(0, NoContext(), octic_coefficients^, work)
    var ellipse: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), RealAlgebraic(octic^, QQ(Rational(-1, 2)), QQ(Rational(-1, 4)), False))]
    var quadric = parse_polynomial("1/4 - 3/4*x0^2", order, work)
    var at_ellipse = algebraic_evaluate_sparse(quadric, ellipse, work)
    c.check(check_algebraic_evaluate_sparse(quadric, ellipse, at_ellipse, work), "direct evaluation whose value interval was narrowed by normalization")
    c.finish()
