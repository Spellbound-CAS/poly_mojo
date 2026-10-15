"""Lazard nullification and minimal valuation witnesses."""
from apn_mojo import Rational
from polynomials import Budget, MonomialOrder, QQ, RealAlgebraic, real_rational, lazard_evaluate, check_lazard_evaluation
from polynomials.algebraic_specialization import _coefficients
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var work = Budget()
    var c = Checker("Lazard evaluation")
    var order = MonomialOrder()
    var zero = real_rational(QQ(Rational(0)), work)
    var point: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), zero.copy()), (Int64(1), zero.copy())]
    var f = parse_polynomial("x0^2*x1*x2 - x0^2*x1 + x0^3", order, work)
    var proof = lazard_evaluate(f, 2, point, work)
    c.check(check_lazard_evaluation(f, 2, point, proof, work), "successive exact Taylor witnesses")
    c.check(len(proof.orders) == 2 and proof.orders[0] == 2 and proof.orders[1] == 1, "valuation is (2,1)")
    var coefficients = _coefficients(proof.evaluations[len(proof.evaluations) - 1], 2, work)
    c.check(len(coefficients) == 2 and coefficients[0].low.value == Rational(-2) and coefficients[1].low.value == Rational(2), "Lazard lift is positive factorial times (z-1)")
    var bad = proof.copy()
    bad.orders[0] = 1
    c.check(not check_lazard_evaluation(f, 2, point, bad, work), "nonminimal or premature valuation rejected")
    c.finish()
