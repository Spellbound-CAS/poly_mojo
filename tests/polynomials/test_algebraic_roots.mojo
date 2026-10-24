"""Algebraic-coefficient root completeness and corrupted witness rejection."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, RealAlgebraic, real_rational, real_root
from polynomials import algebraic_roots, check_algebraic_roots, algebraic_evaluate, check_algebraic_evaluation, algebraic_compare
from support.check import Checker


def main() raises:
    var work = Budget()
    var c = Checker("algebraic roots")
    var qzero = QQ.zero(NoContext())
    var qone = QQ.one(NoContext())
    var zero = real_rational(qzero, work)
    var one = real_rational(qone, work)
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), qzero.copy(), qone.copy()], work)
    var a = real_root(f, qone, QQ(Rational(2)), work)
    var minus = real_root(f, QQ(Rational(-2)), QQ(Rational(-1)), work)
    var coefficients: List[RealAlgebraic] = [minus.copy(), zero.copy(), one.copy()]
    var proof = algebraic_roots(coefficients, work)
    c.check(len(proof.roots) == 2 and len(proof.candidates) == 2, "x^2-sqrt(2) has two real roots")
    c.check(check_algebraic_roots(coefficients, proof, work), "complete fourth-root certificate")
    var evaluation = algebraic_evaluate(coefficients, proof.roots[1], work)
    c.check(evaluation.value.rational and evaluation.value.low.is_zero() and check_algebraic_evaluation(coefficients, proof.roots[1], evaluation, work), "fourth root evaluates exactly to zero")
    var square_coefficients: List[RealAlgebraic] = [zero.copy(), zero.copy(), one.copy()]
    var square = algebraic_evaluate(square_coefficients, proof.roots[1], work)
    c.check(algebraic_compare(square.value, a, work).relation == 0, "positive fourth root squares to sqrt(2)")
    for corruption in range(5):
        var broken = proof.copy()
        if corruption == 0:
            _ = broken.all_roots.certificate.leaves.pop()
        elif corruption == 1:
            broken.accepted[0] = 1
        elif corruption == 2:
            broken.candidate_factors[0].coefficients[0] = broken.candidate_factors[0].coefficients[0].add(qone, work)
        elif corruption == 3:
            broken.evaluations[0].products[0].operation = 0
        else:
            _ = broken.candidates.pop()
        c.check(not check_algebraic_roots(coefficients, broken, work), "corrupted complete-root witness rejected")
    var positive: List[RealAlgebraic] = [a.copy(), zero.copy(), one.copy()]
    var empty = algebraic_roots(positive, work)
    c.check(len(empty.roots) == 0 and len(empty.candidates) == 2 and check_algebraic_roots(positive, empty, work), "conjugate candidates rejected for x^2+sqrt(2)")
    var constant: List[RealAlgebraic] = [a.copy(), zero.copy()]
    var constant_result = algebraic_roots(constant, work)
    c.check(len(constant_result.roots) == 0 and check_algebraic_roots(constant, constant_result, work), "nonzero algebraic constant has no roots")
    c.finish()
