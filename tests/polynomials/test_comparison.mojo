"""Exact comparison cases, independent certificates, and failure atomicity."""
from apn_mojo import Rational
from polynomials import Budget, QQ, DensePoly, MonomialOrder, RealAlgebraic, PolyErrorKind, sturm_sequence, real_rational, algebraic_compare, check_algebraic_comparison
from polynomials._algebraic import _compare
from support.corpus import parse_polynomial
from support.check import Checker


def interval(expression: String, low: Int, high: Int, mut budget: Budget) raises -> RealAlgebraic:
    var f = DensePoly[QQ].from_sparse(parse_polynomial(expression, MonomialOrder(), budget), budget)
    return RealAlgebraic(f.copy(), QQ(Rational(low)), QQ(Rational(high)), sturm_sequence(f, budget), False)


def main() raises:
    var c = Checker("comparison")
    var work = Budget()
    var a = interval("x0^2 - 2", 1, 2, work)
    var others: List[RealAlgebraic] = [interval("x0^2 - 8", 2, 3, work), interval("x0^2 - 3", 1, 2, work), interval("x0^3 - 3*x0^2 - 2*x0 + 6", 1, 2, work), interval("x0^2 - 2", -2, 1, work), real_rational(QQ(Rational(1)), work)]
    var expected: List[Int] = [-1, -1, 0, 1, 1]
    for i in range(len(others)):
        var proof = algebraic_compare(a, others[i], work)
        c.check(proof.relation == expected[i] and _compare(a, others[i], work) == expected[i], "known exact ordering")
        c.check(check_algebraic_comparison(a, others[i], proof, work), "independent comparison check")
        c.check(algebraic_compare(others[i], a, work).relation == -expected[i], "antisymmetry")
        proof.relation = 1 if expected[i] != 1 else -1
        c.check(not check_algebraic_comparison(a, others[i], proof, work), "forged ordering rejected")
    var negative = interval("x0^2 - 2", -2, 1, work)
    var positive = interval("x0^2 - 2", -1, 2, work)
    var distinct = algebraic_compare(negative, positive, work)
    c.check(distinct.relation == -1 and check_algebraic_comparison(negative, positive, distinct, work), "overlapping intervals of conjugate roots stay distinct")
    var repeated = interval("x0^4 - 4*x0^2 + 4", 1, 2, work)
    var equality = algebraic_compare(a, repeated, work)
    c.check(equality.relation == 0 and check_algebraic_comparison(a, repeated, equality, work), "repeated nonminimal definition")
    equality.shared.value().count = 2
    c.check(not check_algebraic_comparison(a, repeated, equality, work), "corrupted shared-root count rejected")
    var point = real_rational(QQ(Rational(0)), work)
    var nonminimal_zero = interval("x0^2", -1, 1, work)
    var rational_equality = algebraic_compare(point, nonminimal_zero, work)
    c.check(rational_equality.relation == 0 and check_algebraic_comparison(point, nonminimal_zero, rational_equality, work), "point and interval equality")
    c.check(algebraic_compare(point, point, work).relation == 0, "rational equality")
    var bad = a.copy()
    bad.sequence.polynomials[0].coefficients[0] = QQ(Rational(-9))
    try:
        _ = algebraic_compare(bad, others[0], work)
        c.check(False, "separated invalid operand accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "fast comparison validates public operands")
    var saved_low = positive.low.copy()
    var saved_high = positive.high.copy()
    for limit in [UInt64(0), UInt64(1), UInt64(100), UInt64(1000), UInt64(10000)]:
        var tiny = Budget(limit)
        try:
            var proof = algebraic_compare(negative, positive, tiny)
            c.check(check_algebraic_comparison(negative, positive, proof, work), "completed limited comparison has full evidence")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "comparison exhaustion remains LIMIT")
        c.check(positive.low == saved_low and positive.high == saved_high, "failed comparison preserves borrowed intervals")
    c.finish()
