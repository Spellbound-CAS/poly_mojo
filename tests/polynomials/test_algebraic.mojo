"""Hand-computed algebraic identities and fail-closed boundary cases."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, sturm_sequence
from polynomials._algebraic import _root, _rational, _compare, _sign, _roots, _copy_real, _comparison
from support.algebraic import _add, _mul, _neg, _inverse, _evaluate, _roots_algebraic, _roots_at
from polynomials.algebraic_types import RealAlgebraic as _Real, _check_real, _check_comparison
from support.check import Checker
from support.corpus import parse_polynomial
from polynomials import MonomialOrder


def q(n: Int) -> QQ:
    return QQ(Rational(n))


def poly(var coefficients: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var values = List[QQ]()
    for n in coefficients:
        values.append(q(n))
    return DensePoly[QQ].from_coefficients(0, NoContext(), values^, work)


def main() raises:
    var work = Budget()
    var c = Checker("algebraic")
    var f = poly([-2, 0, 1], work)
    var a = _root(f, q(1), q(2), work)
    var b = _root(poly([-3, 0, 1], work), q(1), q(2), work)
    var zero = _rational(q(0), work)
    var one = _rational(q(1), work)
    var two = _rational(q(2), work)
    c.check(_compare(_mul(a, a, work), two, work) == 0, "sqrt(2) squared is 2")
    c.check(_compare(_add(a, _neg(a, work), work), zero, work) == 0, "additive inverse")
    c.check(_compare(_mul(a, _inverse(a, work), work), one, work) == 0, "multiplicative inverse")
    var sum = _add(a, b, work)
    c.check(_sign(_evaluate(poly([1, 0, -10, 0, 1], work), sum, work), work) == 0, "sqrt(2)+sqrt(3) annihilated by x^4-10*x^2+1")
    c.check(_compare(a, b, work) < 0 and _compare(b, a, work) > 0, "overlapping intervals compare exactly")
    var same = _root(f.pow(3, work), q(0), q(3), work)
    c.check(_compare(a, same, work) == 0, "repeated defining polynomial gives same value")
    c.check(_check_real(a, work) and _check_real(same, work) and _check_real(zero, work), "independent value witnesses")
    var equality = _comparison(a, same, work)
    c.check(_check_comparison(a, same, equality, work), "independent equality witness")
    equality.left_cofactor.coefficients[0] = q(2)
    c.check(not _check_comparison(a, same, equality, work), "corrupted equality cofactor rejected")
    var ordering = _comparison(a, b, work)
    c.check(_check_comparison(a, b, ordering, work), "independent inequality witness")
    ordering.relation = 1
    c.check(not _check_comparison(a, b, ordering, work), "reversed inequality rejected")
    var malformed = _copy_real(a, work)
    malformed.polynomial.coefficients[0] = q(-3)
    c.check(not _check_real(malformed, work), "corrupted defining polynomial rejected")
    var nonminimal_poly = f.mul(poly([-3, 1], work), work)
    var nonminimal = _Real(nonminimal_poly.copy(), q(1), q(2), sturm_sequence(nonminimal_poly, work), False)
    c.check(_check_real(nonminimal, work), "valid nonminimal defining polynomial")
    c.check(_check_comparison(a, nonminimal, _comparison(a, nonminimal, work), work), "shared-root equality with distinct definitions")
    var zero_poly = poly([0, 0, 1], work)
    var nonminimal_zero = _Real(zero_poly.copy(), q(-1), q(1), sturm_sequence(zero_poly, work), False)
    c.check(_check_comparison(zero, nonminimal_zero, _comparison(zero, nonminimal_zero, work), work), "rational equality with an interval witness")
    c.check(_check_comparison(nonminimal_zero, a, _comparison(nonminimal_zero, a, work), work), "nonminimal rational refinement retains definition")
    var negative = _root(f, q(-2), q(-1), work)
    c.check(_compare(_neg(a, work), negative, work) == 0, "negative conjugate")
    c.check(_sign(a, work) == 1 and _sign(negative, work) == -1 and _sign(zero, work) == 0, "all signs")
    var rational = _root(poly([-1, 3], work), q(0), q(1), work)
    c.check(rational.rational and rational.low.value == Rational(1, 3), "non-dyadic rational collapse")
    var original = _copy_real(a, work)
    for limit in [UInt64(0), UInt64(100)]:
        var tiny = Budget(limit)
        try:
            _ = _add(a, b, tiny)
            c.check(False, "budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "budget fails closed")
    c.check(a.low == original.low and a.high == original.high, "borrowed intervals unchanged")
    for invalid in range(5):
        var low = q(1)
        var high = q(2)
        var constant = poly([1], work)
        try:
            if invalid == 0:
                _ = _inverse(zero, work)
            elif invalid == 1:
                _ = _root(f, high, low, work)
            elif invalid == 2:
                _ = _root(f, q(-2), high, work)
            elif invalid == 3:
                _ = _root(constant, low, high, work)
            else:
                _ = _root(f, low, low, work)
            c.check(False, "invalid algebraic operation accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "invalid algebraic operation rejected")
    var roots = _roots(poly([4, 0, -4, 0, 1], work), work)
    c.check(len(roots) == 2, "distinct roots of repeated polynomial")
    c.check(_compare(roots[0], negative, work) == 0 and _compare(roots[1], a, work) == 0, "roots sorted by exact value")
    var fourth = _roots_algebraic([_neg(a, work), _copy_real(zero, work), _copy_real(one, work)], work)
    c.check(len(fourth) == 2, "x^2-sqrt(2) has two roots")
    for value in fourth:
        c.check(_compare(_mul(value, value, work), a, work) == 0, "fourth roots square to sqrt(2)")
    var no_roots = _roots_algebraic([_copy_real(a, work), _copy_real(zero, work), _copy_real(one, work)], work)
    c.check(len(no_roots) == 0, "x^2+sqrt(2) has no real roots")
    var assigned: List[Tuple[Int64, _Real]] = [(Int64(0), _copy_real(a, work))]
    var at = _roots_at(parse_polynomial("x1^2 - x0", MonomialOrder(), work), 1, assigned, work)
    c.check(len(at) == 2 and _compare(at[1], fourth[1], work) == 0, "roots at an algebraic point")
    var drop = _roots_at(parse_polynomial("x0^2*x1^2 - 2*x1^2 + x1 - 1", MonomialOrder(), work), 1, assigned, work)
    c.check(len(drop) == 1 and _compare(drop[0], one, work) == 0, "leading coefficient vanishes after substitution")
    var vanishing = parse_polynomial("x0^2*x1 - 2*x1", MonomialOrder(), work)
    try:
        _ = _roots_at(vanishing, 1, assigned, work)
        c.check(False, "identically zero specialization accepted as a finite root list")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "identically zero specialization fails closed")
    c.finish()
