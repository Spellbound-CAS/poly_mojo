"""Public algebraic arithmetic and independent corruption checks."""
from std.collections import Optional
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, RealAlgebraic, AlgebraicArithmeticCertificate, RootInterval
from polynomials import real_root, real_rational, algebraic_compare, check_real_algebraic, check_algebraic_comparison
from polynomials import algebraic_add, algebraic_subtract, algebraic_multiply, algebraic_negate, algebraic_inverse, check_algebraic_arithmetic
from support.check import Checker
from polynomials._algebraic import _NormalizationCache, _root_cached


def main() raises:
    var work = Budget()
    var c = Checker("algebraic certificates")
    var zero = QQ.zero(NoContext())
    var one = QQ.one(NoContext())
    var two = QQ(Rational(2))
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), zero.copy(), one.copy()], work)
    var g = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-3)), zero.copy(), one.copy()], work)
    var a = real_root(f, one, two, work)
    var b = real_root(g, one, two, work)
    var cache = _NormalizationCache()
    var combined = f.mul(g, work)
    var first_work = Budget()
    var normalized = _root_cached(combined, QQ(Rational(7, 5)), QQ(Rational(3, 2)), cache, first_work)
    var hit_work = Budget()
    var repeated = _root_cached(combined, QQ(Rational(7, 5)), QQ(Rational(3, 2)), cache, hit_work)
    c.check(check_real_algebraic(repeated, work) and normalized.polynomial == repeated.polynomial and len(cache.polynomials) == 1, "normalization cache reuses factorization")
    c.check(hit_work.remaining() > first_work.remaining(), "normalization cache saves measured work")
    normalized.polynomial.coefficients[0] = zero.copy()
    var intact = _root_cached(combined, QQ(Rational(7, 5)), QQ(Rational(3, 2)), cache, work)
    c.check(check_real_algebraic(intact, work), "returned root witnesses do not alias cached factorizations")
    var conjugate = _root_cached(combined, QQ(Rational(-3, 2)), QQ(Rational(-7, 5)), cache, work)
    c.check(algebraic_compare(conjugate, a, work).relation == -1 and len(cache.polynomials) == 1, "shared polynomial cache keeps conjugate intervals distinct")
    var invalid = False
    try:
        _ = _root_cached(combined, one, two, cache, work)
    except e:
        invalid = e.kind == PolyErrorKind.DOMAIN
    c.check(invalid, "cache hit still rejects a nonisolating interval")
    var comparison = algebraic_compare(a, b, work)
    c.check(comparison.relation == -1 and check_algebraic_comparison(a, b, comparison, work), "public comparison")
    for operation in range(5):
        var result: AlgebraicArithmeticCertificate
        if operation == 0:
            result = algebraic_add(a, b, work)
        elif operation == 1:
            result = algebraic_subtract(a, b, work)
        elif operation == 2:
            result = algebraic_multiply(a, b, work)
        elif operation == 3:
            result = algebraic_negate(a, work)
        else:
            result = algebraic_inverse(a, work)
        var second = Optional[RealAlgebraic]()
        if operation < 3:
            second = b.copy()
        c.check(check_real_algebraic(result.value, work), "arithmetic result is a valid root")
        c.check(check_algebraic_arithmetic(a, second, result, work), "arithmetic certificate accepted")
        var broken = result.copy()
        broken.operation = (operation + 1) % 5
        c.check(not check_algebraic_arithmetic(a, second, broken, work), "changed operation rejected")
        broken = result.copy()
        broken.value_factor.coefficients[0] = broken.value_factor.coefficients[0].add(one, work)
        c.check(not check_algebraic_arithmetic(a, second, broken, work), "changed annihilator quotient rejected")
        if result.selection:
            broken = result.copy()
            broken.selection = RootInterval(QQ(Rational(100)), QQ(Rational(101)))
            c.check(not check_algebraic_arithmetic(a, second, broken, work), "selection window away from the result rejected")
            broken = result.copy()
            broken.selection.value().high = broken.selection.value().low.copy()
            c.check(not check_algebraic_arithmetic(a, second, broken, work), "empty selection window rejected")
        if operation < 3:
            broken = result.copy()
            var elimination = broken.elimination.take()
            elimination.cofactors[0] = elimination.cofactors[0].neg(work)
            broken.elimination = elimination^
            c.check(not check_algebraic_arithmetic(a, second, broken, work), "corrupted elimination identity rejected")
        var tiny = Budget(0)
        try:
            _ = check_algebraic_arithmetic(a, second, result, tiny)
            c.check(False, "checker budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "checker budget enforced")
    var z = real_rational(zero, work)
    var annihilated = algebraic_multiply(a, z, work)
    c.check(annihilated.value.rational and annihilated.value.low.is_zero() and check_algebraic_arithmetic(a, z.copy(), annihilated, work), "zero multiplication exact certificate")
    var twice = algebraic_add(z, real_rational(two, work), work)
    c.check(check_algebraic_arithmetic(z, real_rational(two, work), twice, work), "rational arithmetic certificate")
    var identity = algebraic_multiply(a, a, work)
    c.check(identity.value.rational and identity.value.low == two and check_algebraic_arithmetic(a, a.copy(), identity, work), "algebraic square collapses with a valid proof")
    try:
        _ = algebraic_inverse(z, work)
        c.check(False, "inverse of zero accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "inverse of zero fails closed")
    c.finish()
