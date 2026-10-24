"""Factorization of the binomials x^n - 1 and x^n + 1 into cyclotomic polynomials.

x^n - 1 is the product of Phi_d over the divisors d of n, and x^n + 1 the
product of Phi_d over the divisors d of 2n that do not divide n; every
Phi_d is irreducible over QQ. For d > 1, Phi_d = prod_(e | d) (1 - x^e)^mu(d/e)
as a power series, and its degree is phi(d), so each factor is one pass
per divisor over phi(d) + 1 coefficients: multiplying by 1 - x^e or by the
geometric series 1 / (1 - x^e). Coefficients stay in words; a step that
would leave them aborts, and the general algorithm takes over.
"""
from apn_mojo import Integer
from .coefficients import ZZ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, _allocation, _charge_size, _product

comptime _WORD_LIMIT = Int64(1) << 60
"""Series coefficients beyond this magnitude abort the shortcut."""


def _divisors(n: Int) -> List[Int]:
    var small = List[Int]()
    var large = List[Int]()
    var d = 1
    while d * d <= n:
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
        d += 1
    for i in range(len(large) - 1, -1, -1):
        small.append(large[i])
    return small^


def _mobius(n: Int) -> Int:
    var m = n
    var result = 1
    var p = 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            result = -result
        p += 1
    if m > 1:
        result = -result
    return result


def _totient(n: Int) -> Int:
    var m = n
    var result = n
    var p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def _cyclotomic(d: Int, variable: Int64, mut budget: Budget) raises PolyError -> Optional[DensePoly[ZZ]]:
    """Phi_d, or None when a series coefficient would leave a word."""
    if d == 1:
        return DensePoly[ZZ].from_coefficients(variable, NoContext(), [ZZ(Integer(-1)), ZZ(Integer(1))], budget)
    var degree = _totient(d)
    var divisors = _divisors(d)
    _charge_size(_product(len(divisors), degree + 1), budget)
    _allocation[Int64](degree + 1, budget)
    var c = List[Int64](length=degree + 1, fill=0)
    c[0] = 1
    for e in divisors:
        var mu = _mobius(d // e)
        if mu == 1:
            # c <- c (1 - x^e), from the top down.
            var k = degree
            while k >= e:
                var v = c[k] - c[k - e]
                if v >= _WORD_LIMIT or v <= -_WORD_LIMIT:
                    return None
                c[k] = v
                k -= 1
        elif mu == -1:
            # c <- c / (1 - x^e), from the bottom up.
            for k in range(e, degree + 1):
                var v = c[k] + c[k - e]
                if v >= _WORD_LIMIT or v <= -_WORD_LIMIT:
                    return None
                c[k] = v
    _allocation[ZZ](degree + 1, budget)
    var coefficients = List[ZZ](capacity=degree + 1)
    for k in range(degree + 1):
        coefficients.append(ZZ(Integer(Int(c[k]))))
    return DensePoly[ZZ].from_coefficients(variable, NoContext(), coefficients^, budget)


def _binomial_factors(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Optional[List[DensePoly[ZZ]]]:
    """The cyclotomic factors of f = x^n - 1 or x^n + 1 with n >= 2, or None for any other f."""
    var n = f.degree()
    if n < 2 or not (f.coefficients[n].value == Integer(1)):
        return None
    var minus = f.coefficients[0].value == Integer(-1)
    if not minus and not (f.coefficients[0].value == Integer(1)):
        return None
    for i in range(1, n):
        budget.charge(1)
        if not f.coefficients[i].is_zero():
            return None
    var factors = List[DensePoly[ZZ]]()
    if minus:
        for d in _divisors(n):
            var phi = _cyclotomic(d, f.variable, budget)
            if not phi:
                return None
            factors.append(phi.take())
    else:
        for d in _divisors(2 * n):
            if n % d == 0:
                continue
            var phi = _cyclotomic(d, f.variable, budget)
            if not phi:
                return None
            factors.append(phi.take())
    return factors^
