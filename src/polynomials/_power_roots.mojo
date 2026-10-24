"""Exact roots of integer polynomials and the multiplicity hint that finds
perfect powers."""
from apn_mojo import Integer, Rational
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from .coefficients import ZZ, QQ, NoContext, ModInt, Modulus
from .primes import _residue
from .sparse import SparsePoly, Term
from .dense import DensePoly
from .division import divide_exact
from .certificates import _equal
from .squarefree import squarefree
from ._views import _join
from ._grouping import _without
from ._coefficient_ops import _gcd_coefficient, _exact


def _imul(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    _charge_size(a.magnitude_bit_length() // 64 + b.magnitude_bit_length() // 64 + 1, budget)
    try:
        return a * b
    except err:
        raise PolyError(PolyErrorKind.CAPACITY, String(err))


def _iadd(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    _charge_size(max(a.magnitude_bit_length(), b.magnitude_bit_length()) // 64 + 1, budget)
    try:
        return a + b
    except err:
        raise PolyError(PolyErrorKind.CAPACITY, String(err))


def _ipow(a: Integer, e: Int, mut budget: Budget) raises PolyError -> Integer:
    var result = Integer(1)
    for _ in range(e):
        budget.charge(1)
        result = _imul(result, a, budget)
    return result^


def _integer_root(a: Integer, e: Int, mut budget: Budget) raises PolyError -> Optional[Integer]:
    """The integer r with r^e = a, if one exists (r > 0 for even e)."""
    var negative = a.sign() < 0
    if negative and e % 2 == 0:
        return None
    var magnitude = a.copy()
    if negative:
        magnitude = -magnitude
    var lo = Integer(0)
    var hi: Integer
    try:
        hi = Integer(1) << (magnitude.magnitude_bit_length() // e + 1)
    except err:
        raise PolyError(PolyErrorKind.CAPACITY, String(err))
    while lo < hi:
        budget.charge(1)
        var mid: Integer
        try:
            mid = (lo + hi + Integer(1)) >> 1
        except err:
            raise PolyError(PolyErrorKind.CAPACITY, String(err))
        if _ipow(mid, e, budget) > magnitude:
            try:
                hi = mid - Integer(1)
            except err:
                raise PolyError(PolyErrorKind.CAPACITY, String(err))
        else:
            lo = mid^
    if not (_ipow(lo, e, budget) == magnitude):
        return None
    if negative:
        return -lo
    return lo^


def _gcd_int(a: Int, b: Int) -> Int:
    var x = a
    var y = b
    while y:
        var t = x % y
        x = y
        y = t
    return x


comptime _HINT_PRIME = UInt64(2147483647)
"""The prime 2^31 - 1, whose residue products fit a word."""


def _word_residue(a: Integer, p: UInt64, mut budget: Budget) raises PolyError -> UInt64:
    """a mod p in 0 ..< p, without a division of integers when a fits a word."""
    if a.magnitude_bit_length() > 62:
        return _residue(a, p, budget)
    var x: Int
    try:
        x = Int(a)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    if x >= 0:
        return UInt64(x) % p
    return (p - UInt64(-x) % p) % p


def _power_hint(f: SparsePoly[ZZ], main: Int64, mut budget: Budget) raises PolyError -> Int:
    """The gcd of the multiplicities of f at fixed small points for every
    variable but main, or 1 when that image loses degree or f's degrees
    rule out every power.

    The multiplicities come from the image modulo a prime, which keeps
    those over QQ unless the prime is unlucky; the caller's exact root check
    makes a wrong hint cost only time. An image losing degree modulo the
    prime takes the exact rational path."""
    # Products of residues below 2^31 fit a word.
    var p = _HINT_PRIME
    var image = List[UInt64]()
    var variables = List[Int64]()
    var degrees = List[Int]()
    # The exact path's points, v % 61 + 2, so this image is its reduction.
    var powers = List[List[UInt64]]()
    for _ in range(61):
        powers.append([UInt64(1)])
    for term in f.terms:
        budget.charge(1)
        var value = _word_residue(term.coefficient.value, p, budget)
        var k = 0
        for pair in term.monomial.pairs:
            budget.charge(1)
            var i = 0
            while i < len(variables) and variables[i] != pair[0]:
                i += 1
            if i == len(variables):
                variables.append(pair[0])
                degrees.append(0)
            degrees[i] = max(degrees[i], Int(pair[1]))
            if pair[0] == main:
                k = Int(pair[1])
            else:
                var v = Int(pair[0] % 61)
                while len(powers[v]) <= Int(pair[1]):
                    powers[v].append(powers[v][len(powers[v]) - 1] * UInt64(v + 2) % p)
                value = value * powers[v][Int(pair[1])] % p
        if k >= len(image):
            _allocation[UInt64](_sum(k, 1), budget)
            image.resize(k + 1, 0)
        image[k] = (image[k] + value) % p
    var degree = len(image) - 1
    # f = c h^e needs e to divide its degree in every variable.
    var common = 0
    for d in degrees:
        common = _gcd_int(common, d)
    if common < 2:
        return 1
    if image[degree] == 0:
        return _exact_hint(f, main, budget)
    var ring = Modulus(_prime_value=p)
    _allocation[ModInt](len(image), budget)
    var coefficients = List[ModInt](capacity=len(image))
    for w in image:
        coefficients.append(ModInt(w, ring))
    var blocks = squarefree(DensePoly[ModInt].from_coefficients(main, ring, coefficients^, budget), budget)
    var e = 0
    for entry in blocks.factors:
        budget.charge(1)
        e = _gcd_int(e, entry[1])
    return max(e, 1)


def _exact_hint(f: SparsePoly[ZZ], main: Int64, mut budget: Budget) raises PolyError -> Int:
    """_power_hint from the exact integer image over QQ."""
    var degree = f.degree_in(main, budget)
    _allocation[Integer](_sum(degree, 1), budget)
    var values = List[Integer](length=degree + 1, fill=Integer(0))
    # Distinct small points per variable index: 2, 3, 4, ...
    var powers = List[List[Integer]]()
    for _ in range(61):
        var row = List[Integer]()
        row.append(Integer(1))
        powers.append(row^)
    for term in f.terms:
        budget.charge(1)
        var value = term.coefficient.value.copy()
        var k = 0
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] == main:
                k = Int(pair[1])
            else:
                var v = Int(pair[0] % 61)
                while len(powers[v]) <= Int(pair[1]):
                    powers[v].append(_imul(powers[v][len(powers[v]) - 1], Integer(v + 2), budget))
                value = _imul(value, powers[v][Int(pair[1])], budget)
        values[k] = _iadd(values[k], value, budget)
    if values[degree].sign() == 0:
        return 1
    var coefficients = List[QQ]()
    for v in values:
        coefficients.append(QQ(Rational(v.copy())))
    var image = DensePoly[QQ].from_coefficients(main, NoContext(), coefficients^, budget)
    var blocks = squarefree(image, budget)
    var e = 0
    for entry in blocks.factors:
        budget.charge(1)
        e = _gcd_int(e, entry[1])
    return max(e, 1)


def _primitive_ordered(f: SparsePoly[ZZ], mut budget: Budget) raises PolyError -> SparsePoly[ZZ]:
    """_primitive_zz without re-sorting: the content stops at one, and
    dividing every coefficient by the same scalar keeps the term order."""
    var content = ZZ.zero(NoContext())
    for term in f.terms:
        budget.charge(1)
        content = _gcd_coefficient(content, term.coefficient, budget)
        if content.value == Integer(1):
            break
    if f.terms[len(f.terms) - 1].coefficient.value.sign() < 0:
        content = content.neg(budget)
    if content.value == Integer(1):
        return f._copy(budget)
    var result = SparsePoly[ZZ](f.order.copy(), NoContext())
    _allocation[Term[ZZ]](len(f.terms), budget)
    result.terms = List[Term[ZZ]](capacity=len(f.terms))
    for term in f.terms:
        budget.charge(1)
        _charge_size(len(term.monomial.pairs), budget)
        result.terms.append(Term[ZZ](term.monomial.copy(), _exact(term.coefficient, content, budget)))
    return result^


def _top_coefficients(f: SparsePoly[ZZ], main: Int64, low: Int, mut budget: Budget) raises PolyError -> List[SparsePoly[ZZ]]:
    """The coefficients of main^low .. main^deg in f, in that order."""
    var top = f.degree_in(main, budget)
    var buckets = List[List[Term[ZZ]]]()
    for _ in range(low, top + 1):
        buckets.append(List[Term[ZZ]]())
    for term in f.terms:
        budget.charge(1)
        var k = Int(term.monomial.exponent(main, budget))
        if k < low:
            continue
        _allocation[Term[ZZ]](1, budget)
        buckets[k - low].append(Term[ZZ](_without(term.monomial, main, budget), term.coefficient.copy()))
    var result = List[SparsePoly[ZZ]]()
    for i in range(len(buckets)):
        budget.charge(1)
        var terms = List[Term[ZZ]]()
        swap(terms, buckets[i])
        result.append(SparsePoly[ZZ].from_terms(terms^, f.order, NoContext(), budget))
    return result^


def _power_root(f: SparsePoly[ZZ], e: Int, mut budget: Budget) raises PolyError -> Optional[SparsePoly[ZZ]]:
    """h with h^e = f, if one exists.

    In the main variable, with reversed coefficient sequences g of f and a
    of h, h = f^(1/e) gives h' f = (1/e) h f' and so Miller's recurrence
    k e g_0 a_k = sum_(j=1..k) (j - e (k - j)) g_j a_(k-j); every division is
    exact when the root exists, and the candidate is verified by its power."""
    var variables = f.variables(budget)
    if len(variables) == 0:
        if f.is_zero():
            return None
        var r = _integer_root(f.terms[0].coefficient.value, e, budget)
        if not r:
            return None
        return SparsePoly[ZZ].constant(ZZ(r.take()), f.order, budget)
    var main = variables[len(variables) - 1]
    var n = f.degree_in(main, budget)
    if n % e != 0:
        return None
    var d = n // e
    # g[i] is the coefficient of main^(n - d + i); only the top d + 1 enter.
    var top = _top_coefficients(f, main, n - d, budget)
    var lead = _power_root(top[d], e, budget)
    if not lead:
        return None
    var a = List[SparsePoly[ZZ]]()
    a.append(lead.take())
    for k in range(1, d + 1):
        budget.charge(1)
        var total = SparsePoly[ZZ](f.order.copy(), NoContext())
        for j in range(1, min(k, n) + 1):
            budget.charge(1)
            var weight = j - e * (k - j)
            if weight == 0 or top[d - j].is_zero():
                continue
            total = total.add(top[d - j].mul(a[k - j], budget).scale(ZZ(Integer(weight)), budget), budget)
        var quotient = divide_exact(total, top[d].scale(ZZ(Integer(k * e)), budget), budget)
        if not quotient:
            return None
        a.append(quotient.take())
    var coefficients = List[SparsePoly[ZZ]]()
    for i in range(d + 1):
        coefficients.append(a[d - i]._copy(budget))
    var h = _join(coefficients, main, f, budget)
    if not _equal(h.pow(UInt64(e), budget), f, budget):
        return None
    return h^
