"""Resultants and the complete indexed subresultant family."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum, _product
from .coefficients import CoefficientRing, ZZ, ModInt
from .sparse import SparsePoly
from .dense import DensePoly
from ._views import _variables, _constant
from apn_mojo import Integer
from .primes import prime
from ._nmod_euclid import _Field, _SmallField, _SMALL_PRIME, _resultant, _resultant_words, _resultant_in, _trim_words
from ._univariate_modular import _Words, _hadamard_bits
from ._multivariate_resultant import _integer_resultant, _field_resultant
from .certificates import ResultantResult
from .division import divide_exact
from ._tower import _zero, _lc, _exact_poly
from ._prs import _prs


def subresultants[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> List[SparsePoly[C]]:
    """Compute S_j for every index below the smaller degree, including zeros.

    Parameters:
        C: ZZ or a coefficient field.

    Args:
        f: The first polynomial.
        g: The second polynomial in the same ring and order.
        variable: The elimination variable.
        budget: The shared work allowance.

    Returns:
        The indexed subresultants S_0 through S_(min(deg(f),deg(g))-1).

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    var n = f.degree_in(variable, budget)
    var m = g.degree_in(variable, budget)
    var count = max(min(n, m), 0)
    _allocation[SparsePoly[C]](count, budget)
    var result = List[SparsePoly[C]](length=count, fill=_zero(f, budget))
    if not count:
        return result^
    var sequence = _prs(f, g, variable, False, budget)
    for entry in sequence:
        budget.charge(1)
        var k = entry.polynomial.degree_in(variable, budget)
        result[entry.index] = entry.polynomial._copy(budget)
        if k != entry.index:
            var leading = _lc(entry.polynomial, variable, budget)
            var multiplier = divide_exact(entry.principal, leading, budget)
            if multiplier:
                result[k] = entry.polynomial.mul(multiplier.value(), budget)
            else:
                # A defective PRS scale can be a rational function. The full
                # scaled subresultant is polynomial even when the scale is not.
                result[k] = _exact_poly(entry.polynomial.mul(entry.principal, budget), leading, budget)
    if n < m:
        for j in range(count):
            budget.charge(1)
            if ((n - j) % 2) and ((m - j) % 2):
                result[j] = result[j].neg(budget)
    return result^


def resultant[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> ResultantResult[C]:
    """Compute the resultant and multipliers proving its ideal membership.

    Parameters:
        C: ZZ or a coefficient field.

    Args:
        f: The first polynomial.
        g: The second polynomial in the same ring and order.
        variable: The elimination variable.
        budget: The shared work allowance.

    Returns:
        The Sylvester-convention resultant and its cofactor identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or two nonzero constants in the variable, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    var zero = _zero(f, budget)
    if f.is_zero() or g.is_zero():
        return ResultantResult[C](zero._copy(budget), (zero._copy(budget), zero^))
    var n = f.degree_in(variable, budget)
    var m = g.degree_in(variable, budget)
    if n == 0 and m == 0:
        # E-1 assumes positive degrees. For two constants the empty
        # determinant is 1, which need not lie in their ideal over ZZ.
        raise PolyError(PolyErrorKind.DOMAIN, "resultant identity requires a positive degree")
    if n == 0:
        var u = f.pow(UInt64(m - 1), budget)
        return ResultantResult[C](u.mul(f, budget), (u^, zero^))
    if m == 0:
        var v = g.pow(UInt64(n - 1), budget)
        return ResultantResult[C](v.mul(g, budget), (zero^, v^))
    var sequence = _prs(f, g, variable, True, budget, history=False)
    if not len(sequence) or sequence[len(sequence) - 1].polynomial.degree_in(variable, budget) != 0:
        return ResultantResult[C](zero._copy(budget), (zero._copy(budget), zero^))
    ref last = sequence[len(sequence) - 1]
    var multiplier = divide_exact(last.principal, last.polynomial, budget)
    var u = _zero(f, budget)
    var v = _zero(f, budget)
    if multiplier:
        u = last.u.mul(multiplier.value(), budget)
        v = last.v.mul(multiplier.value(), budget)
    else:
        # Divide only after multiplying: last.principal / last.polynomial
        # need not be polynomial, while both scaled Bezout cofactors are.
        u = _exact_poly(last.u.mul(last.principal, budget), last.polynomial, budget)
        v = _exact_poly(last.v.mul(last.principal, budget), last.polynomial, budget)
    var value = last.principal._copy(budget)
    if n < m and (n % 2) and (m % 2):
        value = value.neg(budget)
        u = u.neg(budget)
        v = v.neg(budget)
    return ResultantResult[C](value^, (u^, v^))


struct _Magnitudes(Movable):
    """Coefficients below 2^64 in magnitude, with their signs, for residues by one Montgomery entry each."""
    var values: List[UInt64]
    var negative: List[Bool]

    def __init__(out self):
        self.values = List[UInt64]()
        self.negative = List[Bool]()

    def __init__(out self, f: DensePoly[ZZ], mut budget: Budget) raises PolyError:
        _allocation[UInt64](len(f.coefficients), budget)
        self.values = List[UInt64](capacity=len(f.coefficients))
        self.negative = List[Bool](capacity=len(f.coefficients))
        for i in range(len(f.coefficients)):
            ref c = f.coefficients[i].value
            self.values.append(UInt64(c._word(0)) | (UInt64(c._word(1)) << 32))
            self.negative.append(c.sign() < 0)

    @staticmethod
    def fits(f: DensePoly[ZZ]) -> Bool:
        for i in range(len(f.coefficients)):
            if f.coefficients[i].value._word_count() > 2:
                return False
        return True

    def residues(self, field: _Field, mut out: List[UInt64]):
        """Every coefficient modulo field.p in Montgomery form, into out; value r2 < 2^64 p, so one reduction suffices."""
        out.clear()
        for i in range(len(self.values)):
            var r = field.enter(self.values[i])
            out.append(field.sub(0, r) if self.negative[i] else r)


def _inverses_mod(primes: List[UInt64], values: List[UInt64], mut budget: Budget) raises PolyError -> List[UInt64]:
    """Standard residues 1/values[i] modulo odd primes[i] below 2^63, for values[i] below primes[i].

    A branchless binary extended gcd on eight lanes at once keeps
    x1 a = u and x2 a = v modulo p. Each step makes u even, by subtracting
    the smaller of u and v from the larger when u is odd, then halves u,
    so the bit lengths of u and v drop by at least one together: after
    bits(a) + bits(p) steps u = 0 and v = gcd(a, p), which must be 1.
    """
    comptime W = 8
    var k = len(primes)
    _allocation[UInt64](k, budget)
    var out = List[UInt64](length=k, fill=0)
    var top = UInt64(0)
    for q in primes:
        top = max(top, q)
    var steps = 0
    while (top >> UInt64(steps)) != 0:
        steps += 1
    steps *= 2
    var zero = SIMD[DType.uint64, W](0)
    var one = SIMD[DType.uint64, W](1)
    var base = 0
    while base < k:
        _charge_size(steps, budget)
        var p = SIMD[DType.uint64, W](3)
        var u = SIMD[DType.uint64, W](1)
        for j in range(min(W, k - base)):
            p[j] = primes[base + j]
            u[j] = values[base + j]
        var v = p
        var x1 = one
        var x2 = zero
        for _ in range(steps):
            var odd = (u & one).ne(zero)
            var exchange = odd & u.lt(v)
            var nu = exchange.select(v, u)
            var nv = exchange.select(u, v)
            var nx1 = exchange.select(x2, x1)
            var nx2 = exchange.select(x1, x2)
            var d = nx1 - nx2
            u = odd.select(nu - nv, nu) >> 1
            x1 = odd.select(nx1.lt(nx2).select(d + p, d), nx1)
            # x1 / 2 mod p; x1 + p stays below 2^64.
            x1 = (x1 + (x1 & one) * p) >> 1
            v = nv
            x2 = nx2
        for j in range(min(W, k - base)):
            if v[j] != 1:
                raise PolyError(PolyErrorKind.DOMAIN, "residue is not invertible")
            out[base + j] = x2[j]
        base += W
    return out^


def _mixed_radix(primes: List[UInt64], digits: List[UInt64], negative: List[Bool], mut budget: Budget) raises PolyError -> Integer:
    """The sum of the signed digits d_i times the products of the earlier primes, by Horner's rule in 64-bit two's complement words.

    Every prime is below 2^63, so the value and each partial sum fit in one
    word per prime plus a sign word; arithmetic modulo 2^(64 L) is exact there.
    """
    var k = len(primes)
    var length = k + 1
    _charge_size(_product(k, length), budget)
    _allocation[UInt64](length, budget)
    var words = List[UInt64](length=length, fill=0)
    var i = k - 1
    while i >= 0:
        var p = UInt128(primes[i])
        var carry = UInt128(0)
        for j in range(length):
            var t = UInt128(words[j]) * p + carry
            words[j] = UInt64(t)
            carry = t >> 64
        var d = digits[i]
        if negative[i]:
            for j in range(length):
                var w = words[j]
                words[j] = w - d
                d = UInt64(1) if w < d else UInt64(0)
                if d == 0:
                    break
        else:
            for j in range(length):
                var w = words[j] + d
                words[j] = w
                d = UInt64(1) if w < d else UInt64(0)
                if d == 0:
                    break
        i -= 1
    var sign = (words[length - 1] >> 63) != 0
    if sign:
        var carry = UInt64(1)
        for j in range(length):
            var w = ~words[j] + carry
            carry = UInt64(1) if (w == 0 and carry == 1) else UInt64(0)
            words[j] = w
    _allocation[UInt32](2 * length, budget)
    var halves = List[UInt32](capacity=2 * length)
    for j in range(length):
        halves.append(UInt32(words[j] & 0xFFFFFFFF))
        halves.append(UInt32(words[j] >> 32))
    return Integer._from_words(halves^, sign)


def _modular_resultant(fd: DensePoly[ZZ], gd: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Integer:
    """The resultant of two integer polynomials of positive degree, by images modulo the cached 62-bit primes up to the Hadamard bound.

    Each image stays a fraction num / den from the word kernel. Garner's
    mixed-radix digit for p is (num - acc den) / (den P) mod p, with acc the
    value of the earlier digits and P the product of the earlier primes,
    both modulo p; the inversions do not depend on the digits, so all of
    them run together first. Balanced digits make the final Horner sum the
    symmetric representative.
    """
    var m = fd.degree()
    var n = gd.degree()
    # |res| <= ||f||^n ||g||^m; the symmetric range needs one more bit.
    var bits = _sum(_sum(_product(n, _hadamard_bits(fd, budget)), _product(m, _hadamard_bits(gd, budget))), 2)
    var narrow = _Magnitudes.fits(fd) and _Magnitudes.fits(gd)
    var fm = _Magnitudes()
    var gm = _Magnitudes()
    var fw = Optional[_Words]()
    var gw = Optional[_Words]()
    if narrow:
        fm = _Magnitudes(fd, budget)
        gm = _Magnitudes(gd, budget)
    else:
        fw = _Words(fd, budget)
        gw = _Words(gd, budget)
    var a = List[UInt64](capacity=m + 1)
    var b = List[UInt64](capacity=n + 1)
    # Each image covers at least 61 bits; a prime that divides a leading coefficient adds one more.
    var expected = bits // 61 + 2
    _allocation[UInt64](_product(expected, 4), budget)
    var fields = List[_Field](capacity=expected)
    var primes = List[UInt64](capacity=expected)
    var numerators = List[UInt64](capacity=expected)
    var denominators = List[UInt64](capacity=expected)
    var covered = 0
    var index = 0
    while covered < bits:
        budget.charge(1)
        var p = prime(index, budget)
        index = _sum(index, 1)
        var field = _Field.of(p)
        if narrow:
            _charge_size(m + n + 2, budget)
            fm.residues(field, a)
            gm.residues(field, b)
        else:
            a = fw.value().residues(field, budget)
            b = gw.value().residues(field, budget)
        if a[m] == 0 or b[n] == 0:
            continue
        var fraction = _resultant_in(a, b, field, budget)
        fields.append(field)
        primes.append(p)
        numerators.append(fraction[0])
        denominators.append(fraction[1])
        var floor_log = 0
        while (p >> UInt64(floor_log + 1)) != 0:
            floor_log += 1
        covered = _sum(covered, floor_log)
    var k = len(primes)
    _charge_size(_product(k, k), budget)
    # w_i = den_i P_i, standard, for the inversions.
    var products = List[UInt64](capacity=k)
    for i in range(k):
        ref field = fields[i]
        var power = field.one()
        for j in range(i):
            power = field.mul(power, field.enter(primes[j]))
        products.append(field.leave(field.mul(denominators[i], power)))
    var inverses = _inverses_mod(primes, products, budget)
    var digits = List[UInt64](capacity=k)
    var negative = List[Bool](capacity=k)
    for i in range(k):
        ref field = fields[i]
        # acc = d_0 + p_0 (d_1 + p_1 (... d_(i-1))) mod p_i.
        var acc = UInt64(0)
        var j = i - 1
        while j >= 0:
            var digit = field.enter(digits[j])
            acc = field.add(field.mul(acc, field.enter(primes[j])), field.sub(0, digit) if negative[j] else digit)
            j -= 1
        var t = field.leave(field.mul(field.sub(numerators[i], field.mul(acc, denominators[i])), field.enter(inverses[i])))
        var high = t > primes[i] >> 1
        digits.append(primes[i] - t if high else t)
        negative.append(high)
    return _mixed_radix(primes, digits, negative, budget)


def _word_resultant(f: DensePoly[ModInt], g: DensePoly[ModInt], mut budget: Budget) raises PolyError -> ModInt:
    _allocation[UInt64](len(f.coefficients) + len(g.coefficients), budget)
    var p = f.ring.value
    if p < _SMALL_PRIME:
        # Straight into the 32-bit Montgomery form, without the 64-bit
        # field's setup or the narrowing and widening around its kernel.
        var small = _SmallField.of(p)
        var a = List[UInt64](capacity=len(f.coefficients))
        for i in range(len(f.coefficients)):
            a.append(small.reduce(f.coefficients[i].value * small.r2))
        var b = List[UInt64](capacity=len(g.coefficients))
        for i in range(len(g.coefficients)):
            b.append(small.reduce(g.coefficients[i].value * small.r2))
        _trim_words(a)
        _trim_words(b)
        if not len(a) or not len(b):
            return ModInt.zero(f.ring)
        return ModInt(small.reduce(_resultant(a, b, small, budget)), f.ring)
    var field = _Field.of(p)
    var a = List[UInt64](capacity=len(f.coefficients))
    for i in range(len(f.coefficients)):
        a.append(field.enter(f.coefficients[i].value))
    var b = List[UInt64](capacity=len(g.coefficients))
    for i in range(len(g.coefficients)):
        b.append(field.enter(g.coefficients[i].value))
    return ModInt(field.leave(_resultant_words(a^, b^, field, budget)), f.ring)


def _only_in[C: CoefficientRing](p: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> Bool:
    """Whether every term of p is a constant or a power of the one variable."""
    for t in p.terms:
        budget.charge(1)
        var k = len(t.monomial.pairs)
        if k > 1 or (k == 1 and t.monomial.pairs[0][0] != variable):
            return False
    return True


def resultant_value[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Compute the resultant alone, without the multipliers of its identity.

    Univariate integer inputs use multi-modular reconstruction up to the
    Hadamard bound, univariate inputs over odd prime fields a word-size
    Euclidean recurrence, multivariate integer and odd prime field inputs
    dense evaluation and interpolation when the specialization grid is
    small enough, and all others the subresultant PRS of `resultant`
    without cofactor propagation.

    Parameters:
        C: ZZ or a coefficient field.

    Args:
        f: The first polynomial.
        g: The second polynomial in the same ring and order.
        variable: The elimination variable.
        budget: The shared work allowance.

    Returns:
        The Sylvester-convention resultant, equal to `resultant(f, g, variable).resultant`.

    Raises:
        PolyError: DOMAIN for incompatible inputs or two nonzero constants in the variable, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    var zero = _zero(f, budget)
    if f.is_zero() or g.is_zero():
        return zero^
    # Univariate inputs in the elimination variable go to the dense kernels
    # after one pass over the terms; their degrees come from the dense forms.
    comptime if C == ZZ:
        if _only_in(f, variable, budget) and _only_in(g, variable, budget):
            var fd = DensePoly[ZZ].from_sparse(rebind[SparsePoly[ZZ]](f), budget, variable=variable)
            var gd = DensePoly[ZZ].from_sparse(rebind[SparsePoly[ZZ]](g), budget, variable=variable)
            if fd.degree() > 0 and gd.degree() > 0:
                var value = ZZ(_modular_resultant(fd, gd, budget))
                if value.is_zero():
                    return zero^
                return _constant(f, rebind[C](value), budget)
    elif C == ModInt:
        if rebind[SparsePoly[ModInt]](f).ring.value > 2 and _only_in(f, variable, budget) and _only_in(g, variable, budget):
            var fd = DensePoly[ModInt].from_sparse(rebind[SparsePoly[ModInt]](f), budget, variable=variable)
            var gd = DensePoly[ModInt].from_sparse(rebind[SparsePoly[ModInt]](g), budget, variable=variable)
            if fd.degree() > 0 and gd.degree() > 0:
                var value = _word_resultant(fd, gd, budget)
                if value.is_zero():
                    return zero^
                return _constant(f, rebind[C](value), budget)
    var n = f.degree_in(variable, budget)
    var m = g.degree_in(variable, budget)
    if n == 0 and m == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "resultant identity requires a positive degree")
    if n == 0:
        return f.pow(UInt64(m), budget)
    if m == 0:
        return g.pow(UInt64(n), budget)
    var variables = _variables(f, g, budget)
    if len(variables) > 1:
        comptime if C == ZZ:
            var image = _integer_resultant(rebind[SparsePoly[ZZ]](f), rebind[SparsePoly[ZZ]](g), variable, variables, budget)
            if image:
                return rebind_var[SparsePoly[C]](image.take())
        elif C == ModInt:
            var image = _field_resultant(rebind[SparsePoly[ModInt]](f), rebind[SparsePoly[ModInt]](g), variable, variables, budget)
            if image:
                return rebind_var[SparsePoly[C]](image.take())
    var sequence = _prs(f, g, variable, False, budget, history=False)
    if not len(sequence) or sequence[len(sequence) - 1].polynomial.degree_in(variable, budget) != 0:
        return zero^
    var value = sequence[len(sequence) - 1].principal._copy(budget)
    if n < m and (n % 2) and (m % 2):
        value = value.neg(budget)
    return value^
