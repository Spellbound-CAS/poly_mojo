"""Multi-modular univariate gcds and resultants over the integers.

Images come from the word-size field kernels in _nmod_euclid. Coefficients are
reconstructed incrementally in symmetric representation: a new prime leaves
a coefficient unchanged exactly when its residue already agrees, so a prime
that changes nothing marks a stable candidate. The gcd is accepted only after
exact division (G-5); the resultant stops at the Hadamard bound.
"""
from apn_mojo import Integer
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum, _product
from .coefficients import ZZ, NoContext
from ._coefficient_ops import _gcd_coefficient, _quotient, _exact
from .dense import DensePoly
from .sparse import SparsePoly
from .primes import prime, _residue
from .coefficient_buffer import _view
from ._ntt import _PRIMES, _PRIME_BITS, _PRIME_COUNT, _series_inverse, _modular_product, _equal_mod, _reconstruct_zz, _bits_zz, _bit_length
from ._nmod_euclid import _Field, _gcd_words, _trim_words


struct _Reconstruction(Movable):
    """Symmetric residues modulo a growing product of primes."""
    var values: List[Integer]
    var modulus: Integer

    def __init__(out self, count: Int, mut budget: Budget) raises PolyError:
        _allocation[Integer](count, budget)
        self.values = List[Integer](length=count, fill=Integer(0))
        self.modulus = Integer(1)

    def headroom(self) -> Int:
        """The bits between the largest value and half the modulus."""
        var top = 0
        for i in range(len(self.values)):
            top = max(top, self.values[i].magnitude_bit_length())
        return self.modulus.magnitude_bit_length() - 1 - top

    def merge(mut self, residues: List[UInt64], p: UInt64, mut budget: Budget) raises PolyError -> Bool:
        """Fold standard residues modulo p into the values; return whether any value changed."""
        var field = _Field.of(p)
        var m = _residue(self.modulus, p, budget)
        if m == 0:
            raise PolyError(PolyErrorKind.DOMAIN, "CRT moduli are not coprime")
        var inverse = field.inverse(field.enter(m))
        var big_p = ZZ(Integer(Int(p)))
        var next_modulus = ZZ(self.modulus.copy()).mul(big_p, budget).value
        var next_half: Integer
        try:
            next_half = next_modulus >> 1
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        var changed = False
        for i in range(len(self.values)):
            budget.charge(1)
            # Zero values with zero residues, most cells of a sparse result, stay zero.
            var zero = self.values[i].magnitude_bit_length() == 0
            if zero and residues[i] == 0:
                continue
            var old = UInt64(0) if zero else _residue(self.values[i], p, budget)
            var difference = field.sub(field.enter(residues[i]), field.enter(old))
            if difference == 0:
                continue
            changed = True
            var correction = field.leave(field.mul(difference, inverse))
            var value = ZZ(self.values[i].copy()).add(ZZ(self.modulus.copy()).mul(ZZ(Integer(Int(correction))), budget), budget).value
            if value > next_half:
                value = ZZ(value^).sub(ZZ(next_modulus.copy()), budget).value
            self.values[i] = value^
        self.modulus = next_modulus^
        return changed


def _residues(f: DensePoly[ZZ], p: UInt64, mut budget: Budget) raises PolyError -> List[UInt64]:
    _allocation[UInt64](len(f.coefficients), budget)
    var out = List[UInt64](capacity=len(f.coefficients))
    for i in range(len(f.coefficients)):
        ref c = f.coefficients[i].value
        var r = _residue(c, p, budget)
        out.append(r)
    return out^


def _entered(values: List[UInt64], field: _Field) -> List[UInt64]:
    var out = List[UInt64](capacity=len(values))
    for v in values:
        out.append(field.enter(v))
    return out^


struct _Words(Movable):
    """The 32-bit magnitude words of a dense integer polynomial's coefficients, flattened once for residues modulo many primes."""
    var words: List[UInt64]
    var offsets: List[Int]
    var negative: List[Bool]
    var longest: Int

    def __init__(out self, f: DensePoly[ZZ], mut budget: Budget) raises PolyError:
        var count = len(f.coefficients)
        _allocation[Int](count + 1, budget)
        self.offsets = List[Int](capacity=count + 1)
        self.negative = List[Bool](capacity=count)
        self.words = List[UInt64]()
        self.longest = 0
        self.offsets.append(0)
        for i in range(count):
            ref c = f.coefficients[i].value
            var n = c._word_count()
            _allocation[UInt64](n, budget)
            for j in range(n):
                self.words.append(UInt64(c._word(j)))
            self.offsets.append(len(self.words))
            self.negative.append(c.sign() < 0)
            self.longest = max(self.longest, n)

    def __init__(out self, p: SparsePoly[ZZ], mut budget: Budget) raises PolyError:
        """The words of a sparse polynomial's coefficients, in term order."""
        var count = len(p.terms)
        _allocation[Int](count + 1, budget)
        self.offsets = List[Int](capacity=count + 1)
        self.negative = List[Bool](capacity=count)
        self.words = List[UInt64]()
        self.longest = 0
        self.offsets.append(0)
        for i in range(count):
            ref c = p.terms[i].coefficient.value
            var n = c._word_count()
            _allocation[UInt64](n, budget)
            for j in range(n):
                self.words.append(UInt64(c._word(j)))
            self.offsets.append(len(self.words))
            self.negative.append(c.sign() < 0)
            self.longest = max(self.longest, n)

    def residues(self, field: _Field, mut budget: Budget, standard: Bool = False) raises PolyError -> List[UInt64]:
        """Every coefficient modulo field.p, in Montgomery form or, when standard, as plain residues.

        With W_j the Montgomery form of 2^(32 j) R, the sum of word_j W_j
        stays below 2^32 p per word, so one reduction yields value R mod p;
        W_j = 2^(32 j) R mod p yields the value itself.
        """
        _charge_size(len(self.words) + self.longest, budget)
        var powers = List[UInt64](capacity=self.longest)
        var step = field.enter(UInt64(1) << 32)
        var w = field.enter(1) if standard else field.r2
        for _ in range(self.longest):
            powers.append(w)
            w = field.mul(w, step)
        _allocation[UInt64](len(self.negative), budget)
        var out = List[UInt64](capacity=len(self.negative))
        var pw = self.words.unsafe_ptr()
        var pp = powers.unsafe_ptr()
        for i in range(len(self.negative)):
            var base = self.offsets[i]
            var count = self.offsets[i + 1] - base
            var word = pw.unsafe_offset(base)
            # Two independent sums; together still below count 2^32 p <= p 2^64.
            var even = UInt128(0)
            var odd = UInt128(0)
            var j = 0
            while j + 1 < count:
                even += UInt128(word.unsafe_offset(j).unsafe_load()) * UInt128(pp.unsafe_offset(j).unsafe_load())
                odd += UInt128(word.unsafe_offset(j + 1).unsafe_load()) * UInt128(pp.unsafe_offset(j + 1).unsafe_load())
                j += 2
            if j < count:
                even += UInt128(word.unsafe_offset(j).unsafe_load()) * UInt128(pp.unsafe_offset(j).unsafe_load())
            var r = field.reduce(even + odd)
            out.append((field.p - r if r else r) if self.negative[i] else r)
        return out^


def _gcd_images(var a: List[UInt64], var b: List[UInt64], gamma: ZZ, field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    """gamma times the monic gcd of the images a, b (Montgomery form) modulo field.p, in standard form."""
    var g = _gcd_words(a^, b^, field, budget)
    var scale = field.enter(_residue(gamma.value, field.p, budget))
    for i in range(len(g)):
        g[i] = field.leave(field.mul(g[i], scale))
    return g^


def _next_gcd_candidate(fd: DensePoly[ZZ], hd: DensePoly[ZZ], mut state: _GcdSearch, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """Advance the prime search to the next stable candidate gamma·g/lc(g), not yet primitive.

    A constant image modulo a prime dividing neither leading coefficient
    proves the gcd constant, so the constant candidate returns at once.
    """
    while True:
        budget.charge(1)
        var p = prime(state.index, budget)
        state.index = _sum(state.index, 1)
        # Both degrees must survive reduction, as required by G-5's proof.
        var field = _Field.of(p)
        var a = state.fw.residues(field, budget)
        var b = state.hw.residues(field, budget)
        if a[len(a) - 1] == 0 or b[len(b) - 1] == 0:
            continue
        var image = _gcd_images(a^, b^, state.gamma, field, budget)
        var degree = len(image) - 1
        if degree > state.best:
            continue
        if degree == 0:
            return DensePoly[ZZ].from_coefficients(fd.variable, NoContext(), [ZZ(Integer(1))], budget)
        if degree < state.best:
            state.best = degree
            state.values = _Reconstruction(degree + 1, budget)
            state.primes = 0
        var changed = state.values.merge(image, p, budget)
        state.primes = _sum(state.primes, 1)
        if state.primes > 1 and not changed:
            _allocation[ZZ](state.best + 1, budget)
            var coefficients = List[ZZ](capacity=state.best + 1)
            for i in range(state.best + 1):
                budget.charge(1)
                coefficients.append(ZZ(state.values.values[i].copy()))
            var dense = DensePoly[ZZ].from_coefficients(fd.variable, NoContext(), coefficients^, budget)
            if dense.degree() == state.best:
                return dense^


struct _GcdSearch(Movable):
    var gamma: ZZ
    var best: Int
    var index: Int
    var primes: Int
    var values: _Reconstruction
    var fw: _Words
    """The operands' coefficient words, flattened once for their images modulo every prime."""
    var hw: _Words

    def __init__(out self, var gamma: ZZ, best: Int, fd: DensePoly[ZZ], hd: DensePoly[ZZ], mut budget: Budget) raises PolyError:
        self.fw = _Words(fd, budget)
        self.hw = _Words(hd, budget)
        self.gamma = gamma^
        self.best = best
        self.index = 0
        self.primes = 0
        self.values = _Reconstruction(0, budget)


def _hadamard_bits(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Int:
    """An upper bound on log2 of the Euclidean norm of the coefficient vector."""
    var top = 0
    for i in range(len(f.coefficients)):
        budget.charge(1)
        top = max(top, f.coefficients[i].value.magnitude_bit_length())
    var count_bits = 0
    while (1 << count_bits) < len(f.coefficients):
        count_bits += 1
    # ||f||_2 <= sqrt(len) max|f_i| < 2^(top + ceil(log2 len) / 2).
    return top + (count_bits + 1) // 2


comptime _SCHOOLBOOK = 2048
"""Quotient-by-divisor term counts up to which certification divides by long division modulo each prime."""


def _dense_content(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> ZZ:
    """The nonnegative gcd of the coefficients."""
    var content = ZZ(Integer(0))
    for i in range(len(f.coefficients)):
        budget.charge(1)
        content = _gcd_coefficient(content, f.coefficients[i], budget)
        if content.value == Integer(1):
            break
    return content^


def _divide_scalar(f: DensePoly[ZZ], c: ZZ, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(f.coefficients), budget)
    var values = List[ZZ](capacity=len(f.coefficients))
    for i in range(len(f.coefficients)):
        budget.charge(1)
        values.append(_exact(f.coefficients[i], c, budget))
    return DensePoly[ZZ]._from_canonical(f.variable, NoContext(), values^)


def _quotient_images(words: _Words, divisor: List[UInt64], inverse: List[UInt64], field: _Field, mut budget: Budget) raises PolyError -> Optional[List[UInt64]]:
    """The quotient of the dividend by the divisor modulo field.p, or None when a remainder is left.

    inverse is the reversed divisor's series inverse to at least the quotient's length.
    """
    var q = field.p
    var a = words.residues(field, budget, standard=True)
    var m = len(a) - len(divisor) + 1
    _allocation[UInt64](2 * m, budget)
    var a_rev = List[UInt64](capacity=m)
    for i in range(m):
        a_rev.append(a[len(a) - 1 - i])
    var head = List[UInt64](capacity=m)
    for i in range(m):
        head.append(inverse[i])
    var q_rev = _modular_product(a_rev, head, m, q, False, budget)
    var quotient = List[UInt64](capacity=m)
    for i in range(m):
        quotient.append(q_rev[m - 1 - i])
    if not _equal_mod(a, quotient, divisor, q, budget):
        return None
    return quotient^


def _schoolbook_quotient(words: _Words, divisor: List[UInt64], field: _Field, mut budget: Budget) raises PolyError -> Optional[List[UInt64]]:
    """The quotient of the dividend by the Montgomery-form divisor modulo field.p, as standard residues, or None when a remainder is left.

    The divisor's leading residue is nonzero. When the remainder vanishes the
    quotient is the one _quotient_images finds, so either may serve a prime.
    """
    var r = words.residues(field, budget)
    var dg = len(divisor)
    var m = len(r) - dg + 1
    _charge_size(_product(m, dg), budget)
    _allocation[UInt64](m, budget)
    var quotient = List[UInt64](length=m, fill=0)
    var lead = field.inverse(divisor[dg - 1])
    var pr = r.unsafe_ptr()
    var pd = divisor.unsafe_ptr()
    # Every index below stays within len(r) = m + dg - 1.
    for k in range(m - 1, -1, -1):
        var c = field.mul(pr.unsafe_offset(k + dg - 1).unsafe_load(), lead)
        quotient[k] = field.leave(c)
        if c == 0:
            continue
        var row = pr.unsafe_offset(k)
        for j in range(dg - 1):
            var cell = row.unsafe_offset(j)
            cell.unsafe_store(field.sub(cell.unsafe_load(), field.mul(c, pd.unsafe_offset(j).unsafe_load())))
    for j in range(dg - 1):
        if r[j] != 0:
            return None
    return quotient^


def _bounded(quotients: List[List[UInt64]], k: Int, count: Int, top_dividend: Int, top_divisor: Int, mut budget: Budget) raises PolyError -> Bool:
    """Whether the symmetric reconstruction Q of the quotient images certifies divisor Q = dividend.

    Q times the divisor agrees with the dividend modulo the product M of the
    first k transform primes; once |dividend| and |divisor| |Q|_1 both lie
    below M/2 the congruence is an identity. Every transform prime exceeds
    2^61, so M > 2^(61 k).
    """
    var candidate = _reconstruct_zz[ZZ](quotients, k, count, budget)
    var top = 0
    for i in range(len(candidate)):
        top = max(top, candidate[i].value.magnitude_bit_length())
    var room = _PRIME_BITS * k - 1
    return top_dividend + 1 <= room and top_divisor + top + _bit_length(count) <= room


def _transform_divides(f: DensePoly[ZZ], fw: _Words, h: DensePoly[ZZ], hw: _Words, g: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Optional[Bool]:
    """Whether g divides both f and h, of degree at least deg g, from quotients modulo transform primes, or None when they cannot decide.

    Each prime's quotients must leave no remainder, and their
    reconstructions must pass the size bound of _bounded; no product is
    formed. Both dividends share the divisor's residues and series inverse.
    """
    var gw = _Words(g, budget)
    var primes = materialize[_PRIMES]()
    var top_f = _bits_zz(_view(f.coefficients))
    var top_h = _bits_zz(_view(h.coefficients))
    var top_g = _bits_zz(_view(g.coefficients))
    var dg = len(g.coefficients)
    var longest = max(len(f.coefficients), len(h.coefficients))
    var span = longest - dg + 1
    var k = (max(top_f, top_h) + _bit_length(longest) + 2 + _PRIME_BITS - 1) // _PRIME_BITS
    var quotients_f = List[List[UInt64]]()
    var quotients_h = List[List[UInt64]]()
    while k <= _PRIME_COUNT:
        for i in range(len(quotients_f), k):
            var field = _Field.of(primes[i])
            var b = gw.residues(field, budget, standard=True)
            if b[dg - 1] == 0:
                return None  # the divisor's degree drops modulo this prime
            var qf: Optional[List[UInt64]]
            var qh: Optional[List[UInt64]]
            if _product(span, dg) <= _SCHOOLBOOK:
                # Small operands: long division costs less than the transforms.
                var bm = gw.residues(field, budget)
                qf = _schoolbook_quotient(fw, bm, field, budget)
                if not qf:
                    return Optional(False)
                qh = _schoolbook_quotient(hw, bm, field, budget)
                if not qh:
                    return Optional(False)
            else:
                var head = List[UInt64](capacity=min(span, dg))
                for j in range(min(span, dg)):
                    head.append(b[dg - 1 - j])
                var inverse = _series_inverse(head, span, primes[i], budget)
                qf = _quotient_images(fw, b, inverse, field, budget)
                if not qf:
                    return Optional(False)
                qh = _quotient_images(hw, b, inverse, field, budget)
                if not qh:
                    return Optional(False)
            quotients_f.append(qf.take())
            quotients_h.append(qh.take())
        if _bounded(quotients_f, k, len(f.coefficients) - dg + 1, top_f, top_g, budget) and _bounded(quotients_h, k, len(h.coefficients) - dg + 1, top_h, top_g, budget):
            return Optional(True)
        if k == _PRIME_COUNT:
            break
        k = min(2 * k, _PRIME_COUNT)
    return None


def _divides_both(f: DensePoly[ZZ], fw: _Words, h: DensePoly[ZZ], hw: _Words, g: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Bool:
    """Whether g divides both f and h, of degree at least deg g: by transform primes when they decide, else by long division."""
    # Long division wins only for short divisors of wide coefficients.
    var wide = max(_bits_zz(_view(f.coefficients)), _bits_zz(_view(h.coefficients))) > 640
    if g.degree() >= 31 or (g.degree() >= 4 and not wide):
        var fast = _transform_divides(f, fw, h, hw, g, budget)
        if fast:
            return fast.value()
    return _dense_divides(f, g, budget) and _dense_divides(h, g, budget)


def _dense_divides(f: DensePoly[ZZ], g: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Bool:
    """Whether g divides f over the integers, by dense exact long division."""
    if f.is_zero():
        return True
    var n = g.degree()
    if f.degree() < n:
        return False
    ref lc = g.coefficients[n]
    var r = f.coefficients.to_list(budget)
    for k in range(f.degree(), n - 1, -1):
        budget.charge(1)
        if r[k].is_zero():
            continue
        var q = _quotient(r[k], lc, budget)
        if not q:
            return False
        ref digit = q.value()
        var base = k - n
        for i in range(n):
            budget.charge(1)
            if not g.coefficients[i].is_zero():
                r[base + i] = r[base + i].sub(digit.mul(g.coefficients[i], budget), budget)
    for i in range(n):
        budget.charge(1)
        if not r[i].is_zero():
            return False
    return True


def _univariate_gcd(f: DensePoly[ZZ], h: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """The positive-leading gcd of two nonzero integer polynomials, integer content included."""
    var cf = _dense_content(f, budget)
    var ch = _dense_content(h, budget)
    var common = _gcd_coefficient(cf, ch, budget)
    var pf = _divide_scalar(f, cf, budget)
    var ph = _divide_scalar(h, ch, budget)
    var constant = DensePoly[ZZ].from_coefficients(f.variable, NoContext(), [common.copy()], budget)
    if pf.degree() == 0 or ph.degree() == 0:
        return constant^
    var gamma = _gcd_coefficient(pf.coefficients[pf.degree()], ph.coefficients[ph.degree()], budget)
    var fw = _Words(pf, budget)
    var hw = _Words(ph, budget)
    var best = min(pf.degree(), ph.degree()) + 1
    var values = _Reconstruction(0, budget)
    var used = 0
    var index = 0
    while True:
        budget.charge(1)
        var p = prime(index, budget)
        index = _sum(index, 1)
        var field = _Field.of(p)
        var a = fw.residues(field, budget)
        var b = hw.residues(field, budget)
        # Both degrees must survive reduction, as required by G-5's proof.
        if a[len(a) - 1] == 0 or b[len(b) - 1] == 0:
            continue
        var image = _gcd_words(a^, b^, field, budget)
        var degree = len(image) - 1
        if degree > best:
            continue
        if degree == 0:
            return constant^
        if degree < best:
            best = degree
            values = _Reconstruction(degree + 1, budget)
            used = 0
        var scale = field.enter(_residue(gamma.value, p, budget))
        for i in range(len(image)):
            image[i] = field.leave(field.mul(image[i], scale))
        var changed = values.merge(image, p, budget)
        used = _sum(used, 1)
        # Exact division decides; it is tried once the images agree or every
        # reconstructed coefficient sits 16 bits inside the symmetric range,
        # which an incomplete reconstruction of two or more random
        # coefficients does with probability below 2^-32.
        if (used > 1 and not changed) or values.headroom() >= 16:
            _allocation[ZZ](best + 1, budget)
            var coefficients = List[ZZ](capacity=best + 1)
            for i in range(best + 1):
                budget.charge(1)
                coefficients.append(ZZ(values.values[i].copy()))
            var candidate = DensePoly[ZZ].from_coefficients(pf.variable, NoContext(), coefficients^, budget)
            if candidate.degree() != best:
                continue
            var primitive = _divide_scalar(candidate, _dense_content(candidate, budget), budget)
            if primitive.coefficients[primitive.degree()].value.sign() < 0:
                primitive = primitive.scale(ZZ(Integer(-1)), budget)
            if _divides_both(pf, fw, ph, hw, primitive, budget):
                return primitive.scale(common, budget)
