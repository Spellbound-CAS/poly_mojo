"""Traces of Gröbner bases over QQ from traced F4 over word primes.

Modulo a prime that divides no numerator or denominator of the inputs, the
traced F4 run (`_f4_trace`) makes the same choices as one over QQ when the
prime is lucky, and its trace is then the image of a rational trace: the
same steps, sources and monomials, with each coefficient the residue of a
rational number. The runs' traces are kept while their structure agrees with
the first one's, their coefficients combined by the Chinese remainder
theorem, and each coefficient recovered as the fraction a/b with |a| and b
at most the square root of half the modulus, then confirmed by the next
prime. Nothing here is trusted: a
candidate is accepted only by an exact replay against the certified basis
(`groebner._qq_f4_trace`), and an unlucky prime or a premature
reconstruction costs only time.
"""
from apn_mojo import Integer, Rational, div_rem_euclid, isqrt, gcd as integer_gcd
from .coefficients import ModInt, Modulus, NoContext, QQ, _work
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .primes import _residue
from ._packed_groebner import _Packed, _Term, _equal
from ._rational_groebner import _integer, _reconstruct
from ._f4_trace import _TracedBasis
from ._f4 import _LANES, _Lanes
from ._integer_f4 import _magnitude


def _admissible[W: Int](inputs: List[List[_Term[QQ, W]]], p: UInt64, mut budget: Budget) raises PolyError -> Bool:
    """Whether p divides no numerator or denominator of the inputs' coefficients."""
    for f in inputs:
        for t in f:
            budget.charge(1)
            if _residue(t.coefficient.value.numerator(), p, budget) == 0 or _residue(t.coefficient.value.denominator(), p, budget) == 0:
                return False
    return True


def _images[W: Int](inputs: List[List[_Term[QQ, W]]], ring: Modulus, mut budget: Budget) raises PolyError -> List[List[_Term[ModInt, W]]]:
    """The inputs modulo the ring's prime, which is admissible for them."""
    var p = ring.value
    _allocation[List[_Term[ModInt, W]]](len(inputs), budget)
    var result = List[List[_Term[ModInt, W]]](capacity=len(inputs))
    for f in inputs:
        _allocation[_Term[ModInt, W]](len(f), budget)
        var terms = List[_Term[ModInt, W]](capacity=len(f))
        for t in f:
            budget.charge(1)
            var numerator = ModInt(_residue(t.coefficient.value.numerator(), p, budget), ring)
            var denominator = ModInt(_residue(t.coefficient.value.denominator(), p, budget), ring)
            terms.append(_Term[ModInt, W](t.monomial, numerator.mul(denominator.inverse(budget), budget)))
        result.append(terms^)
    return result^


def _remainder(limbs: List[UInt64], p: UInt64) -> UInt64:
    """The nonnegative number with these limbs, least significant first, modulo p."""
    var r = UInt128(0)
    for l in range(len(limbs) - 1, -1, -1):
        r = ((r << 64) | UInt128(limbs.unsafe_get(l))) % UInt128(p)
    return UInt64(r)


def _add_scaled(mut x: List[UInt64], m: List[UInt64], t: UInt64):
    """x += t m, for nonnegative numbers in limbs."""
    while len(x) < len(m) + 1:
        x.append(0)
    var carry = UInt128(0)
    for l in range(len(m)):
        var s = UInt128(x.unsafe_get(l)) + UInt128(m.unsafe_get(l)) * UInt128(t) + carry
        x.unsafe_get(l) = UInt64(s & UInt128(UInt64.MAX))
        carry = s >> 64
    var l = len(m)
    while carry:
        if l == len(x):
            x.append(0)
        var s = UInt128(x.unsafe_get(l)) + carry
        x.unsafe_get(l) = UInt64(s & UInt128(UInt64.MAX))
        carry = s >> 64
        l += 1


@always_inline
def _mod_words(v: _Lanes, primes: _Lanes, reciprocals: SIMD[DType.float64, _LANES]) -> _Lanes:
    """Each lane of v, below 2^62, modulo its prime below 2^31, as `_f4._mod_lanes`."""
    var q = (v.cast[DType.float64]() * reciprocals).cast[DType.uint64]()
    var r = v - q * primes
    r = (r >> 63).ne(0).select(r + primes, r)
    return r.ge(primes).select(r - primes, r)


@always_inline
def _half_limbs(r: _Lanes, w: UInt64, primes: _Lanes, reciprocals: SIMD[DType.float64, _LANES], shift: _Lanes) -> _Lanes:
    """(r * 2^64 + w) modulo each lane's prime, by w's two half-limbs; shift is 2^32 modulo each."""
    var high = _mod_words(r * shift + _Lanes(w >> 32), primes, reciprocals)
    return _mod_words(high * shift + _Lanes(w & 0xFFFFFFFF), primes, reciprocals)


def _remainders(limbs: List[UInt64], primes: _Lanes, reciprocals: SIMD[DType.float64, _LANES], shift: _Lanes) -> _Lanes:
    """The number with these limbs, least significant first, modulo each lane's prime; shift is 2^32 modulo each.

    Half-limbs from the top: r * 2^32 + d becomes r * shift + d, below 2^62.
    Each step waits on the last one's reduction, so a long number is cut
    into four blocks whose Horner chains run side by side, joined at the
    end by 2^(64 block).
    """
    var n = len(limbs)
    if n < 16:
        var r = _Lanes(0)
        for l in range(n - 1, -1, -1):
            r = _half_limbs(r, limbs.unsafe_get(l), primes, reciprocals, shift)
        return r
    var block = (n + 3) // 4
    # Only the top block can run past the limbs; its missing limbs are zeros.
    var top = n - 3 * block
    var r0 = _Lanes(0)
    var r1 = _Lanes(0)
    var r2 = _Lanes(0)
    var r3 = _Lanes(0)
    for j in range(block - 1, -1, -1):
        r0 = _half_limbs(r0, limbs.unsafe_get(j), primes, reciprocals, shift)
        r1 = _half_limbs(r1, limbs.unsafe_get(block + j), primes, reciprocals, shift)
        r2 = _half_limbs(r2, limbs.unsafe_get(2 * block + j), primes, reciprocals, shift)
        r3 = _half_limbs(r3, limbs.unsafe_get(3 * block + j) if j < top else UInt64(0), primes, reciprocals, shift)
    # 2^(64 block) is shift^(2 block).
    var power = _Lanes(1)
    var base = shift
    var e = 2 * block
    while e:
        if e & 1:
            power = _mod_words(power * base, primes, reciprocals)
        base = _mod_words(base * base, primes, reciprocals)
        e >>= 1
    var r = _mod_words(r3 * power + r2, primes, reciprocals)
    r = _mod_words(r * power + r1, primes, reciprocals)
    return _mod_words(r * power + r0, primes, reciprocals)


def _lcm_guess(guess: Integer, denominator: Integer, bound: Integer) raises PolyError -> Integer:
    """The least common multiple of a step's denominators so far, or the latest alone once that exceeds the bound."""
    try:
        if guess.sign() <= 0:
            return denominator.copy()
        var l = guess * (denominator // integer_gcd(guess, denominator))
        if l > bound:
            return denominator.copy()
        return l^
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _guessed(value: Integer, denominator: Integer, modulus: Integer, bound: Integer, mut budget: Budget) raises PolyError -> Optional[Rational]:
    """The fraction n/d, in lowest terms, congruent to value with |n| and d at most the bound, for the given d, if there is one: then it is the fraction `_reconstruct` finds."""
    if denominator.sign() <= 0:
        return None
    _work(modulus.magnitude_bit_length() // 8 + 1, denominator.magnitude_bit_length() // 8 + 1, True, budget)
    try:
        if denominator > bound:
            return None
        var n = div_rem_euclid(value * denominator, modulus)[1]
        if n > bound:
            n = n - modulus
            if -n > bound:
                return None
        return Rational(n^, denominator.copy())
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _mul_word(mut x: List[UInt64], w: UInt64):
    """x *= w, for a nonnegative number in limbs."""
    var carry = UInt128(0)
    for l in range(len(x)):
        var s = UInt128(x.unsafe_get(l)) * UInt128(w) + carry
        x.unsafe_get(l) = UInt64(s & UInt128(UInt64.MAX))
        carry = s >> 64
    if carry:
        x.append(UInt64(carry))


def _add_word(mut x: List[UInt64], w: UInt64):
    """x += w, for a nonnegative number in limbs."""
    var carry = UInt128(w)
    var l = 0
    while carry:
        if l == len(x):
            x.append(0)
        var s = UInt128(x.unsafe_get(l)) + carry
        x.unsafe_get(l) = UInt64(s & UInt128(UInt64.MAX))
        carry = s >> 64
        l += 1


def _add_scaled_at(mut x: List[UInt64], m: List[UInt64], t: UInt64, shift: Int):
    """x += t m 2^(64 shift), for nonnegative numbers in limbs."""
    while len(x) < len(m) + shift + 1:
        x.append(0)
    var carry = UInt128(0)
    for l in range(len(m)):
        var s = UInt128(x.unsafe_get(l + shift)) + UInt128(m.unsafe_get(l)) * UInt128(t) + carry
        x.unsafe_get(l + shift) = UInt64(s & UInt128(UInt64.MAX))
        carry = s >> 64
    var l = len(m) + shift
    while carry:
        if l == len(x):
            x.append(0)
        var s = UInt128(x.unsafe_get(l)) + carry
        x.unsafe_get(l) = UInt64(s & UInt128(UInt64.MAX))
        carry = s >> 64
        l += 1


struct _Lifted[W: Int](Movable):
    """The agreed structure of the traces so far and their coefficients combined.

    Each coefficient is open, lifted further by every prime; reconstructed,
    a fraction that the next prime must confirm; or confirmed, a fraction
    no later prime touches. Lifts are kept as limbs and extended in place,
    one word remainder and one multiply-add per prime. Open coefficients
    are reconstructed at prime counts that grow geometrically, so each
    costs a few attempts beyond the primes its height needs.
    """
    var sources: List[Int]
    var monomials: List[_Packed[Self.W]]
    var starts: List[Int]
    var basis_sources: List[Int]
    var values: List[List[UInt64]]
    """Each unconfirmed coefficient's lift below the modulus, in limbs."""
    var fractions: List[QQ]
    var state: List[UInt8]
    """0 open, 1 reconstructed, 2 confirmed."""
    var unconfirmed: Int
    var modulus: List[UInt64]
    var primes: Int
    var next_try: Int

    def __init__(out self, run: _TracedBasis[Self.W], p: UInt64):
        var n = len(run.traced.coefficients)
        self.values = List[List[UInt64]](capacity=n)
        for c in run.traced.coefficients:
            self.values.append([c])
        self.fractions = List[QQ](length=n, fill=QQ.zero(NoContext()))
        self.state = List[UInt8](length=n, fill=0)
        self.unconfirmed = n
        self.sources = run.traced.sources.copy()
        self.monomials = run.traced.monomials.copy()
        self.starts = run.traced.starts.copy()
        self.basis_sources = run.sources.copy()
        self.modulus = [p]
        self.primes = 1
        self.next_try = 2

    def agrees(self, run: _TracedBasis[Self.W]) -> Bool:
        """Whether a run's trace has this structure."""
        ref traced = run.traced
        if traced.starts != self.starts or traced.sources != self.sources or run.sources != self.basis_sources:
            return False
        for k in range(len(self.monomials)):
            if not _equal(traced.monomials[k], self.monomials[k]):
                return False
        return True

    def combine(mut self, run: _TracedBasis[Self.W], p: UInt64, mut budget: Budget) raises PolyError:
        """Fold in a run that agrees."""
        var ring = Modulus(_prime_value=p)
        var inverse = ModInt(_remainder(self.modulus, p), ring).inverse(budget).value
        budget.charge(UInt64(self.unconfirmed * (2 * len(self.modulus) + 1)))
        for k in range(len(self.values)):
            if self.state[k] == 2:
                continue
            var r = run.traced.coefficients[k]
            if self.state[k] == 1:
                # A premature reconstruction is reopened and lifted further.
                ref fraction = self.fractions[k].value
                var denominator = _residue(fraction.denominator(), p, budget)
                if denominator != 0 and ModInt(_residue(fraction.numerator(), p, budget), ring).mul(ModInt(denominator, ring).inverse(budget), budget).value == r:
                    self.state[k] = 2
                    self.unconfirmed -= 1
                    self.values[k] = List[UInt64]()
                    continue
                self.state[k] = 0
            var x = _remainder(self.values[k], p)
            var t = UInt64((UInt128((r + p - x) % p) * UInt128(inverse)) % UInt128(p))
            if t:
                _add_scaled(self.values[k], self.modulus, t)
        var next = List[UInt64]()
        _add_scaled(next, self.modulus, p)
        while len(next) and next[len(next) - 1] == 0:
            _ = next.pop()
        self.modulus = next^
        self.primes += 1
        if self.primes >= self.next_try:
            self.next_try = self.primes + max(1, self.primes // 4)
            var modulus = _integer(self.modulus.unsafe_ptr(), len(self.modulus))
            var bound: Integer
            try:
                bound = isqrt(modulus // Integer(2))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            for k in range(len(self.values)):
                if self.state[k] != 0:
                    continue
                var fraction = _reconstruct(_integer(self.values[k].unsafe_ptr(), len(self.values[k])), modulus, bound, budget)
                if fraction:
                    self.fractions[k] = QQ(fraction.take())
                    self.state[k] = 1

    def combine_lanes(mut self, residues: List[UInt32], live: Int, primes: List[UInt64], mut budget: Budget) raises PolyError:
        """Fold in the first `live` primes of a batch at once, each with this structure: slot k's coefficient modulo primes[l] is residues[_LANES * k + l].

        The lift is the one `combine` reaches prime by prime, by the Chinese
        remainder theorem: t_l = (r_l - x) / M modulo each prime, Garner's
        mixed radix makes them one t below P, their product, and x grows by
        M t in one multiply-add. A reconstructed fraction is confirmed when
        every prime of the batch agrees with it. Reconstruction is tried
        once the prime count reaches the same geometric schedule, on each
        step's open coefficients in order up to the first that fails.
        """
        if live == 0:
            return
        # Per prime: the inverse of M modulo it; Garner's radices, the
        # products of the earlier primes modulo every lane's prime, and the
        # inverse of each prime's own radix.
        var lane_primes = _Lanes(primes[live - 1])
        for l in range(live):
            lane_primes[l] = primes[l]
        var reciprocals = SIMD[DType.float64, _LANES](1.0) / lane_primes.cast[DType.float64]()
        var shifts = (_Lanes(1) << 32) % lane_primes
        var inverses = _remainders(self.modulus, lane_primes, reciprocals, shifts)
        for l in range(live):
            inverses[l] = ModInt(inverses[l], Modulus(_prime_value=primes[l])).inverse(budget).value
        var radices = List[_Lanes](capacity=live)
        var radix = _Lanes(1)
        var garner = _Lanes(1)
        for l in range(live):
            radices.append(radix)
            garner[l] = ModInt(radix[l], Modulus(_prime_value=primes[l])).inverse(budget).value
            radix = _mod_words(radix * (_Lanes(primes[l]) % lane_primes), lane_primes, reciprocals)
        # P in limbs, as Horner's rule builds it.
        var big = List[UInt64]()
        big.append(1)
        for l in range(live):
            _mul_word(big, primes[l])
        budget.charge(UInt64(live * self.unconfirmed * (2 * len(self.modulus) + 1)))
        var mixed = List[UInt64](length=live, fill=0)
        var scaled = List[UInt64]()
        var limbs = List[UInt64]()
        for k in range(len(self.values)):
            if self.state[k] == 2:
                continue
            if self.state[k] == 1:
                # n / d agrees with residue r modulo p when d is a unit and
                # n = r d there: no inverses, every lane at once.
                ref fraction = self.fractions[k].value
                limbs.clear()
                _magnitude(fraction.numerator(), limbs)
                var numerators = _remainders(limbs, lane_primes, reciprocals, shifts)
                if fraction.numerator().sign() < 0:
                    numerators = numerators.eq(0).select(numerators, lane_primes - numerators)
                limbs.clear()
                _magnitude(fraction.denominator(), limbs)
                var denominators = _remainders(limbs, lane_primes, reciprocals, shifts)
                var r = residues.unsafe_ptr().unsafe_offset(_LANES * k).unsafe_load[width=_LANES]().cast[DType.uint64]()
                var products = _mod_words(r * denominators, lane_primes, reciprocals)
                var agrees = True
                for l in range(live):
                    if denominators[l] == 0 or products[l] != numerators[l]:
                        agrees = False
                        break
                if agrees:
                    self.state[k] = 2
                    self.unconfirmed -= 1
                    self.values[k] = List[UInt64]()
                    continue
                self.state[k] = 0
            var xs = _remainders(self.values[k], lane_primes, reciprocals, shifts)
            var r = residues.unsafe_ptr().unsafe_offset(_LANES * k).unsafe_load[width=_LANES]().cast[DType.uint64]()
            var d = r + lane_primes - xs
            d = d.ge(lane_primes).select(d - lane_primes, d)
            var t = _mod_words(d * inverses, lane_primes, reciprocals)
            var nonzero = False
            for l in range(live):
                if t[l]:
                    nonzero = True
            if not nonzero:
                continue
            # Garner: t = m_0 + p_0 (m_1 + p_1 (m_2 + ...)); acc holds, in
            # each lane, the digits so far times their radices.
            var acc = _Lanes(0)
            for l in range(live):
                var p = primes[l]
                var m = ((t[l] + p - acc[l]) % p) * garner[l] % p
                mixed[l] = m
                acc = _mod_words(acc + _Lanes(m) * radices[l], lane_primes, reciprocals)
            scaled.clear()
            scaled.append(mixed[live - 1])
            for l in range(live - 2, -1, -1):
                _mul_word(scaled, primes[l])
                _add_word(scaled, mixed[l])
            for j in range(len(scaled)):
                if scaled[j]:
                    _add_scaled_at(self.values[k], self.modulus, scaled[j], j)
        var next = List[UInt64]()
        for j in range(len(big)):
            if big[j]:
                _add_scaled_at(next, self.modulus, big[j], j)
        while len(next) and next[len(next) - 1] == 0:
            _ = next.pop()
        self.modulus = next^
        self.primes += live
        if self.primes >= self.next_try:
            self.next_try = self.primes + max(1, self.primes // 4)
            var modulus = _integer(self.modulus.unsafe_ptr(), len(self.modulus))
            var bound: Integer
            try:
                bound = isqrt(modulus // Integer(2))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            # A step's coefficients mostly share their denominators' factors:
            # the least common multiple of those found in the step so far is
            # tried first, at the cost of one product.
            var step = 0
            var guess = Integer(0)
            var failed = -1
            for k in range(len(self.values)):
                while self.starts[step + 1] <= k:
                    step += 1
                    guess = Integer(0)
                if step == failed:
                    continue
                if self.state[k] != 0:
                    guess = _lcm_guess(guess, self.fractions[k].value.denominator(), bound)
                    continue
                var value = _integer(self.values[k].unsafe_ptr(), len(self.values[k]))
                var fraction = _guessed(value, guess, modulus, bound, budget)
                if not fraction:
                    fraction = _reconstruct(value, modulus, bound, budget)
                if not fraction:
                    # The step's open coefficients after it wait for the
                    # next attempt: an attempt costs its successes and one
                    # failure per step. Steps differ in height, so a tall
                    # one no longer holds back the short ones after it.
                    failed = step
                    continue
                guess = _lcm_guess(guess, fraction.value().denominator(), bound)
                self.fractions[k] = QQ(fraction.take())
                self.state[k] = 1

    def confirmed(self) -> Bool:
        return self.unconfirmed == 0
