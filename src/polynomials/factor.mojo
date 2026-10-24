"""Univariate integer and rational factorization by Hensel–Zassenhaus (Z-6)."""
from apn_mojo import Integer
from .coefficients import CoefficientRing, ZZ, QQ, ModInt, Modulus, NoContext, _pow, _prime
from .dense import DensePoly
from ._subsets import _Recombination
from ._coefficient_ops import _copy as _copy_coefficient
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _charge_size, _product
from .certificates import Factorization
from .squarefree import squarefree
from .finite_field import factor_mod_p, _ModularPattern, _factor_pattern, _pattern_factors
from .primes import prime, _residue, _symmetric
from ._dense_ops import _copy, _constant, _monic, _derivative, _gcd
from ._images import _integer_dense as _integer, _rational_dense as _rational, _mod_image
from ._coefficient_ops import _quotient
from ._hensel import _HenselState, _hensel, _hensel_start, _hensel_extend, _hensel_resume, _truncate, _trial_divide
from ._integer_ops import _content, _primitive, _divide_scalar
from ._lattice_factor import _van_hoeij_lifted, _lattice_precision, _may_divide
from ._gcd import _modular_primitive
from ._univariate_modular import _GcdSearch, _next_gcd_candidate
from ._coefficient_ops import _gcd_coefficient
from .monomial import MonomialOrder
from . import _nmod
from ._cyclotomic import _binomial_factors
from ._ntt import _ntt_divide_exact
from .coefficient_buffer import _view


def _bound(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> ZZ:
    # A factor h of degree d with cofactor q reconstructs as lc(q) h, and
    # |lc(q) h_j| <= binomial(d, j) |lc(q)| M(h) <= binomial(d, j) M(q) M(h)
    # = binomial(d, j) M(f) <= 2^d ||f||_2 <= 2^d ||f||_1 for the Mahler
    # measure M, so lc(f) does not enter. Recombination always reconstructs
    # the side whose degree is at most half the original degree, and twice
    # the bound keeps the symmetric residues exact.
    var norm = ZZ.zero(NoContext())
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        norm = norm.add(c.neg(budget) if c.value.sign() < 0 else _copy_coefficient(c, budget), budget)
    var power = _pow(ZZ(Integer(2)), UInt64(f.degree() // 2 + 1), budget)
    return norm.mul(power, budget)


def _subset_degrees(degrees: List[Int], total: Int, mut budget: Budget) raises PolyError -> List[Bool]:
    """Which degrees 0..total are sums of a subset of the modular factor degrees."""
    _allocation[Bool](_sum(total, 1), budget)
    var reachable = List[Bool](length=total + 1, fill=False)
    reachable[0] = True
    var top = 0
    for d in degrees:
        _charge_size(_sum(top, 1), budget)
        var i = top
        while i >= 0:
            if reachable[i] and i + d <= total:
                reachable[i + d] = True
            i -= 1
        top = min(total, top + d)
    return reachable^


def _recombine(f: DensePoly[ZZ], lifted: List[DensePoly[ZZ]], modulus: ZZ, allowed: List[Bool], mut budget: Budget, prime: UInt64 = 0) raises PolyError -> List[DensePoly[ZZ]]:
    var partial = _recombine_until(f, lifted, modulus, allowed, 0, budget, prime)
    var result = List[DensePoly[ZZ]]()
    swap(result, partial.factors)
    if partial.remaining.degree() > 0:
        _allocation[DensePoly[ZZ]](1, budget)
        result.append(_copy(partial.remaining, budget))
    return result^


@fieldwise_init
struct _Partial(Movable):
    """Zassenhaus recombination so far: the factors found, the cofactor left,
    and the indices of the lifted factors it still holds when the search
    stopped at its subset limit (empty when it finished)."""
    var factors: List[DensePoly[ZZ]]
    var remaining: DensePoly[ZZ]
    var active: List[Int]


def _recombine_until(f: DensePoly[ZZ], lifted: List[DensePoly[ZZ]], modulus: ZZ, allowed: List[Bool], limit: Int, mut budget: Budget, prime: UInt64 = 0) raises PolyError -> _Partial:
    """Zassenhaus recombination; with a positive limit it stops before any
    subset larger than it while _LATTICE_SUBSETS lifted factors remain."""
    # `allowed[d]` is False for degrees that no factor of f can have, because
    # some prime's factor degrees cannot sum to d.
    var search = _Recombination(len(lifted), budget)
    var remaining = _copy(f, budget)
    var result = List[DensePoly[ZZ]]()
    var fractions = List[Float64]()
    var total = Float64(0)
    var tolerance = Float64(1)
    var fresh = False
    # Below 2^62 the constant-term product modulo M takes words: its symmetric
    # value must not exceed |lc(f) f(0)| in magnitude to divide it, which
    # rejects nearly every wrong subset before the exact test.
    # When M = p^k with k <= 4 and p < 2^28, residues modulo M are k base-p
    # digits whose products stay in words: the constant-term product of a
    # subset then costs word operations, and its symmetric value must not
    # exceed |lc(f) f(0)| in magnitude to divide it.
    var power = _small_power(modulus.value, prime)
    var small = power > 0
    var constants = List[_Digits]()
    var lead_word = _Digits()
    var target_word = UInt128(0)
    var m_word = UInt128(0)
    var bounded = False
    # The same holds at x = 1 and x = -1: lc(q) h(x) divides lc(f) f(x).
    # Those digits are prepared once a subset first passes the constant term.
    var at = List[List[_Digits]]()
    var at_target = List[UInt128](length=_POINTS, fill=0)
    var at_bounded = List[Bool](length=_POINTS, fill=False)
    var at_fresh = False
    var chosen = List[Bool]()
    var reciprocal = UInt64(0)
    if small:
        m_word = _as_u128(modulus.value)
        reciprocal = UInt64(0xFFFFFFFFFFFFFFFF) // prime
        _allocation[_Digits](len(lifted), budget)
        for i in range(len(lifted)):
            constants.append(_Digits(lifted[i].coefficients[0].value, prime, power, budget))
    while search.ready():
        budget.charge(1)
        if limit > 0 and len(search.selection.indices) > limit and len(search.active) >= _LATTICE_SUBSETS:
            return _Partial(result^, remaining^, search.active.copy())
        if not fresh:
            fractions = _subleading_fractions(remaining, lifted, modulus, budget)
            total = 0
            for i in range(len(search.active)):
                budget.charge(1)
                total += fractions[search.active[i]]
            tolerance = _subleading_tolerance(remaining, modulus, budget)
            if small:
                lead_word = _Digits(remaining.coefficients[remaining.degree()].value, prime, power, budget)
                var top = _copy_coefficient(remaining.coefficients[remaining.degree()], budget)
                var target = top.mul(remaining.coefficients[0], budget).value
                bounded = target.sign() != 0 and target.magnitude_bit_length() <= 126
                if bounded:
                    target_word = _as_u128(target)
                at_fresh = False
            fresh = True
        var degree = 0
        for i in search.selection.indices:
            budget.charge(1)
            degree = _sum(degree, lifted[search.active[i]].degree())
        var complement = degree > remaining.degree() // 2
        _allocation[Bool](len(search.active), budget)
        chosen.resize(len(search.active), complement)
        for i in range(len(chosen)):
            chosen[i] = complement
        for i in search.selection.indices:
            budget.charge(1)
            chosen[i] = not complement
        var chosen_degree = remaining.degree() - degree if complement else degree
        if chosen_degree < len(allowed) and not allowed[chosen_degree]:
            search.advance(budget)
            continue
        # The x^(k-1) coefficient of lc(q) * h is lc(f) times the sum of the
        # roots of h, so its residue lies within the tolerance of a multiple
        # of the modulus. Compare the fractional parts first.
        if tolerance < 0.5:
            var part = Float64(0)
            for i in search.selection.indices:
                budget.charge(1)
                part += fractions[search.active[i]]
            if complement:
                part = total - part
            part -= _floor(part)
            if min(part, 1 - part) > tolerance:
                search.advance(budget)
                continue
        # A true factor h with cofactor q reconstructs as lc(q) * h exactly, so
        # its constant term lc(q) * h(0) divides lc(f) * f(0). Test that first.
        if bounded:
            var t = lead_word.copy()
            for i in range(len(search.active)):
                if chosen[i]:
                    t = t.times(constants[search.active[i]], prime, reciprocal, power)
            if not _divides_word(t.value(prime, power), target_word, m_word):
                search.advance(budget)
                continue
        if small and len(lifted) >= _POINT_FACTORS:
            if len(at) == 0:
                _allocation[_Digits](_product(len(lifted), _POINTS), budget)
                for point in range(_POINTS):
                    var values = List[_Digits](capacity=len(lifted))
                    for i in range(len(lifted)):
                        values.append(_Digits(_at_unit(lifted[i], point, budget).value, prime, power, budget))
                    at.append(values^)
            if not at_fresh:
                var top = _copy_coefficient(remaining.coefficients[remaining.degree()], budget)
                for point in range(_POINTS):
                    var value = top.mul(_at_unit(remaining, point, budget), budget).value
                    at_bounded[point] = value.sign() != 0 and value.magnitude_bit_length() <= 126
                    if at_bounded[point]:
                        at_target[point] = _as_u128(value)
                at_fresh = True
            var rejected = False
            for point in range(_POINTS):
                if not at_bounded[point]:
                    continue
                var t = lead_word.copy()
                for i in range(len(search.active)):
                    if chosen[i]:
                        t = t.times(at[point][search.active[i]], prime, reciprocal, power)
                if not _divides_word(t.value(prime, power), at_target[point], m_word):
                    rejected = True
                    break
            if rejected:
                search.advance(budget)
                continue
        var lead = _copy_coefficient(remaining.coefficients[remaining.degree()], budget)
        if not remaining.coefficients[0].is_zero():
            var trailing = _copy_coefficient(lead, budget)
            for i in range(len(search.active)):
                budget.charge(1)
                if chosen[i]:
                    trailing = ZZ(_symmetric(trailing.mul(lifted[search.active[i]].coefficients[0], budget).value, modulus.value, budget))
            var target = lead.mul(remaining.coefficients[0], budget)
            if trailing.is_zero() or not _quotient(target, trailing, budget):
                search.advance(budget)
                continue
        var candidate = _constant(remaining, lead, budget)
        for i in range(len(search.active)):
            budget.charge(1)
            if chosen[i]:
                candidate = _truncate(candidate.mul(lifted[search.active[i]], budget), modulus, budget)
        candidate = _primitive(candidate, budget)
        var quotient = Optional[DensePoly[ZZ]](None)
        if candidate.degree() > 0 and candidate.degree() < remaining.degree():
            if candidate.coefficients[0].is_zero() or Bool(_quotient(remaining.coefficients[0], candidate.coefficients[0], budget)):
                quotient = _trial_divide(remaining, candidate, budget)
        if quotient:
            _allocation[DensePoly[ZZ]](1, budget)
            if complement:
                # The selected (smaller-cardinality) side is irreducible:
                # every proper subset was tested earlier. Reconstructing
                # its lower-degree complement only supplies the quotient.
                result.append(quotient.take())
                remaining = candidate^
            else:
                result.append(candidate^)
                remaining = quotient.take()
            search.accept(budget)
            fresh = False
        else:
            search.advance(budget)
    return _Partial(result^, remaining^, List[Int]())


comptime _DIGIT_PRIME = UInt64(1) << 28
comptime _POINTS = 2
"""Evaluation points x = 1 and x = -1 of the word-sized divisibility filter."""
comptime _POINT_FACTORS = 8
"""Lifted factor counts from which the evaluation filter pays for its digits."""


def _at_unit(g: DensePoly[ZZ], point: Int, mut budget: Budget) raises PolyError -> ZZ:
    """g(1) for point 0, g(-1) for point 1."""
    var total = ZZ.zero(NoContext())
    for i in range(len(g.coefficients)):
        budget.charge(1)
        if point == 1 and i % 2 == 1:
            total = total.sub(g.coefficients[i], budget)
        else:
            total = total.add(g.coefficients[i], budget)
    return total^


@always_inline
def _divides_word(v: UInt128, target: UInt128, m: UInt128) -> Bool:
    """Whether the residue v modulo m can be a divisor of a value of magnitude target.

    A true divisor has magnitude at most target; below m / 2 its symmetric
    residue is itself, and below m it is nonzero.
    """
    var magnitude = v if v <= m / 2 else m - v
    if magnitude > target:
        return False
    return v != 0 or target >= m


comptime _DIGITS = 4


struct _Digits(ImplicitlyCopyable, Movable):
    """A residue modulo p^k (k <= 4, p < 2^28) as k base-p digits, lowest first."""
    var d: SIMD[DType.uint64, _DIGITS]

    def __init__(out self):
        self.d = SIMD[DType.uint64, _DIGITS](0)

    def __init__(out self, value: Integer, p: UInt64, k: Int, mut budget: Budget) raises PolyError:
        self.d = SIMD[DType.uint64, _DIGITS](0)
        var modulus = Integer(1)
        try:
            for _ in range(k):
                modulus = modulus * Integer(Int(p))
            var r = value % modulus
            if r.sign() < 0:
                r = r + modulus
            # r < p^k < 2^112 fits 128 bits; peel its digits there.
            var v = _as_u128(r)
            for i in range(k):
                self.d[i] = UInt64(v % UInt128(p))
                v = v / UInt128(p)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    @always_inline
    def times(self, other: Self, p: UInt64, reciprocal: UInt64, k: Int) -> Self:
        """The product modulo p^k: digit convolution below p^k, sums below 4 p^2 < 2^58, then carries.

        reciprocal = floor((2^64 - 1) / p) estimates each carry within one.
        """
        var out = Self()
        var carry = UInt64(0)
        for t in range(k):
            var c = carry
            for i in range(t + 1):
                c += self.d[i] * other.d[t - i]
            var q = UInt64((UInt128(c) * UInt128(reciprocal)) >> 64)
            var r = c - q * p
            if r >= p:
                r -= p
                q += 1
            out.d[t] = r
            carry = q
        return out

    @always_inline
    def value(self, p: UInt64, k: Int) -> UInt128:
        var v = UInt128(0)
        for i in range(k - 1, -1, -1):
            v = v * UInt128(p) + UInt128(self.d[i])
        return v


def _small_power(modulus: Integer, p: UInt64) -> Int:
    """k when modulus = p^k with 1 <= k <= 4 and 2 < p < 2^28, else 0."""
    if p < 3 or p >= _DIGIT_PRIME:
        return 0
    var power = Integer(1)
    try:
        for k in range(1, _DIGITS + 1):
            power = power * Integer(Int(p))
            if power == modulus:
                return k
    except:
        return 0
    return 0


def _as_u128(value: Integer) -> UInt128:
    var out = UInt128(0)
    for w in range(value._word_count() - 1, -1, -1):
        out = (out << 32) | UInt128(value._word(w))
    return out


def _floor(x: Float64) -> Float64:
    var t = Float64(Int(x))
    return t - 1 if t > x else t


def _fraction(r: Integer, modulus: Integer) raises PolyError -> Float64:
    """r / modulus for 0 <= r < modulus, within 2^-52 of the exact ratio."""
    try:
        var shift = max(0, modulus.magnitude_bit_length() - 62)
        return Float64(Int(r >> shift)) / Float64(Int(modulus >> shift))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _subleading_fractions(f: DensePoly[ZZ], lifted: List[DensePoly[ZZ]], modulus: ZZ, mut budget: Budget) raises PolyError -> List[Float64]:
    """lc(f) * (subleading coefficient of each monic lifted factor) mod the
    modulus, as a fraction of the modulus."""
    ref lead = f.coefficients[f.degree()]
    _allocation[Float64](len(lifted), budget)
    var result = List[Float64](capacity=len(lifted))
    for g in lifted:
        budget.charge(1)
        var r = _symmetric(lead.mul(g.coefficients[g.degree() - 1], budget).value, modulus.value, budget)
        if r.sign() < 0:
            r = ZZ(r).add(modulus, budget).value
        result.append(_fraction(r, modulus.value))
    return result^


def _subleading_tolerance(f: DensePoly[ZZ], modulus: ZZ, mut budget: Budget) raises PolyError -> Float64:
    """A bound on |lc(f) * (sum of the roots of any factor)| over the modulus,
    plus the rounding slack of summing fractions.

    The roots of f with |a| >= 1 have product at most M(f)/|lc(f)| <= ||f||_1
    / |lc(f)|, and numbers at least 1 sum to at most their count less one plus
    their product, so the sum of |a| over all roots is at most (d - 1) +
    ||f||_1 / |lc(f)|."""
    ref lead = f.coefficients[f.degree()]
    var bound = ZZ.zero(NoContext())
    for i in range(len(f.coefficients)):
        budget.charge(1)
        ref c = f.coefficients[i]
        bound = bound.add(c.neg(budget) if c.value.sign() < 0 else _copy_coefficient(c, budget), budget)
    var magnitude = lead.neg(budget) if lead.value.sign() < 0 else _copy_coefficient(lead, budget)
    bound = bound.add(magnitude.mul(ZZ(Integer(f.degree() - 1)), budget), budget)
    if bound.value.magnitude_bit_length() >= modulus.value.magnitude_bit_length():
        return 1
    return _fraction(bound.value, modulus.value) + Float64(2.0) ** -40


comptime _PRIME_TRIALS = 3
"""Good primes whose modular factorizations are compared before lifting."""
comptime _LATTICE_FACTORS = 20
"""Modular factors from which recombination takes the lattice."""
comptime _LATTICE_SUBSETS = 12
"""Modular factors from which Zassenhaus recombination tries only subsets up
to _SMALL_SUBSETS before the lattice takes the rest."""
comptime _SMALL_SUBSETS = 2
comptime _LATTICE_TRIAL_DEGREE = 128
"""Degree from which one prime with the lattice ahead ends the trials:
measured, a trial there costs more than its fewer factors save."""
comptime _FEW_TRIALS_DEGREE = 10
"""Degrees up to which two good primes are compared; lifting and recombination there cost less than a third trial."""


comptime _SMALL_PRIME_COUNT = 48
comptime _SMALL_PRIMES: Array[UInt64, _SMALL_PRIME_COUNT] = [
    268435399, 268435367, 268435361, 268435337, 268435331, 268435313, 268435291, 268435273, 268435243, 268435183, 268435171, 268435157, 268435147, 268435133, 268435129, 268435121, 268435109, 268435091, 268435067, 268435043, 268435039, 268435033, 268435019, 268435009, 268435007, 268434997, 268434979, 268434977, 268434961, 268434949, 268434941, 268434937, 268434857, 268434841, 268434827, 268434821, 268434787, 268434781, 268434779, 268434773, 268434731, 268434721, 268434713, 268434707, 268434703, 268434697, 268434659, 268434623
]
"""The largest primes below 2^28, descending: the sequence _next_small_prime walks."""


def _next_small_prime(previous: UInt64, mut budget: Budget) raises PolyError -> UInt64:
    """The largest prime below `previous`; starting from 2^28 this walks the
    descending primes whose residue products take the SIMD path, and
    quadratic lifting makes their size irrelevant. The first steps come from
    a table instead of primality tests."""
    budget.charge(1)
    var table = materialize[_SMALL_PRIMES]()
    if previous <= (UInt64(1) << 28) and previous > table[_SMALL_PRIME_COUNT - 1]:
        for i in range(_SMALL_PRIME_COUNT):
            if table[i] < previous:
                return table[i]
    var candidate = previous - 1
    if candidate % 2 == 0:
        candidate -= 1
    while True:
        budget.charge(1)
        if candidate < 3:
            raise PolyError(PolyErrorKind.CAPACITY, "small prime sequence exhausted")
        if _prime(candidate, budget):
            return candidate
        candidate -= 2


def _squarefree_factors(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    if f.degree() <= 1:
        _allocation[DensePoly[ZZ]](1, budget)
        return [_copy(f, budget)]
    # x^n - 1 and x^n + 1 split into known cyclotomic irreducibles.
    var binomial = _binomial_factors(f, budget)
    if binomial:
        return binomial.take()
    # Factor modulo several good primes. Each one's factor degrees bound the
    # degrees of true factors; lifting uses the prime with fewest factors.
    var p = UInt64(1) << 28
    var trials = 0
    _allocation[Bool](_sum(f.degree(), 1), budget)
    var allowed = List[Bool](length=f.degree() + 1, fill=True)
    # Trial primes need only their factor degrees; the chosen prime alone is split.
    var best_count = 0
    var best_image = Optional[DensePoly[ModInt]](None)
    var best_pattern = Optional[_ModularPattern](None)
    var prime_trials = 2 if f.degree() <= _FEW_TRIALS_DEGREE else _PRIME_TRIALS
    while trials < prime_trials:
        budget.charge(1)
        p = _next_small_prime(p, budget)
        if _residue(f.coefficients[f.degree()].value, p, budget) == 0:
            continue
        var ctx = Modulus(_prime_value=p)
        var image = _monic(_mod_image(f, ctx, budget), budget)
        if _gcd(image, _derivative(image, budget), budget).degree() > 0:
            continue
        trials += 1
        var pattern = _factor_pattern(image, budget)
        var degrees = pattern.degrees.copy()
        if len(degrees) == 1:
            _allocation[DensePoly[ZZ]](1, budget)
            return [_copy(f, budget)]
        var reachable = _subset_degrees(degrees, f.degree(), budget)
        var proper = False
        for d in range(f.degree() + 1):
            allowed[d] = allowed[d] and reachable[d]
            if d > 0 and d < f.degree() and allowed[d]:
                proper = True
        if not proper:
            _allocation[DensePoly[ZZ]](1, budget)
            return [_copy(f, budget)]
        var lattice = len(degrees) >= _LATTICE_SUBSETS and f.degree() >= _LATTICE_TRIAL_DEGREE
        if lattice:
            # A trial costs about degree^2 word operations and lifting about
            # precision * degree per factor, so a prime pays for the factors
            # it saves only while the lattice's precision exceeds the degree.
            var precision = _lattice_precision(f, len(degrees), ZZ(Integer(Int(p))), budget)
            lattice = Bool(precision) and precision.value().value.magnitude_bit_length() < f.degree()
        if best_count == 0 or len(degrees) < best_count:
            best_count = len(degrees)
            best_image = image^
            best_pattern = pattern^
        # The lattice is ahead, and another prime would trim only a few
        # factors from it, for less than its own trial costs.
        if lattice:
            break
    var best = _pattern_factors(best_pattern.value(), best_image.value(), budget)
    # Measured crossover: exact LLL loses on 8/16-way splitting but
    # wins at 20. Peel singleton groups before paying for the lattice.
    if len(best) >= _LATTICE_FACTORS:
        return _large_recombine(f, best, allowed, budget)
    var lifting = _hensel(f, best, _bound(f, budget), budget)
    if len(best) < _LATTICE_SUBSETS:
        return _recombine(f, lifting[0], lifting[1], allowed, budget, prime=best[0].ring.value)
    # Small subsets find the factors of a product cheaply; an irreducible
    # remainder over many modular factors would need every subset up to half
    # of them, which the lattice replaces by one reduction.
    var partial = _recombine_until(f, lifting[0], lifting[1], allowed, _SMALL_SUBSETS, budget, prime=best[0].ring.value)
    var result = List[DensePoly[ZZ]]()
    swap(result, partial.factors)
    if len(partial.active) == 0:
        if partial.remaining.degree() > 0:
            _allocation[DensePoly[ZZ]](1, budget)
            result.append(_copy(partial.remaining, budget))
        return result^
    var active_modular = List[DensePoly[ModInt]]()
    var active_lifted = List[DensePoly[ZZ]]()
    for i in partial.active:
        budget.charge(1)
        _allocation[DensePoly[ZZ]](1, budget)
        active_modular.append(_copy(best[i], budget))
        active_lifted.append(_copy(lifting[0][i], budget))
    var resumed = _hensel_resume(active_lifted^, _copy_coefficient(lifting[1], budget))
    for g in _van_hoeij_lifted(partial.remaining, active_modular, _bound(partial.remaining, budget), resumed, budget):
        budget.charge(1)
        _allocation[DensePoly[ZZ]](1, budget)
        result.append(_copy(g, budget))
    return result^


def _large_recombine(f: DensePoly[ZZ], modular: List[DensePoly[ModInt]], mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    _allocation[Bool](_sum(f.degree(), 1), budget)
    return _large_recombine(f, modular, List[Bool](length=f.degree() + 1, fill=True), budget)


def _large_recombine(f: DensePoly[ZZ], modular: List[DensePoly[ModInt]], allowed: List[Bool], mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    var state = _hensel_start(modular, budget)
    # Peel and run the lattice at the lattice's own precision; factors that
    # need more are lifted further only when reconstruction asks for it.
    var bound = _bound(f, budget)
    var precision = _lattice_precision(f, len(modular), ZZ(Integer(Int(modular[0].ring.value))), budget)
    # Lift to the lattice's first precision even past the bound: a lattice
    # after peeling would otherwise lift again in steps that miss it. The
    # precision is a power of p, which the lattice may use as its modulus,
    # and the lift goes strictly past its bound.
    var first = _copy_coefficient(bound, budget)
    if precision:
        first = precision.value().sub(ZZ.one(NoContext()), budget)
    _hensel_extend(f, modular, first, state, budget)
    var lifting = (state.lifted.copy(), _copy_coefficient(state.modulus, budget))
    var remaining = _copy(f, budget)
    var result = List[DensePoly[ZZ]]()
    var active_modular = List[DensePoly[ModInt]]()
    var active_lifted = List[DensePoly[ZZ]]()
    for i in range(len(modular)):
        budget.charge(1)
        var candidate = _primitive(_truncate(lifting[0][i].scale(remaining.coefficients[remaining.degree()], budget), lifting[1], budget), budget)
        var quotient = Optional[DensePoly[ZZ]](None)
        if _may_divide(remaining, candidate, budget):
            quotient = _trial_divide(remaining, candidate, budget)
        if candidate.degree() > 0 and candidate.degree() == lifting[0][i].degree() and quotient:
            remaining = quotient.take()
            _allocation[DensePoly[ZZ]](1, budget)
            result.append(candidate^)
        else:
            _allocation[DensePoly[ModInt]](1, budget)
            _allocation[DensePoly[ZZ]](1, budget)
            active_modular.append(_copy(modular[i], budget))
            active_lifted.append(_copy(lifting[0][i], budget))
    if remaining.degree() > 0:
        var rest: List[DensePoly[ZZ]]
        if len(active_modular) >= _LATTICE_FACTORS:
            if len(active_modular) == len(modular):
                # Nothing peeled: keep the factor tree for further lifting.
                rest = _van_hoeij_lifted(remaining, active_modular, _bound(remaining, budget), state, budget)
            else:
                var resumed = _hensel_resume(active_lifted^, _copy_coefficient(lifting[1], budget))
                rest = _van_hoeij_lifted(remaining, active_modular, _bound(remaining, budget), resumed, budget)
        else:
            var remaining_bound = _bound(remaining, budget)
            if lifting[1].value > remaining_bound.value:
                rest = _recombine(remaining, active_lifted, lifting[1], allowed, budget)
            else:
                var full = _hensel(remaining, active_modular, remaining_bound, budget)
                rest = _recombine(remaining, full[0], full[1], allowed, budget)
        for g in rest:
            budget.charge(1)
            _allocation[DensePoly[ZZ]](1, budget)
            result.append(_copy(g, budget))
    return result^


def _less[C: CoefficientRing](a: DensePoly[C], b: DensePoly[C], mut budget: Budget) raises PolyError -> Bool:
    budget.charge(1)
    if a.degree() != b.degree():
        return a.degree() < b.degree()
    for i in range(a.degree(), -1, -1):
        budget.charge(1)
        var delta = a.coefficients[i].sub(b.coefficients[i], budget)
        if not delta.is_zero():
            comptime if C == ZZ:
                return rebind[ZZ](delta).value.sign() < 0
            elif C == QQ:
                return rebind[QQ](delta).value.sign() < 0
    return False


def _sort[C: CoefficientRing](mut factors: List[Tuple[DensePoly[C], Int]], template: DensePoly[C], mut budget: Budget) raises PolyError:
    for i in range(1, len(factors)):
        budget.charge(1)
        var current = (DensePoly[C](template.variable, template.ring.copy()), 0)
        swap(current, factors[i])
        var j = i
        while j > 0 and _less(current[0], factors[j - 1][0], budget):
            var moved = (DensePoly[C](template.variable, template.ring.copy()), 0)
            swap(moved, factors[j - 1])
            swap(moved, factors[j])
            j -= 1
        swap(current, factors[j])


def _modular_squarefree(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Bool:
    # A prime that keeps the degree and leaves f coprime to f' proves f
    # squarefree over QQ: a repeated primitive factor would survive reduction.
    var p = UInt64(1) << 28
    for _ in range(2):
        budget.charge(1)
        p = _next_small_prime(p, budget)
        if _residue(f.coefficients[f.degree()].value, p, budget) == 0:
            continue
        var field = _nmod._field(p)
        _allocation[UInt64](len(f.coefficients), budget)
        var image = List[UInt64](capacity=len(f.coefficients))
        for i in range(len(f.coefficients)):
            image.append(_nmod._to(_residue(f.coefficients[i].value, p, budget), field))
        if len(_nmod._gcd(image, _nmod._derivative(image, field, budget), field, budget)) == 1:
            return True
    return False


def _integer_gcd(a: DensePoly[ZZ], b: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    # The primitive gcd with positive leading coefficient, by modular images.
    if b.is_zero():
        return _primitive(a, budget)
    if a.is_zero():
        return _primitive(b, budget)
    if a.degree() <= 0 or b.degree() <= 0:
        return _constant(a, ZZ.one(NoContext()), budget)
    var order = MonomialOrder()
    var g = _modular_primitive(_primitive(a, budget).to_sparse(order, budget), _primitive(b, budget).to_sparse(order, budget), a.variable, budget)
    return _primitive(DensePoly[ZZ].from_sparse(g, budget, variable=a.variable), budget)


comptime _TRANSFORM_QUOTIENT = 16
"""Divisor degrees from which Yun's exact quotients go through transforms."""


def _divide_integer(a: DensePoly[ZZ], b: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Optional[DensePoly[ZZ]]:
    """a / b when b divides a over ZZ, else None."""
    if b.degree() >= _TRANSFORM_QUOTIENT and a.degree() >= b.degree():
        var fast = _ntt_divide_exact(_view(a.coefficients), _view(b.coefficients), budget)
        if fast:
            var exact = fast.take()
            if not exact:
                return None
            return DensePoly[ZZ]._from_canonical(a.variable, a.ring, exact.take())
    return _trial_divide(a, b, budget)


def _exact_integer(a: DensePoly[ZZ], b: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    var quotient = _divide_integer(a, b, budget)
    if not quotient:
        raise PolyError(PolyErrorKind.DOMAIN, "inexact integer squarefree quotient")
    return quotient.take()


def _gcd_cofactors(a: DensePoly[ZZ], b: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Tuple[DensePoly[ZZ], DensePoly[ZZ], DensePoly[ZZ]]:
    """(g, a / g, b / g) with g the gcd of _integer_gcd.

    The candidates are those of the modular primitive gcd; each is checked by
    dividing a and b themselves, so the accepted check leaves the cofactors.
    """
    if a.degree() <= 0 or b.degree() <= 0:
        var g = _integer_gcd(a, b, budget)
        return (g.copy(), _exact_integer(a, g, budget), _exact_integer(b, g, budget))
    var fd = _primitive(a, budget)
    var hd = _primitive(b, budget)
    var gamma = _gcd_coefficient(fd.coefficients[fd.degree()], hd.coefficients[hd.degree()], budget)
    var state = _GcdSearch(gamma^, min(fd.degree(), hd.degree()) + 1, fd, hd, budget)
    while True:
        budget.charge(1)
        var dense = _next_gcd_candidate(fd, hd, state, budget)
        if dense.degree() == 0:
            var one = _constant(a, ZZ.one(NoContext()), budget)
            return (one^, _copy(a, budget), _copy(b, budget))
        var g = _primitive(dense, budget)
        var left = _divide_integer(a, g, budget)
        if not left:
            continue
        var right = _divide_integer(b, g, budget)
        if right:
            return (g^, left.take(), right.take())


def _integer_squarefree(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> List[Tuple[DensePoly[ZZ], Int]]:
    """Primitive squarefree blocks of a primitive polynomial by Yun's recurrence over ZZ.

    By Gauss's lemma every quotient below is exact over ZZ, so modular gcds
    replace rational remainder sequences.
    """
    var blocks = List[Tuple[DensePoly[ZZ], Int]]()
    if f.degree() <= 0:
        return blocks^
    _allocation[Tuple[DensePoly[ZZ], Int]](1, budget)
    if _modular_squarefree(f, budget):
        blocks.append((_primitive(f, budget), 1))
        return blocks^
    var derivative = _derivative(f, budget)
    var first = _gcd_cofactors(f, derivative, budget)
    var b = first[1].copy()
    var c = first[2].copy()
    var multiplicity = 1
    while b.degree() > 0:
        budget.charge(1)
        # Exactly b = prod_{j>=i} a_j and c = b * sum_{j>=i} (j-i+1) a_j'/a_j;
        # b and c must share their scale, so neither is normalized.
        var d = c.sub(_derivative(b, budget), budget)
        var split = _gcd_cofactors(b, d, budget)
        if split[0].degree() > 0:
            _allocation[Tuple[DensePoly[ZZ], Int]](1, budget)
            blocks.append((_copy(split[0], budget), multiplicity))
        b = split[1].copy()
        c = split[2].copy()
        multiplicity = _sum(multiplicity, 1)
    return blocks^


def factor_integer(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Factorization[ZZ]:
    """Factor an integer polynomial into primitive factors with positive leading coefficients.

    Args:
        f: The integer polynomial, including zero or constants.
        budget: The shared work allowance.

    Returns:
        Signed content and factors sorted by degree, then leading-to-constant coefficients.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    var unit = _content(f, budget)
    var factors = List[Tuple[DensePoly[ZZ], Int]]()
    if f.degree() <= 0:
        return Factorization[ZZ](unit^, factors^)
    var primitive = _divide_scalar(f, unit, budget)
    var blocks = _integer_squarefree(primitive, budget)
    for entry in blocks:
        budget.charge(1)
        var irreducibles = _squarefree_factors(entry[0], budget)
        for i in range(len(irreducibles)):
            budget.charge(1)
            _allocation[Tuple[DensePoly[ZZ], Int]](1, budget)
            var irreducible = DensePoly[ZZ](f.variable, NoContext())
            swap(irreducible, irreducibles[i])
            factors.append((irreducible^, entry[1]))
    _sort(factors, f, budget)
    return Factorization[ZZ](unit^, factors^)


def factor(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> Factorization[QQ]:
    """Factor a rational polynomial into monic irreducibles with multiplicities.

    Args:
        f: The rational polynomial, including zero or constants.
        budget: The shared work allowance.

    Returns:
        The leading coefficient and factors in the deterministic coefficient order.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    var unit = QQ.zero(f.ring)
    if not f.is_zero():
        unit = _copy_coefficient(f.coefficients[f.degree()], budget)
    var factors = List[Tuple[DensePoly[QQ], Int]]()
    if f.degree() > 0:
        var integer_factors = factor_integer(_integer(f, budget), budget)
        for entry in integer_factors.factors:
            budget.charge(1)
            _allocation[Tuple[DensePoly[QQ], Int]](1, budget)
            factors.append((_monic(_rational(entry[0], budget), budget), entry[1]))
        _sort(factors, f, budget)
    return Factorization[QQ](unit^, factors^)
