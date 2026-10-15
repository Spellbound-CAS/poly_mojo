"""Independent irreducibility verification for rational-function denominators.

Prime reductions handle the common case. Kronecker's integer interpolation
argument supplies a finite fallback, including irreducibles with no irreducible
prime reduction. Gaussian factors use a squarefree shifted conjugate norm.
The checker never factors the source denominator again.
"""
from std.memory import ArcPointer
from apn_mojo import Integer, Rational, ExactComplex, div_rem_euclid, lcm
from polynomials import Budget, PolyError, PolyErrorKind, CoefficientField, QQ, QQi, ZZ, ModInt, Modulus, NoContext, DensePoly, Factorization, check_factorization, divmod, evaluate_dense, MultipointPlan, taylor_shift, derivative_dense
from ._common import _allocation, _charge_size, _sum, _product, _grow
from ._algebraic_value import Real, Complex
from . import _algebraic_factor


@fieldwise_init
struct _Irreducibility(Copyable, Movable):
    # prime == 0 selects exhaustive integer interpolation. A linear factor
    # needs neither a prime nor a shift. Shifts apply only to QQi norms.
    var prime: UInt64
    var shift: Int


@fieldwise_init
struct _Factors[C: CoefficientField](Copyable, Movable):
    var value: Factorization[Self.C]
    var certificates: List[_Irreducibility]


def _copy[C: CoefficientField](p: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    _allocation[C](len(p.coefficients), budget)
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        _charge_size(c.size(), budget)
    return p.copy()


def _remainder[C: CoefficientField](p: DensePoly[C], q: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    var division = divmod(p, q, budget)
    var result = DensePoly[C](p.variable, p.ring.copy())
    swap(result, division.remainder)
    return result^


def _coprime[C: CoefficientField](p: DensePoly[C], q: DensePoly[C], mut budget: Budget) raises PolyError -> Bool:
    var a = _copy(p, budget)
    var b = _copy(q, budget)
    while not b.is_zero():
        budget.charge(1)
        var remainder = _remainder(a, b, budget)
        a = b^
        b = remainder^
    return a.degree() == 0


def _work(a: Integer, b: Integer, mut budget: Budget) raises PolyError:
    _charge_size(_product(_sum(Int(a.magnitude_bit_length() // 8), 1), _sum(Int(b.magnitude_bit_length() // 8), 1)), budget)


def _integer(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # Keep integer-valued coefficients in QQ so interpolation and division use
    # the same public field operations. Gauss's lemma needs no primitive form.
    var common = Integer(1)
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        _work(common, c.value.denominator(), budget)
        try:
            common = lcm(common, c.value.denominator())
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    _charge_size(Int(common.magnitude_bit_length() // 8) + 1, budget)
    return p.scale(QQ(Rational(common)), budget)


def _residue(value: Integer, ctx: Modulus, mut budget: Budget) raises PolyError -> ModInt:
    var modulus = Integer(ctx.value)
    _work(value, modulus, budget)
    try:
        var pair = div_rem_euclid(value, modulus)
        return ModInt(UInt64(Int(pair[1])), ctx)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _power(p: DensePoly[ModInt], exponent: UInt64, modulus: DensePoly[ModInt], mut budget: Budget) raises PolyError -> DensePoly[ModInt]:
    _allocation[ModInt](1, budget)
    var out = DensePoly[ModInt].from_coefficients(p.variable, p.ring, [ModInt.one(p.ring)], budget)
    var base = _copy(p, budget)
    var n = exponent
    while n:
        budget.charge(1)
        if n & 1:
            out = _remainder(out.mul(base, budget), modulus, budget)
        n >>= 1
        if n:
            base = _remainder(base.mul(base, budget), modulus, budget)
    return out^


def _modular(p: DensePoly[QQ], prime: UInt64, mut budget: Budget) raises PolyError -> Bool:
    # A degree-n reducible polynomial has an irreducible factor of degree at
    # most floor(n/2). Such a degree-d factor divides x^(prime^d)-x.
    var ctx = Modulus(prime, budget)
    _allocation[ModInt](len(p.coefficients), budget)
    var coefficients = List[ModInt](capacity=len(p.coefficients))
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        var numerator = _residue(c.value.numerator(), ctx, budget)
        var denominator = _residue(c.value.denominator(), ctx, budget)
        if denominator.is_zero():
            return False
        coefficients.append(numerator.mul(denominator.inverse(budget), budget))
    var image = DensePoly[ModInt].from_coefficients(p.variable, ctx, coefficients^, budget)
    if image.degree() != p.degree():
        return False
    _allocation[ModInt](2, budget)
    var x = DensePoly[ModInt].from_coefficients(p.variable, ctx, [ModInt.zero(ctx), ModInt.one(ctx)], budget)
    var power = _copy(x, budget)
    for _ in range(p.degree() // 2):
        budget.charge(1)
        power = _power(power, prime, image, budget)
        if not _coprime(image, power.sub(x, budget), budget):
            return False
    return True


def _divisors(value: QQ, positive_only: Bool, mut budget: Budget) raises PolyError -> List[QQ]:
    var n = ZZ(value.value.numerator())
    if n.value.sign() < 0:
        n = n.neg(budget)
    var result = List[QQ]()
    var d = ZZ.one(NoContext())
    while True:
        budget.charge(1)
        _work(n.value, d.value, budget)
        var quotient = Integer(0)
        var remainder = Integer(0)
        try:
            var pair = div_rem_euclid(n.value, d.value)
            swap(quotient, pair[0])
            swap(remainder, pair[1])
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        _work(d.value, quotient, budget)
        if d.value > quotient:
            break
        if remainder == 0:
            _grow(result, budget)
            _charge_size(d.size(), budget)
            result.append(QQ(Rational(d.value)))
            if not positive_only:
                _grow(result, budget)
                result.append(QQ(Rational(d.neg(budget).value)))
            if d.value != quotient:
                _grow(result, budget)
                _charge_size(Int(quotient.magnitude_bit_length() // 8) + 1, budget)
                result.append(QQ(Rational(quotient)))
                if not positive_only:
                    _grow(result, budget)
                    _charge_size(Int(quotient.magnitude_bit_length() // 8) + 1, budget)
                    result.append(QQ(Rational(-quotient)))
        d = d.add(ZZ.one(NoContext()), budget)
    return result^


def _kronecker(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> Bool:
    # Every rational proper factor has, up to sign, an integer factor of degree
    # <= n/2. Its values at h+1 distinct integers divide the source values.
    # Fix a positive first value, enumerate the Cartesian product iteratively,
    # and test every interpolant. Store only divisor lists and one candidate.
    var integer = _integer(p, budget)
    var count = _sum(p.degree() // 2, 1)
    _allocation[QQ](count, budget)
    _allocation[List[QQ]](count, budget)
    var points = List[QQ](capacity=count)
    var divisors = List[List[QQ]](capacity=count)
    for i in range(count):
        budget.charge(1)
        var x = (i + 1) // 2 if i % 2 else -(i // 2)
        var point = QQ(Rational(x))
        var value = evaluate_dense(integer, point, budget)
        if value.is_zero():
            return False
        divisors.append(_divisors(value, i == 0, budget))
        points.append(point^)
    var plan = MultipointPlan[QQ].prepare_interpolation(p.variable, NoContext(), points, budget)
    _allocation[Int](count, budget)
    var indices = List[Int](length=count, fill=0)
    while True:
        budget.charge(1)
        _allocation[QQ](count, budget)
        var values = List[QQ](capacity=count)
        for i in range(count):
            _charge_size(divisors[i][indices[i]].size(), budget)
            values.append(divisors[i][indices[i]].copy())
        var candidate = plan.interpolate(values, budget)
        var integral = candidate.degree() > 0
        for coefficient_index in range(len(candidate.coefficients)):
            ref c = candidate.coefficients[coefficient_index]
            _charge_size(c.size(), budget)
            integral = integral and c.value.denominator() == 1
        if integral and _remainder(integer, candidate, budget).is_zero():
            return False
        var i = count - 1
        while i >= 0:
            budget.charge(1)
            indices[i] += 1
            if indices[i] < len(divisors[i]):
                break
            indices[i] = 0
            i -= 1
        if i < 0:
            return True


def _rational_certificate(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> _Irreducibility:
    # Finite prime search is a choice of verification route, never evidence of
    # irreducibility by itself. Some irreducibles split modulo every prime.
    _allocation[UInt64](8, budget)
    var primes: List[UInt64] = [2, 3, 5, 7, 11, 13, 17, 19]
    for prime in primes:
        if _modular(p, prime, budget):
            return _Irreducibility(prime, 0)
    return _Irreducibility(0, 0)


def _norm(p: DensePoly[QQi], shift: Int, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    var translated = taylor_shift(p, QQi(ExactComplex(0, shift)), budget)
    _allocation[QQ](len(translated.coefficients), budget)
    _allocation[QQ](len(translated.coefficients), budget)
    var real = List[QQ](capacity=len(translated.coefficients))
    var imag = List[QQ](capacity=len(translated.coefficients))
    for coefficient_index in range(len(translated.coefficients)):
        ref c = translated.coefficients[coefficient_index]
        _charge_size(c.size(), budget)
        real.append(QQ(c.value.real()))
        imag.append(QQ(c.value.imag()))
    var a = DensePoly[QQ].from_coefficients(p.variable, NoContext(), real^, budget)
    var b = DensePoly[QQ].from_coefficients(p.variable, NoContext(), imag^, budget)
    return a.mul(a, budget).add(b.mul(b, budget), budget)


def _certificate[C: CoefficientField](p: DensePoly[C], mut budget: Budget) raises PolyError -> _Irreducibility:
    if p.degree() == 1:
        comptime if C == ModInt:
            return _Irreducibility(rebind[DensePoly[ModInt]](p).ring.value, 0)
        else:
            return _Irreducibility(0, 0)
    comptime if C == ModInt:
        return _Irreducibility(rebind[DensePoly[ModInt]](p).ring.value, 0)
    elif C == QQ:
        return _rational_certificate(rebind[DensePoly[QQ]](p), budget)
    else:
        comptime assert C == QQi
        var shift = 0
        while True:
            budget.charge(1)
            var norm = _norm(rebind[DensePoly[QQi]](p), shift, budget)
            if _coprime(norm, derivative_dense(norm, budget), budget):
                var certificate = _rational_certificate(norm, budget)
                certificate.shift = shift
                return certificate^
            shift = _sum(shift, 1)


def _check_irreducible[C: CoefficientField](p: DensePoly[C], proof: _Irreducibility, mut budget: Budget) raises PolyError -> Bool:
    p.check_invariants(budget)
    if p.degree() < 1 or proof.shift < 0:
        return False
    comptime if C == ModInt:
        return proof.shift == 0 and proof.prime == rebind[DensePoly[ModInt]](p).ring.value and _finite_irreducible(rebind[DensePoly[ModInt]](p), budget)
    else:
        if p.degree() == 1:
            return proof.prime == 0 and proof.shift == 0
        var rational: DensePoly[QQ]
        comptime if C == QQ:
            if proof.shift != 0:
                return False
            rational = _copy(rebind[DensePoly[QQ]](p), budget)
        else:
            comptime assert C == QQi
            rational = _norm(rebind[DensePoly[QQi]](p), proof.shift, budget)
        if proof.prime == 0:
            return _kronecker(rational, budget)
        try:
            return _modular(rational, proof.prime, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                return False
            raise e^

def _retain_factors[C: CoefficientField](var factors: Factorization[C], mut budget: Budget) raises PolyError -> ArcPointer[_Factors[C]]:
    var certificates = List[_Irreducibility]()
    comptime if C == QQ or C == QQi or C == ModInt:
        _allocation[_Irreducibility](len(factors.factors), budget)
        certificates.reserve(len(factors.factors))
        for entry in factors.factors:
            certificates.append(_certificate(entry[0], budget))
    _allocation[_Factors[C]](1, budget)
    return ArcPointer(_Factors[C](factors^, certificates^))


def _check_factors[C: CoefficientField](p: DensePoly[C], proof: _Factors[C], mut budget: Budget) raises PolyError -> Bool:
    comptime if C == Real or C == Complex:
        return len(proof.certificates) == 0 and _algebraic_factor._check(p, proof.value, budget)
    else:
        if len(proof.certificates) != len(proof.value.factors) or len(proof.value.factors) > max(0, p.degree()):
            return False
        var total = 0
        for i in range(len(proof.value.factors)):
            budget.charge(1)
            ref factor = proof.value.factors[i][0]
            var multiplicity = proof.value.factors[i][1]
            factor.check_invariants(budget)
            var d = factor.degree()
            if factor.variable != p.variable or d < 1 or multiplicity <= 0 or multiplicity > (p.degree() - total) // d:
                return False
            total += d * multiplicity
            if not factor.coefficients[d].sub(C.one(p.ring), budget).is_zero():
                return False
            for j in range(i):
                if factor.sub(proof.value.factors[j][0], budget).is_zero():
                    return False
            if not _check_irreducible(factor, proof.certificates[i], budget):
                return False
        return total == p.degree() and check_factorization(p, proof.value, budget)


def _finite_irreducible(p: DensePoly[ModInt], mut budget: Budget) raises PolyError -> Bool:
    if p.degree() < 1:
        return False
    var ctx = p.ring
    _allocation[ModInt](2, budget)
    var x = DensePoly[ModInt].from_coefficients(p.variable, ctx, [ModInt.zero(ctx), ModInt.one(ctx)], budget)
    var power = _remainder(x, p, budget)
    for _ in range(p.degree() // 2):
        power = _power(power, ctx.value, p, budget)
        if not _coprime(p, power.sub(x, budget), budget):
            return False
    return True
