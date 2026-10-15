"""Field-generic rational-function algorithms using public polynomial operations."""
from polynomials import ModInt, factor_finite_field
from polynomials import CoefficientField, DensePoly, SparsePoly, Budget, PolyError, PolyErrorKind, Factorization, QQ, QQi, factor_gaussian, gcd, extended_gcd, divmod, factor, evaluate, derivative, substitute_many, taylor_shift
from polynomials import CoefficientBuffer
from ._common import _allocation, _charge_size, _sum
from ._algebraic_value import Real, Complex
from . import _algebraic_factor
from std.memory import ArcPointer
from ._rational_certificates import _Factors, _retain_factors, _check_factors


@fieldwise_init
struct RationalFunction[C: CoefficientField](Copyable, Movable):
    var numerator: SparsePoly[Self.C]
    var denominator: SparsePoly[Self.C]

    @staticmethod
    def _reduce(p: SparsePoly[Self.C], q: SparsePoly[Self.C], mut common: SparsePoly[Self.C], mut budget: Budget) raises PolyError -> Self:
        _same(p, q, budget)
        if q.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "zero rational-function denominator")
        var result = gcd(p, q, budget)
        ref scalar = result.cofactors[1].terms[len(result.cofactors[1].terms) - 1].coefficient
        var inverse = scalar.inverse(budget)
        common = result.gcd.scale(scalar, budget)
        return Self(result.cofactors[0].scale(inverse, budget), result.cofactors[1].scale(inverse, budget))


@fieldwise_init
struct Cancellation[C: CoefficientField](Copyable, Movable):
    var common: SparsePoly[Self.C]
    var cofactors: Tuple[SparsePoly[Self.C], SparsePoly[Self.C]]


def _same[C: CoefficientField](a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError:
    a.check_invariants(budget)
    b.check_invariants(budget)
    if a.order != b.order or a.ring != b.ring:
        raise PolyError(PolyErrorKind.DOMAIN, "incompatible polynomial contexts or orders")


def _copy[C: CoefficientField](p: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    return DensePoly[C].from_coefficients(p.variable, p.ring, _slice(p.coefficients, 0, len(p.coefficients), budget), budget)


def _slice[C: CoefficientField](a: CoefficientBuffer[C], begin: Int, end: Int, mut budget: Budget) raises PolyError -> List[C]:
    _allocation[C](end - begin, budget)
    var out = List[C](capacity=end - begin)
    for i in range(begin, end):
        _charge_size(a[i].size(), budget)
        out.append(a[i].copy())
    return out^


def _constant[C: CoefficientField](p: DensePoly[C], value: C, mut budget: Budget) raises PolyError -> DensePoly[C]:
    _allocation[C](1, budget)
    _charge_size(value.size(), budget)
    return DensePoly[C].from_coefficients(p.variable, p.ring, List[C](length=1, fill=value), budget)


def _exact[C: CoefficientField](a: DensePoly[C], b: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    var result = divmod(a, b, budget)
    if not result.remainder.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "nonexact polynomial division")
    var out = DensePoly[C](a.variable, a.ring.copy())
    swap(out, result.quotient)
    return out^


def _remainder[C: CoefficientField](a: DensePoly[C], b: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    var result = divmod(a, b, budget)
    var out = DensePoly[C](a.variable, a.ring.copy())
    swap(out, result.remainder)
    return out^


def _factor[C: CoefficientField](a: DensePoly[C], mut budget: Budget) raises PolyError -> Factorization[C]:
    comptime if C == Real or C == Complex:
        return _algebraic_factor._factor(a, budget)
    elif C == ModInt:
        return rebind_var[Factorization[C]](factor_finite_field(rebind[DensePoly[ModInt]](a), budget))
    elif C == QQi:
        return rebind_var[Factorization[C]](factor_gaussian(rebind[DensePoly[QQi]](a), budget))
    else:
        comptime assert C == QQ
        return rebind_var[Factorization[C]](factor(rebind[DensePoly[QQ]](a), budget))


def _values[C: CoefficientField](a: List[C], length: Int, ctx: C.Context, mut budget: Budget) raises PolyError -> List[C]:
    _allocation[C](length, budget)
    var result = List[C](capacity=length)
    for i in range(length):
        if i < len(a):
            _charge_size(a[i].size(), budget)
            result.append(a[i].copy())
        else:
            result.append(C.zero(ctx))
    return result^


def _divide_values[C: CoefficientField](a: List[C], b: List[C], length: Int, ctx: C.Context, mut budget: Budget) raises PolyError -> List[C]:
    if length < 0 or not len(b) or b[0].is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "invalid formal division")
    var result = _values(a, length, ctx, budget)
    var inverse = b[0].inverse(budget)
    for i in range(length):
        for j in range(1, min(i + 1, len(b))):
            result[i] = result[i].sub(b[j].mul(result[i-j], budget), budget)
        result[i] = result[i].mul(inverse, budget)
    return result^


def _mul_values[C: CoefficientField](a: List[C], b: List[C], length: Int, ctx: C.Context, mut budget: Budget) raises PolyError -> List[C]:
    var result = _values(List[C](), length, ctx, budget)
    for i in range(min(length, len(a))):
        for j in range(min(length-i, len(b))):
            result[i+j] = result[i+j].add(a[i].mul(b[j], budget), budget)
    return result^


def _validate[C: CoefficientField](f: RationalFunction[C], mut budget: Budget) raises PolyError:
    f.numerator.check_invariants(budget)
    f.denominator.check_invariants(budget)
    _same(f.numerator, f.denominator, budget)
    if f.denominator.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "zero rational-function denominator")


def _normalized[C: CoefficientField](p: SparsePoly[C], q: SparsePoly[C], mut budget: Budget) raises PolyError -> RationalFunction[C]:
    var common = SparsePoly[C].zero(p.order, p.ring, budget)
    return RationalFunction[C]._reduce(p, q, common, budget)


def _canonical[C: CoefficientField](f: RationalFunction[C], mut budget: Budget) raises PolyError -> RationalFunction[C]:
    _validate(f, budget)
    return _normalized(f.numerator, f.denominator, budget)


def _univariate[C: CoefficientField](f: RationalFunction[C], variable: Int64, mut budget: Budget) raises PolyError -> Tuple[DensePoly[C], DensePoly[C]]:
    _validate(f, budget)
    return (DensePoly[C].from_sparse(f.numerator, budget, variable=variable), DensePoly[C].from_sparse(f.denominator, budget, variable=variable))


def rational_evaluate[C: CoefficientField](f: RationalFunction[C], assignments: List[Tuple[Int64, C]], mut budget: Budget) raises PolyError -> C:

    var reduced = _canonical(f, budget)
    var denominator = evaluate(reduced.denominator, assignments, budget)
    return evaluate(reduced.numerator, assignments, budget).mul(denominator.inverse(budget), budget)


def rational_derivative[C: CoefficientField](f: RationalFunction[C], variable: Int64, mut budget: Budget) raises PolyError -> RationalFunction[C]:

    _validate(f, budget)
    var numerator = derivative(f.numerator, variable, budget).mul(f.denominator, budget).sub(f.numerator.mul(derivative(f.denominator, variable, budget), budget), budget)
    var denominator = f.denominator.mul(f.denominator, budget)
    return _normalized(numerator, denominator, budget)


def rational_substitute[C: CoefficientField](f: RationalFunction[C], assignments: List[Tuple[Int64, SparsePoly[C]]], mut budget: Budget) raises PolyError -> RationalFunction[C]:

    var reduced = _canonical(f, budget)
    var numerator = substitute_many(reduced.numerator, assignments, budget)
    var denominator = substitute_many(reduced.denominator, assignments, budget)
    return _normalized(numerator, denominator, budget)


def _homogenize[C: CoefficientField](poly: DensePoly[C], numerator: DensePoly[C], denominator: DensePoly[C], degree: Int, mut budget: Budget) raises PolyError -> DensePoly[C]:
    var result = DensePoly[C](poly.variable, poly.ring.copy())
    var power = _constant(poly, C.one(poly.ring), budget)
    for k in range(degree, -1, -1):
        result = result.mul(numerator, budget)
        if k < len(poly.coefficients):
            result = result.add(power.scale(poly.coefficients[k], budget), budget)
        if k:
            power = power.mul(denominator, budget)
    return result^


def rational_compose[C: CoefficientField](outer: RationalFunction[C], inner: RationalFunction[C], variable: Int64, mut budget: Budget) raises PolyError -> RationalFunction[C]:

    var f = _canonical(outer, budget)
    var g = _canonical(inner, budget)
    _same(f.numerator, g.numerator, budget)
    var a = _univariate(f, variable, budget)
    var b = _univariate(g, variable, budget)
    var degree = max(a[0].degree(), a[1].degree())
    var numerator = _homogenize(a[0], b[0], b[1], degree, budget)
    var denominator = _homogenize(a[1], b[0], b[1], degree, budget)
    return _normalized(numerator.to_sparse(f.numerator.order, budget), denominator.to_sparse(f.numerator.order, budget), budget)


@fieldwise_init
struct PartialFractionTerm[C: CoefficientField](Copyable, Movable):

    var numerator: DensePoly[Self.C]

    var factor: DensePoly[Self.C]

    var exponent: Int



@fieldwise_init
struct PartialFractions[C: CoefficientField](Copyable, Movable):

    var polynomial: DensePoly[Self.C]

    var terms: List[PartialFractionTerm[Self.C]]

    var factors: ArcPointer[_Factors[Self.C]]



def partial_fractions[C: CoefficientField](f: RationalFunction[C], variable: Int64, mut budget: Budget) raises PolyError -> PartialFractions[C]:

    var reduced = _canonical(f, budget)
    var pair = _univariate(reduced, variable, budget)
    var division = divmod(pair[0], pair[1], budget)
    var factors = _factor(pair[1], budget)
    var terms = List[PartialFractionTerm[C]]()
    for entry in factors.factors:
        var power = entry[0].pow(UInt64(entry[1]), budget)
        var rest = _exact(pair[1], power, budget)
        var bezout = extended_gcd(rest, power, budget)
        var inverse = DensePoly[C].from_sparse(bezout.bezout.value()[0], budget, variable=variable)
        var local = _remainder(division.remainder.mul(inverse, budget), power, budget)
        for exponent in range(entry[1], 0, -1):
            var digit = divmod(local, entry[0], budget)
            if not digit.remainder.is_zero():
                _allocation[PartialFractionTerm[C]](1, budget)
                terms.append(PartialFractionTerm[C](_copy(digit.remainder, budget), _copy(entry[0], budget), exponent))
            swap(local, digit.quotient)
    var polynomial = DensePoly[C](variable, pair[0].ring.copy())
    swap(polynomial, division.quotient)
    return PartialFractions[C](polynomial^, terms^, _retain_factors(factors^, budget))


def check_partial_fractions[C: CoefficientField](f: RationalFunction[C], variable: Int64, decomposition: PartialFractions[C], mut budget: Budget) raises PolyError -> Bool:

    var pair = _univariate(f, variable, budget)
    if not _check_factors(pair[1], decomposition.factors[], budget):
        return False
    if len(decomposition.terms) > max(0, pair[1].degree()):
        return False
    decomposition.polynomial.check_invariants(budget)
    if decomposition.polynomial.variable != variable:
        return False
    var total = decomposition.polynomial.mul(pair[1], budget)
    for i in range(len(decomposition.terms)):
        budget.charge(1)
        ref term = decomposition.terms[i]
        if term.numerator.is_zero():
            return False
        var found = False
        for entry in decomposition.factors[].value.factors:
            budget.charge(1)
            if term.factor.sub(entry[0], budget).is_zero():
                found = term.exponent > 0 and term.exponent <= entry[1]
                break
        if not found:
            return False
        for j in range(i):
            budget.charge(1)
            if term.exponent == decomposition.terms[j].exponent and term.factor.sub(decomposition.terms[j].factor, budget).is_zero():
                return False
        term.factor.check_invariants(budget)
        term.numerator.check_invariants(budget)
        if term.exponent <= 0 or term.factor.variable != variable or term.numerator.variable != variable or term.factor.degree() < 1 or term.numerator.degree() >= term.factor.degree():
            return False
        _charge_size(term.factor.coefficients[term.factor.degree()].size(), budget)
        if term.factor.coefficients[term.factor.degree()] != C.one(pair[0].ring):
            return False
        if term.exponent > pair[1].degree() // term.factor.degree():
            return False
        var power = term.factor.pow(UInt64(term.exponent), budget)
        var division = divmod(pair[1], power, budget)
        if not division.remainder.is_zero():
            return False
        total = total.add(term.numerator.mul(division.quotient, budget), budget)
    for coefficient_index in range(len(total.coefficients)):
        ref coefficient = total.coefficients[coefficient_index]
        _charge_size(coefficient.size(), budget)
    for coefficient_index in range(len(pair[0].coefficients)):
        ref coefficient = pair[0].coefficients[coefficient_index]
        _charge_size(coefficient.size(), budget)
    return total == pair[0]


def rational_poles[C: CoefficientField](f: RationalFunction[C], variable: Int64, mut budget: Budget) raises PolyError -> Factorization[C]:

    var reduced = _canonical(f, budget)
    var pair = _univariate(reduced, variable, budget)
    return _factor(pair[1], budget)


@fieldwise_init
struct LaurentExpansion[C: CoefficientField](Copyable, Movable):

    var variable: Int64

    var center: Self.C

    var valuation: Int

    var coefficients: List[Self.C]



def _local[C: CoefficientField](f: RationalFunction[C], variable: Int64, center: C, mut budget: Budget) raises PolyError -> Tuple[List[C], List[C], Int]:
    var pair = _univariate(f, variable, budget)
    var numerator = taylor_shift(pair[0], center, budget)
    var denominator = taylor_shift(pair[1], center, budget)
    var ni = 0
    var di = 0
    while ni < len(numerator.coefficients):
        _charge_size(numerator.coefficients[ni].size(), budget)
        if not numerator.coefficients[ni].is_zero():
            break
        ni += 1
    while di < len(denominator.coefficients):
        _charge_size(denominator.coefficients[di].size(), budget)
        if not denominator.coefficients[di].is_zero():
            break
        di += 1
    var valuation = ni - di if not numerator.is_zero() else 0
    return (_slice(numerator.coefficients, ni, len(numerator.coefficients), budget), _slice(denominator.coefficients, di, len(denominator.coefficients), budget), valuation)


def laurent_expansion[C: CoefficientField](f: RationalFunction[C], variable: Int64, center: C, terms: Int, mut budget: Budget) raises PolyError -> LaurentExpansion[C]:

    var local = _local(f, variable, center, budget)
    var coefficients = _divide_values[C](local[0], local[1], terms, center.context(), budget)
    _charge_size(center.size(), budget)
    return LaurentExpansion[C](variable, center.copy(), local[2], coefficients^)


def check_laurent_expansion[C: CoefficientField](f: RationalFunction[C], expansion: LaurentExpansion[C], mut budget: Budget) raises PolyError -> Bool:

    var local = _local(f, expansion.variable, expansion.center, budget)
    if local[2] != expansion.valuation:
        return False
    var length = len(expansion.coefficients)
    var product = _mul_values[C](local[1], expansion.coefficients, length, expansion.center.context(), budget)
    var actual = _values[C](product, length, expansion.center.context(), budget)
    var expected = _values[C](local[0], length, expansion.center.context(), budget)
    for i in range(length):
        _charge_size(_sum(actual[i].size(), expected[i].size()), budget)
        if actual[i] != expected[i]:
            return False
    return True


def rational_residue[C: CoefficientField](f: RationalFunction[C], variable: Int64, center: C, mut budget: Budget) raises PolyError -> C:

    var local = _local(f, variable, center, budget)
    if local[2] >= 0:
        return C.zero(center.context())
    var coefficients = _divide_values[C](local[0], local[1], -local[2], center.context(), budget)
    return coefficients.pop()
