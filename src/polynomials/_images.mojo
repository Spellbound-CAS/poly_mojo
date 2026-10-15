"""Coefficient-ring images independent of gcd and factorization algorithms."""
from apn_mojo import Integer, Rational, lcm as integer_lcm
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from .coefficients import ZZ, QQ, ModInt, Modulus, NoContext, _work
from .sparse import SparsePoly, Term
from .dense import DensePoly
from .primes import _residue


def _denominator(common: Integer, value: QQ, mut budget: Budget) raises PolyError -> Integer:
    _work(common.magnitude_bit_length() // 8 + 1, value.size(), True, budget)
    try:
        return integer_lcm(common, value.value.denominator())
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _scaled_integer(value: QQ, common: Integer, mut budget: Budget) raises PolyError -> ZZ:
    _work(value.size(), common.magnitude_bit_length() // 8 + 1, True, budget)
    try:
        return ZZ((value.value * Rational(common)).to_integer_exact())
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _integer_dense(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    var common = Integer(1)
    for coefficient_index in range(len(f.coefficients)):
        ref value = f.coefficients[coefficient_index]
        budget.charge(1)
        if not value.is_zero():
            common = _denominator(common, value, budget)
    _allocation[ZZ](len(f.coefficients), budget)
    var values = List[ZZ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref value = f.coefficients[coefficient_index]
        budget.charge(1)
        values.append(ZZ.zero(NoContext()) if value.is_zero() else _scaled_integer(value, common, budget))
    return DensePoly[ZZ]._from_canonical(f.variable, NoContext(), values^)


def _rational_dense(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    _allocation[QQ](len(f.coefficients), budget)
    var values = List[QQ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref value = f.coefficients[coefficient_index]
        budget.charge(1)
        _charge_size(value.size(), budget)
        values.append(QQ(Rational(value.value)))
    return DensePoly[QQ]._from_canonical(f.variable, NoContext(), values^)


def _integer_image(p: SparsePoly[QQ], mut budget: Budget) raises PolyError -> SparsePoly[ZZ]:
    var common = Integer(1)
    for term in p.terms:
        budget.charge(1)
        common = _denominator(common, term.coefficient, budget)
    _allocation[Term[ZZ]](len(p.terms), budget)
    var terms = List[Term[ZZ]](capacity=len(p.terms))
    for term in p.terms:
        budget.charge(1)
        _charge_size(len(term.monomial.pairs), budget)
        terms.append(Term[ZZ](term.monomial.copy(), _scaled_integer(term.coefficient, common, budget)))
    var result = SparsePoly[ZZ].zero(p.order, NoContext(), budget)
    result.terms = terms^
    return result^


def _rational_image(p: SparsePoly[ZZ], mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    _allocation[Term[QQ]](len(p.terms), budget)
    var terms = List[Term[QQ]](capacity=len(p.terms))
    for term in p.terms:
        budget.charge(1)
        _charge_size(_sum(len(term.monomial.pairs), term.coefficient.size()), budget)
        terms.append(Term[QQ](term.monomial.copy(), QQ(Rational(term.coefficient.value))))
    var result = SparsePoly[QQ].zero(p.order, NoContext(), budget)
    result.terms = terms^
    return result^


def _mod_image(p: DensePoly[ZZ], ctx: Modulus, mut budget: Budget) raises PolyError -> DensePoly[ModInt]:
    _allocation[ModInt](len(p.coefficients), budget)
    var coefficients = List[ModInt](capacity=len(p.coefficients))
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        budget.charge(1)
        coefficients.append(ModInt(_residue(c.value, ctx.value, budget), ctx))
    return DensePoly[ModInt].from_coefficients(p.variable, ctx, coefficients^, budget)


def _sparse_mod_image(f: SparsePoly[ZZ], ctx: Modulus, mut budget: Budget) raises PolyError -> SparsePoly[ModInt]:
    _allocation[Term[ModInt]](len(f.terms), budget)
    var terms = List[Term[ModInt]]()
    for term in f.terms:
        budget.charge(1)
        _charge_size(len(term.monomial.pairs), budget)
        var value = ModInt(_residue(term.coefficient.value, ctx.value, budget), ctx)
        if not value.is_zero():
            terms.append(Term[ModInt](term.monomial.copy(), value^))
    var result = SparsePoly[ModInt].zero(f.order, ctx, budget)
    result.terms = terms^
    return result^
