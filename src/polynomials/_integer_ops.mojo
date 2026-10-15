"""Dense integer content, scalar division and primitive normalization."""
from apn_mojo import Integer
from .coefficients import ZZ, _work
from .dense import DensePoly
from .errors import Budget, PolyError, _allocation
from ._coefficient_ops import _exact as _coefficient_exact, _gcd_coefficient
from ._dense_ops import _copy


def _divide_scalar(f: DensePoly[ZZ], scalar: ZZ, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(f.coefficients), budget)
    var values = List[ZZ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        values.append(_coefficient_exact(c, scalar, budget))
    return DensePoly[ZZ]._from_canonical(f.variable, f.ring, values^)


def _content(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> ZZ:
    var result = ZZ.zero(f.ring)
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        result = _gcd_coefficient(result, c, budget)
        _work(result.size(), 1, False, budget)
        if result.value == Integer(1):
            break
    if not f.is_zero() and f.coefficients[f.degree()].value.sign() < 0:
        result = result.neg(budget)
    return result^


def _primitive(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    if f.is_zero():
        return _copy(f, budget)
    return _divide_scalar(f, _content(f, budget), budget)
