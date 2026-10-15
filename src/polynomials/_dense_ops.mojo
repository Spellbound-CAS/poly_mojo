"""Shared budgeted dense views and field operations."""
from .coefficients import CoefficientRing, CoefficientField
from .dense import DensePoly, _slice
from ._coefficient_ops import _copy as _copy_coefficient
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .evaluate import _derivative_dense as _derivative
from .division import _dense_divide
from ._euclidean import _extended


def _copy[C: CoefficientRing](f: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    return DensePoly[C]._from_canonical(f.variable, f.ring, _slice(_view(f.coefficients), 0, len(f.coefficients), budget))


def _constant[C: CoefficientRing](f: DensePoly[C], c: C, mut budget: Budget) raises PolyError -> DensePoly[C]:
    _allocation[C](1, budget)
    return DensePoly[C].from_coefficients(f.variable, f.ring, [_copy_coefficient(c, budget)], budget)


def _one[C: CoefficientRing](f: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    return _constant(f, C.one(f.ring), budget)


def _monic[C: CoefficientField](f: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    if f.is_zero():
        return _copy(f, budget)
    return f.scale(f.coefficients[f.degree()].inverse(budget), budget)


def _gcd[C: CoefficientField](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    var result = _extended(f, g, budget, carry=False)
    var d = DensePoly[C](f.variable, f.ring.copy())
    swap(d, result[0])
    return d^


def _remainder[C: CoefficientField](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    var result = _dense_divide(f, g, False, budget, quotient=False)
    var r = DensePoly[C](f.variable, f.ring.copy())
    swap(r, result.remainder)
    return r^


def _exact[C: CoefficientField](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    var result = _dense_divide(f, g, False, budget)
    if not result.remainder.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "inexact dense field quotient")
    var q = DensePoly[C](f.variable, f.ring.copy())
    swap(q, result.quotient)
    return q^
