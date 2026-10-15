"""Shared field Euclidean loop, with optional Bézout propagation."""
from .coefficients import CoefficientField
from .errors import Budget, PolyError
from .dense import DensePoly
from .division import _dense_divide


def _extended[C: CoefficientField](left: DensePoly[C], right: DensePoly[C], mut budget: Budget, carry: Bool = True) raises PolyError -> Tuple[DensePoly[C], DensePoly[C], DensePoly[C]]:
    left._same(right)
    var f = left.scale(C.one(left.ring), budget)
    var h = right.scale(C.one(right.ring), budget)
    var u0 = DensePoly[C].from_coefficients(left.variable, left.ring, [C.one(left.ring)], budget)
    var v0 = DensePoly[C](left.variable, left.ring.copy())
    var u1 = DensePoly[C](left.variable, left.ring.copy())
    var v1 = DensePoly[C].from_coefficients(left.variable, left.ring, [C.one(left.ring)], budget)
    while not h.is_zero():
        budget.charge(1)
        var division = _dense_divide(f, h, False, budget, quotient=carry)
        var u = DensePoly[C](left.variable, left.ring.copy())
        var v = DensePoly[C](left.variable, left.ring.copy())
        if carry:
            u = u0.sub(division.quotient.mul(u1, budget), budget)
            v = v0.sub(division.quotient.mul(v1, budget), budget)
        swap(f, h)
        swap(h, division.remainder)
        u0 = u1^
        v0 = v1^
        u1 = u^
        v1 = v^
    if not f.is_zero():
        var scale = f.coefficients[len(f.coefficients) - 1].inverse(budget)
        f = f.scale(scale, budget)
        if carry:
            u0 = u0.scale(scale, budget)
            v0 = v0.scale(scale, budget)
    return (f^, u0^, v0^)
