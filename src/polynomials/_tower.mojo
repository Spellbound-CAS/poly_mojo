"""Internal views of sparse polynomials as univariates over smaller rings."""
from .errors import Budget, PolyError, PolyErrorKind
from .coefficients import CoefficientRing
from .sparse import SparsePoly
from .division import divide_exact


from ._views import _constant, _one, _zero, _split, _shift, _lc, _variables


def _exact_poly[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var quotient = divide_exact(a, b, budget)
    if not quotient:
        raise PolyError(PolyErrorKind.DOMAIN, "internal polynomial division is not exact")
    return quotient.take()


def _pseudo[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, carry: Bool, mut budget: Budget) raises PolyError -> Tuple[SparsePoly[C], SparsePoly[C], SparsePoly[C]]:
    # Return q, r, multiplier with multiplier*f = q*g+r.
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "zero pseudo-divisor")
    var n = Int(g.degree_in(variable, budget))
    var m = Int(f.degree_in(variable, budget))
    var delta = 0 if f.is_zero() else max(m - n + 1, 0)
    var remaining = delta
    var leading = _lc(g, variable, budget)
    var q = _zero(f, budget)
    var r = f._copy(budget)
    while not r.is_zero():
        budget.charge(1)
        var degree = Int(r.degree_in(variable, budget))
        if degree < n:
            break
        var factor = _shift(_lc(r, variable, budget), variable, UInt32(degree - n), budget)
        if carry:
            q = q.mul(leading, budget).add(factor, budget)
        r = r.mul(leading, budget).sub(factor.mul(g, budget), budget)
        remaining -= 1
    var suffix = leading.pow(UInt64(remaining), budget)
    var multiplier = _one(f, budget)
    if carry:
        q = q.mul(suffix, budget)
        multiplier = leading.pow(UInt64(delta), budget)
    return (q^, r.mul(suffix, budget), multiplier^)
