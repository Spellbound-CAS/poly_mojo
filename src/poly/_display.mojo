"""Named-variable display through the public, budgeted core formatter."""
from polynomials import PolynomialRing, CoefficientRing, SparsePoly, Budget, PolyError, PolyErrorKind
from ._common import _allocation, _charge_size


def _format[C: CoefficientRing](source: SparsePoly[C], width: Int, requested: List[String], mut budget: Budget) raises PolyError -> String:
    if len(requested) and len(requested) != width:
        raise PolyError(PolyErrorKind.DOMAIN, "provide one display name per array axis")
    _allocation[String](width, budget)
    var names = List[String](capacity=width)
    for i in range(width):
        if len(requested):
            _charge_size(requested[i].byte_length(), budget)
            names.append(requested[i].copy())
        else:
            budget.charge(32)
            names.append("x" if width == 1 else String("x", i))
    var ring = PolynomialRing[C].create(names, source.ring, budget)
    return ring.format(source, budget)
