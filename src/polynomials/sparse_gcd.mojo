"""Certified sparse modular gcds over the integers and rationals."""
from .errors import Budget, PolyError, PolyErrorKind
from .coefficients import CoefficientRing, ZZ, QQ
from .sparse import SparsePoly
from .monomial import MonomialOrder
from .certificates import GcdResult
from ._sparse_gcd import _sparse_gcd
from ._gcd import _normalized
from ._images import _integer_image, _rational_image
from ._views import _variables, _zero
from ._tower import _exact_poly
from .gcd import coprimality_certificate


def _certified_result[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], var common: SparsePoly[C], mut budget: Budget) raises PolyError -> GcdResult[C]:
    var left = _zero(f, budget)
    var right = _zero(g, budget)
    if common.is_zero():
        return GcdResult[C](common^, (left^, right^), None, None)
    left = _exact_poly(f, common, budget)
    right = _exact_poly(g, common, budget)
    var proof = coprimality_certificate(left, right, budget)
    if not proof:
        raise PolyError(PolyErrorKind.DOMAIN, "reconstructed gcd lacks maximality evidence")
    return GcdResult[C](common^, (left^, right^), None, proof^)


def sparse_gcd_integer(f: SparsePoly[ZZ], g: SparsePoly[ZZ], mut budget: Budget) raises PolyError -> GcdResult[ZZ]:
    """Compute an integer gcd using modular support discovery and sparse interpolation.

    Args:
        f: The first sparse integer polynomial.
        g: The second polynomial in the same stored order.
        budget: The shared work allowance.

    Returns:
        The normalized gcd, cofactors, and explicit coprimality evidence when nonzero.

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    var order = MonomialOrder(MonomialOrder.LEX, _variables(f, g, budget), budget)
    var common = _normalized(_sparse_gcd(f.reorder(order, budget), g.reorder(order, budget), budget).reorder(f.order, budget), budget)
    return _certified_result(f, g, common^, budget)


def sparse_gcd(f: SparsePoly[QQ], g: SparsePoly[QQ], mut budget: Budget) raises PolyError -> GcdResult[QQ]:
    """Compute a monic rational gcd by sparse modular reconstruction.

    Args:
        f: The first sparse rational polynomial.
        g: The second polynomial in the same stored order.
        budget: The shared work allowance.

    Returns:
        The monic gcd with cofactors and independently checkable maximality evidence.

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    var order = MonomialOrder(MonomialOrder.LEX, _variables(f, g, budget), budget)
    var common = _normalized(_rational_image(_sparse_gcd(_integer_image(f.reorder(order, budget), budget), _integer_image(g.reorder(order, budget), budget), budget), budget).reorder(f.order, budget), budget)
    return _certified_result(f, g, common^, budget)
