"""Fixed section prefixes shared by CAD producers and independent checkers."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .algebraic_types import RealAlgebraic, _check_real, _copy_real, _same_root, _check_comparison
from .algebraic import algebraic_compare
from .cad_types import CADStackCertificate
from ._cad_lifting import _position


def _validate_fixed(fixed: List[RealAlgebraic], width: Int, mut budget: Budget) raises PolyError:
    if len(fixed) > width:
        raise PolyError(PolyErrorKind.DOMAIN, "CAD fixed prefix exceeds the variable count")
    for value in fixed:
        if not _check_real(value, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid fixed CAD coordinate")


def _copy_fixed(fixed: List[RealAlgebraic], mut budget: Budget) raises PolyError -> List[RealAlgebraic]:
    _allocation[RealAlgebraic](len(fixed), budget)
    var result = List[RealAlgebraic](capacity=len(fixed))
    for value in fixed:
        result.append(_copy_real(value, budget))
    return result^


def _selection(stack: CADStackCertificate, fixed: List[RealAlgebraic], mut budget: Budget) raises PolyError -> Tuple[Int, Int]:
    if stack.level < 0:
        raise PolyError(PolyErrorKind.DOMAIN, "negative CAD stack level")
    if stack.level >= len(fixed):
        return (0, len(stack.samples))
    var position = _position(stack.roots, fixed[stack.level], budget)
    # A fixed prefix denotes a zero-dimensional base. A defining polynomial
    # must occur in the projection; a point inside a sector is not sufficient.
    if position % 2 == 0 or position >= len(stack.samples):
        raise PolyError(PolyErrorKind.DOMAIN, "fixed CAD coordinate is not a certified section")
    # Most samples refine the supplied witness directly. Other valid witnesses
    # can use a different defining polynomial or a narrower isolating interval.
    # Bind the value with an independently checked equality in those cases.
    if not _same_root(stack.samples[position], fixed[stack.level], budget):
        var equality = algebraic_compare(stack.samples[position], fixed[stack.level], budget)
        if equality.relation != 0 or not _check_comparison(stack.samples[position], fixed[stack.level], equality, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "fixed CAD section has a different exact value")
    return (position, position + 1)
