"""Certified feasibility and cylindrical coverings for polynomial sign conjunctions."""
from .errors import Budget, PolyError, _allocation
from .coefficients import QQ
from .sparse import SparsePoly
from .cad_types import CADCoveringCertificate
from .cad import _cad, _check_cad
from .algebraic_types import RealAlgebraic
from ._cad_common import _masks, _violation


def cad_covering(inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CADCoveringCertificate:
    """Decide a conjunction using a satisfying cell or a complete violating-cell covering.

    Args:
        inputs: The rational constraint polynomials.
        masks: Allowed sign bits per polynomial: one for negative, two for zero, four for positive.
        variables: Distinct variables in base-to-top cylindrical order.
        budget: The shared work allowance.
        lazard: Whether to use Lazard projection and lifting.

    Returns:
        A satisfying terminal cell, or one violated constraint for every cell covering real space.

    Raises:
        PolyError: UNSUPPORTED for disallowed McCallum nullification, DOMAIN for invalid masks or inputs, LIMIT or CAPACITY on work limits.
    """
    return _cad_covering(inputs, masks, variables, List[RealAlgebraic](), budget, lazard)


def _cad_covering(inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], fixed: List[RealAlgebraic], mut budget: Budget, lazard: Bool) raises PolyError -> CADCoveringCertificate:
    _masks(inputs, masks, budget)
    var decomposition = _cad(inputs, variables, fixed, budget, lazard)
    var violations = List[Int]()
    for i in range(len(decomposition.cells)):
        budget.charge(1)
        if decomposition.cells[i].level != len(variables) - 1:
            continue
        var violation = _violation(decomposition.cells[i], masks, budget)
        if violation == -1:
            return CADCoveringCertificate(decomposition^, i, List[Int]())
        _allocation[Int](1, budget)
        violations.append(violation)
    return CADCoveringCertificate(decomposition^, -1, violations^)


def check_cad_covering(inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CADCoveringCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check a feasibility witness or a complete sign-conflict covering of real space.

    Args:
        inputs: The original polynomial constraints.
        masks: The original allowed-sign masks.
        certificate: The proposed satisfying cell or complete violating-cell covering.
        budget: The shared work allowance.

    Returns:
        Whether the stated feasibility result follows from the certified decomposition.

    Raises:
        PolyError: DOMAIN for invalid masks, LIMIT or CAPACITY on work limits.
    """
    from .errors import PolyErrorKind
    _masks(inputs, masks, budget)
    try:
        return _check_cad_covering(inputs, masks, certificate, List[RealAlgebraic](), budget)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


def _check_cad_covering(inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CADCoveringCertificate, fixed: List[RealAlgebraic], mut budget: Budget) raises PolyError -> Bool:
    _masks(inputs, masks, budget)
    ref decomposition = certificate.decomposition
    if not _check_cad(inputs, decomposition, fixed, budget):
        return False
    var level = len(decomposition.projection.variables) - 1
    if certificate.witness_cell != -1:
        var i = certificate.witness_cell
        return i >= 0 and i < len(decomposition.cells) and decomposition.cells[i].level == level and len(certificate.violations) == 0 and _violation(decomposition.cells[i], masks, budget) == -1
    var cursor = 0
    for cell in decomposition.cells:
        budget.charge(1)
        if cell.level != level:
            continue
        if cursor >= len(certificate.violations):
            return False
        var i = certificate.violations[cursor]
        if i < 0 or i >= len(masks) or masks[i] & (UInt8(1) << UInt8(cell.signs[i].relation + 1)) != 0:
            return False
        cursor += 1
    return cursor == len(certificate.violations)
