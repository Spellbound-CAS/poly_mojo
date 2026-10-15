"""Resumable root searches with checked reuse of prior local certificates."""
from .coefficients import CoefficientRing, QQ
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .certificates import SturmSequence
from ._dense_ops import _copy
from .numeric_domains import NumericDomain, _numeric_cost
from .numeric_dense import NumericDensePoly
from .numeric_bounds import NumericBox, _box_charge
from .numeric_options import NumericOptions, NumericOutcome
from .numeric_evaluate import EnclosureEvaluation, _exact_dense_boxes, _numeric_dense_boxes
from ._numeric_roots import _root_degree, _rational_polynomial, _derivative_boxes
from .numeric_real_roots import RealRootCell, RealRootResult, _check_real, _real_roots, _real_cell, _refine_real_cells
from .numeric_complex_roots import ComplexRootRegion, ComplexRootResult, _check_complex, _complex_roots, _representative, _refine_complex_regions


def _charge_evaluation[D: NumericDomain](proof: EnclosureEvaluation[D], mut budget: Budget) raises PolyError:
    _allocation[D.Scalar](len(proof._value), budget)
    _allocation[NumericBox](len(proof.steps), budget)
    for value in proof._value:
        _charge_size(_numeric_cost[D](value), budget)
    for step in proof.steps:
        _box_charge(step, budget)


def _copy_cell(cell: RealRootCell, mut budget: Budget) raises PolyError -> RealRootCell:
    _box_charge(cell.interval, budget)
    _allocation[type_of(cell.evaluations[0])](len(cell.evaluations), budget)
    for proof in cell.evaluations:
        _charge_evaluation(proof, budget)
    return cell.copy()


def _copy_region(region: ComplexRootRegion, mut budget: Budget) raises PolyError -> ComplexRootRegion:
    _box_charge(region.rectangle, budget)
    _allocation[type_of(region.boundary[0])](len(region.boundary), budget)
    for segment in region.boundary:
        _box_charge(segment.start, budget)
        _box_charge(segment.end, budget)
        _charge_evaluation(segment.evaluation, budget)
    return region.copy()


def _copy_sturm(sequence: SturmSequence, mut budget: Budget) raises PolyError -> SturmSequence:
    _allocation[DensePoly[QQ]](len(sequence.polynomials), budget)
    _allocation[DensePoly[QQ]](len(sequence.quotients), budget)
    var polynomials = List[DensePoly[QQ]](capacity=len(sequence.polynomials))
    var quotients = List[DensePoly[QQ]](capacity=len(sequence.quotients))
    for p in sequence.polynomials:
        polynomials.append(_copy(p, budget))
    for q in sequence.quotients:
        quotients.append(_copy(q, budget))
    return SturmSequence(polynomials^, quotients^)


def _resume_real(coefficients: List[NumericBox], previous: RealRootResult, previous_options: NumericOptions, options: NumericOptions, mut budget: Budget, complete_only: Bool) raises PolyError -> RealRootResult:
    previous_options._validate()
    options._validate()
    if not len(previous.cells):
        # No evidence exists to reuse. Retrying an explicitly inconclusive empty
        # result recomputes all degree/family obligations with the new allowance.
        if previous.complete or previous.accurate or previous.isolated or previous.sequence or (previous.outcome != NumericOutcome.UNCERTAIN_DEGREE and previous.outcome != NumericOutcome.INPUT_UNCERTAINTY):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid empty real-root result")
        return _real_roots(coefficients, options, budget, complete_only)
    if not _check_real(coefficients, previous, previous_options, budget, allow_partial=True):
        raise PolyError(PolyErrorKind.DOMAIN, "previous real-root evidence is invalid")
    var exact = _rational_polynomial(coefficients, budget)
    var degree = _root_degree(coefficients, budget)
    var derivative = _derivative_boxes(coefficients, budget)
    var sequence = Optional[SturmSequence]()
    if previous.sequence:
        sequence = _copy_sturm(previous.sequence.value(), budget)
    _allocation[RealRootCell](len(previous.cells), budget)
    var pending = List[RealRootCell](capacity=len(previous.cells))
    var i = len(previous.cells) - 1
    while i >= 0:
        budget.charge(1)
        if previous.cells[i].count >= 0:
            pending.append(_copy_cell(previous.cells[i], budget))
            i -= 1
        else:
            # Merge adjacent unresolved cells so a root on a previous cut can
            # become interior again. Certified neighboring boundaries stay fixed.
            var last = i
            while i > 0 and previous.cells[i - 1].count < 0:
                budget.charge(1)
                i -= 1
            var interval = NumericBox.real_interval(previous.cells[i].interval.real_low, previous.cells[last].interval.real_high, budget)
            pending.append(_real_cell(coefficients, derivative, exact, sequence, interval, options, budget))
            i -= 1
    return _refine_real_cells(coefficients, derivative, exact, sequence^, pending^, degree.value(), options, budget, complete_only)


def _resume_complex(coefficients: List[NumericBox], previous: ComplexRootResult, previous_options: NumericOptions, options: NumericOptions, mut budget: Budget, complete_only: Bool) raises PolyError -> ComplexRootResult:
    previous_options._validate()
    options._validate()
    if not previous.complete and not len(previous.regions):
        if previous.accurate or previous.isolated or (previous.outcome != NumericOutcome.UNCERTAIN_DEGREE and previous.outcome != NumericOutcome.UNRESOLVED_BOUNDARY):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid empty complex-root result")
        return _complex_roots(coefficients, options, budget, complete_only)
    if not _check_complex(coefficients, previous, previous_options, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "previous complex-root evidence is invalid")
    var degree = _root_degree(coefficients, budget)
    var representative = _representative(coefficients, budget)
    _allocation[ComplexRootRegion](len(previous.regions), budget)
    var pending = List[ComplexRootRegion](capacity=len(previous.regions))
    for i in range(len(previous.regions) - 1, -1, -1):
        pending.append(_copy_region(previous.regions[i], budget))
    var subdivisions = 0
    return _refine_complex_regions(coefficients, representative, pending^, degree.value(), options, subdivisions, budget, complete_only)


def refine_real_roots_exact[C: CoefficientRing](poly: DensePoly[C], previous: RealRootResult, previous_options: NumericOptions, options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> RealRootResult:
    """Resume real refinement after validating the previous root evidence.

    Parameters:
        C: ZZ or QQ.

    Args:
        poly: The original polynomial or coefficient family.
        previous: The prior result; certified satisfactory regions are retained.
        previous_options: Options describing the previous result's flags.
        options: New targets and an additional subdivision allowance.
        budget: The shared work allowance.
        complete_only: Raise unless every new requested guarantee succeeds.

    Returns:
        A refined cover; empty inconclusive results restart, and unresolved real cells merge before retrying.

    Raises:
        PolyError: DOMAIN for invalid prior evidence, LIMIT or CAPACITY on work limits.
    """
    return _resume_real(_exact_dense_boxes(poly, budget), previous, previous_options, options, budget, complete_only)


def refine_real_roots_numeric[D: NumericDomain](poly: NumericDensePoly[D], previous: RealRootResult, previous_options: NumericOptions, options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> RealRootResult:
    """Resume real refinement after validating the previous root evidence.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        poly: The original polynomial or coefficient family.
        previous: The prior result; certified satisfactory regions are retained.
        previous_options: Options describing the previous result's flags.
        options: New targets and an additional subdivision allowance.
        budget: The shared work allowance.
        complete_only: Raise unless every new requested guarantee succeeds.

    Returns:
        A refined cover; empty inconclusive results restart, and unresolved real cells merge before retrying.

    Raises:
        PolyError: DOMAIN for invalid prior evidence, LIMIT or CAPACITY on work limits.
    """
    return _resume_real(_numeric_dense_boxes(poly, budget), previous, previous_options, options, budget, complete_only)


def refine_complex_roots_exact[C: CoefficientRing](poly: DensePoly[C], previous: ComplexRootResult, previous_options: NumericOptions, options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> ComplexRootResult:
    """Resume complex refinement after validating the previous root evidence.

    Parameters:
        C: ZZ, QQ, or QQi.

    Args:
        poly: The original polynomial or coefficient family.
        previous: The prior result; certified satisfactory regions are retained.
        previous_options: Options describing the previous result's flags.
        options: New targets and an additional subdivision allowance.
        budget: The shared work allowance.
        complete_only: Raise unless every new requested guarantee succeeds.

    Returns:
        A refined cover; empty inconclusive results restart, and unresolved real cells merge before retrying.

    Raises:
        PolyError: DOMAIN for invalid prior evidence, LIMIT or CAPACITY on work limits.
    """
    return _resume_complex(_exact_dense_boxes(poly, budget), previous, previous_options, options, budget, complete_only)


def refine_complex_roots_numeric[D: NumericDomain](poly: NumericDensePoly[D], previous: ComplexRootResult, previous_options: NumericOptions, options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> ComplexRootResult:
    """Resume complex refinement after validating the previous root evidence.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        poly: The original polynomial or coefficient family.
        previous: The prior result; certified satisfactory regions are retained.
        previous_options: Options describing the previous result's flags.
        options: New targets and an additional subdivision allowance.
        budget: The shared work allowance.
        complete_only: Raise unless every new requested guarantee succeeds.

    Returns:
        A refined cover; empty inconclusive results restart, and unresolved real cells merge before retrying.

    Raises:
        PolyError: DOMAIN for invalid prior evidence, LIMIT or CAPACITY on work limits.
    """
    return _resume_complex(_numeric_dense_boxes(poly, budget), previous, previous_options, options, budget, complete_only)
