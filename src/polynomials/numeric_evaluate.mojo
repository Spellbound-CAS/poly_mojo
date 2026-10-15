"""Ball evaluation with independent rational-bound certificates and bounded refinement."""
from apn_mojo import BallContext, Ball, ComplexBall
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .coefficients import CoefficientRing
from .dense import DensePoly
from .sparse import SparsePoly
from .monomial import Monomial
from .numeric_domains import NumericDomain, BallDomain, ComplexBallDomain, _numeric_apply, _context_cost, _numeric_cost
from .numeric_dense import NumericDensePoly
from .numeric_sparse import NumericSparsePoly
from .numeric_bounds import NumericBox, numeric_bounds, numeric_value, _box_valid, _box_add, _box_mul, _box_charge, _evaluation_bounds
from .numeric_convert import _exact_box, _require_exact_embedding
from .numeric_options import NumericOptions, NumericOutcome, _accurate, _next_precision
from ._numeric_program import _NumericProgram, _dense_program, _sparse_program
from ._numeric_storage import _NumericStorage, _storage_plan


def _enclosure_value[D: NumericDomain](values: List[D.Scalar], mut budget: Budget) raises PolyError -> D.Scalar:
    if len(values) != 1:
        raise PolyError(PolyErrorKind.DOMAIN, "malformed enclosure value storage")
    _charge_size(_numeric_cost[D](values[0]), budget)
    return values[0].copy()


def _enclosure_header[D: NumericDomain](precision: Int, values: List[D.Scalar]) -> Bool:
    if len(values) != 1 or precision < 2:
        return False
    comptime if D == BallDomain:
        return rebind[Ball](values[0]).midpoint().precision() == precision
    elif D == ComplexBallDomain:
        var scalar = rebind[ComplexBall](values[0])
        return scalar.real().midpoint().precision() == precision and scalar.imag().midpoint().precision() == precision
    else:
        return False


def _enclosure_result[D: NumericDomain](values: List[D.Scalar], bounds: NumericBox, mut budget: Budget) raises PolyError -> Bool:
    try:
        return numeric_bounds[D](_enclosure_value[D](values, budget), budget).contains(bounds, budget)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


@always_inline
def _instruction_value[D: NumericDomain, compact: Bool](program: _NumericProgram, storage: _NumericStorage, values: List[D.Scalar], i: Int, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    ref instruction = program.instructions[i]
    if instruction.operation == 0:
        return numeric_value[D](program.inputs[instruction.left], context, budget)
    comptime if not compact:
        if instruction.left < 0 or instruction.right < 0 or instruction.left >= i or instruction.right >= i:
            raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical program operand")
    return _numeric_apply[D](values[storage.slots[instruction.left] if compact else instruction.left], values[storage.slots[instruction.right] if compact else instruction.right], 0 if instruction.operation == 1 else 2, context, budget)


@fieldwise_init
struct EnclosureEvaluation[D: NumericDomain](Copyable, Movable):
    """A ball value with an exact endpoint witness for every evaluation step.

    Parameters:
        D: BallDomain or ComplexBallDomain.
    """
    var precision: Int
    """The working midpoint precision in bits."""
    # A one-element owning list avoids Mojo 1.1 partial-field lifetime bugs
    # for structs with directly stored associated scalar types.
    var _value: List[Self.D.Scalar]

    def value(self, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Copy the reported enclosure.

        Args:
            budget: The shared work allowance for copying the scalar.

        Returns:
            The enclosing scalar, independently owned by the caller.

        Raises:
            PolyError: DOMAIN for malformed storage, LIMIT or CAPACITY on work limits.
        """
        return _enclosure_value[Self.D](self._value, budget)
    var steps: List[NumericBox]
    """Bounds in the deterministic evaluation program's instruction order."""


@fieldwise_init
struct AdaptiveEvaluation[D: NumericDomain](Copyable, Movable):
    """An accuracy outcome accompanied by a checkable enclosure.

    Parameters:
        D: BallDomain or ComplexBallDomain.
    """
    var outcome: NumericOutcome
    """SUCCESS only when the requested width was established."""
    var evaluation: EnclosureEvaluation[Self.D]
    """The last rigorous enclosure, including for an inconclusive accuracy request."""


def _require_enclosure[D: NumericDomain]() raises PolyError:
    comptime if D != BallDomain and D != ComplexBallDomain:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "enclosure evaluation requires a ball domain")


def _evaluate_program[D: NumericDomain](program: _NumericProgram, context: D.Context, mut budget: Budget) raises PolyError -> EnclosureEvaluation[D]:
    if len(program.instructions) < 64:
        return _evaluate_scheduled[D, False](program, context, budget)
    if len(program.aliases):
        return _evaluate_scheduled[D, True, True](program, context, budget)
    return _evaluate_scheduled[D, True](program, context, budget)


def _evaluate_scheduled[D: NumericDomain, compact: Bool, reuse: Bool = False](program: _NumericProgram, context: D.Context, mut budget: Budget) raises PolyError -> EnclosureEvaluation[D]:
    _require_enclosure[D]()
    var precision = _context_cost[D](context)
    if program.result < 0 or program.result >= len(program.instructions):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical program result")
    var storage: _NumericStorage
    comptime if compact:
        storage = _storage_plan[reuse](program, budget)
    else:
        storage = _NumericStorage(List[Int](), len(program.instructions))
    _allocation[D.Scalar](storage.count, budget)
    _allocation[NumericBox](len(program.instructions), budget)
    var values = List[D.Scalar](capacity=storage.count)
    var bounds = List[NumericBox](capacity=len(program.instructions))
    for i in range(len(program.instructions)):
        budget.charge(1)
        comptime if reuse:
            if program.aliases[i] >= 0:
                _box_charge(bounds[program.aliases[i]], budget)
                bounds.append(bounds[program.aliases[i]].copy())
                continue
        var value = _instruction_value[D, compact](program, storage, values, i, context, budget)
        bounds.append(_evaluation_bounds[D](value, budget))
        comptime if compact:
            var slot = storage.slots[i]
            if slot == len(values):
                values.append(value^)
            else:
                values[slot] = value^
        else:
            values.append(value^)
    _allocation[D.Scalar](1, budget)
    var result = storage.slots[program.result] if compact else program.result
    # Transfer the owned result out of the workspace; copying an arbitrary-
    # precision scalar here adds allocation work to every small root evaluation.
    var value = values.pop()
    if result < len(values):
        swap(value, values[result])
    return EnclosureEvaluation[D](precision, [value^], bounds^)


def _check_program[D: NumericDomain](program: _NumericProgram, certificate: EnclosureEvaluation[D], mut budget: Budget) raises PolyError -> Bool:
    _require_enclosure[D]()
    if not _enclosure_header[D](certificate.precision, certificate._value) or len(certificate.steps) != len(program.instructions):
        return False
    for i in range(len(program.instructions)):
        budget.charge(1)
        ref instruction = program.instructions[i]
        if not _box_valid(certificate.steps[i], budget):
            return False
        var expected: NumericBox
        if instruction.operation == 0:
            expected = program.inputs[instruction.left].copy()
        elif instruction.operation == 1:
            expected = _box_add(certificate.steps[instruction.left], certificate.steps[instruction.right], budget)
        else:
            expected = _box_mul(certificate.steps[instruction.left], certificate.steps[instruction.right], budget)
        if not certificate.steps[i].contains(expected, budget):
            return False
    return _enclosure_result[D](certificate._value, certificate.steps[program.result], budget)


def _enclosure_context[D: NumericDomain](bits: Int, ceiling: Int) raises PolyError -> D.Context:
    _require_enclosure[D]()
    try:
        comptime if D == BallDomain or D == ComplexBallDomain:
            return rebind_var[D.Context](BallContext(bits, max_precision=ceiling))
        else:
            raise Error("not an enclosure domain")
    except e:
        raise PolyError(PolyErrorKind.DOMAIN, String(e))


def _adaptive_program[D: NumericDomain](program: _NumericProgram, options: NumericOptions, mut budget: Budget) raises PolyError -> AdaptiveEvaluation[D]:
    options._validate()
    var bits = options.initial_precision
    while True:
        budget.charge(1)
        var proof = _evaluate_program[D](program, _enclosure_context[D](bits, options.max_precision), budget)
        if _accurate(numeric_bounds[D](proof.value(budget), budget), options, budget):
            return AdaptiveEvaluation[D](NumericOutcome.SUCCESS, proof^)
        if bits == options.max_precision:
            return AdaptiveEvaluation[D](NumericOutcome.PRECISION_CEILING, proof^)
        bits = _next_precision(bits, options.max_precision)


def _exact_dense_boxes[C: CoefficientRing](poly: DensePoly[C], mut budget: Budget) raises PolyError -> List[NumericBox]:
    _require_exact_embedding[C]()
    poly.check_invariants(budget)
    _allocation[NumericBox](len(poly.coefficients), budget)
    var values = List[NumericBox](capacity=len(poly.coefficients))
    for value in poly.coefficients:
        values.append(_exact_box(value, budget))
    return values^


def _numeric_dense_boxes[D: NumericDomain](poly: NumericDensePoly[D], mut budget: Budget) raises PolyError -> List[NumericBox]:
    _allocation[NumericBox](len(poly.coefficients), budget)
    var values = List[NumericBox](capacity=len(poly.coefficients))
    for value in poly.coefficients:
        values.append(numeric_bounds[D](value, budget))
    return values^


def _assignment_boxes[D: NumericDomain](assignments: List[Tuple[Int64, D.Scalar]], mut budget: Budget) raises PolyError -> List[Tuple[Int64, NumericBox]]:
    _allocation[Tuple[Int64, NumericBox]](len(assignments), budget)
    var result = List[Tuple[Int64, NumericBox]](capacity=len(assignments))
    for item in assignments:
        result.append((item[0], numeric_bounds[D](item[1], budget)))
    return result^


def _exact_sparse_program[D: NumericDomain, C: CoefficientRing](poly: SparsePoly[C], assignments: List[Tuple[Int64, D.Scalar]], mut budget: Budget) raises PolyError -> _NumericProgram:
    _require_exact_embedding[C]()
    poly.check_invariants(budget)
    _allocation[NumericBox](len(poly.terms), budget)
    _allocation[Monomial](len(poly.terms), budget)
    var coefficients = List[NumericBox]()
    var monomials = List[Monomial]()
    for term in poly.terms:
        _charge_size(len(term.monomial.pairs), budget)
        coefficients.append(_exact_box(term.coefficient, budget))
        monomials.append(term.monomial.copy())
    return _sparse_program(coefficients, monomials, _assignment_boxes[D](assignments, budget), budget)


def _numeric_sparse_program[D: NumericDomain, E: NumericDomain](poly: NumericSparsePoly[E], assignments: List[Tuple[Int64, D.Scalar]], mut budget: Budget) raises PolyError -> _NumericProgram:
    _allocation[NumericBox](len(poly.terms), budget)
    _allocation[Monomial](len(poly.terms), budget)
    var coefficients = List[NumericBox]()
    var monomials = List[Monomial]()
    for term in poly.terms:
        poly.order.validate(term.monomial, budget)
        coefficients.append(numeric_bounds[E](term.coefficient, budget))
        monomials.append(term.monomial.copy())
    return _sparse_program(coefficients, monomials, _assignment_boxes[D](assignments, budget), budget)


def enclose_dense[D: NumericDomain, C: CoefficientRing](poly: DensePoly[C], argument: D.Scalar, context: D.Context, mut budget: Budget) raises PolyError -> EnclosureEvaluation[D]:
    """Evaluate an exact ZZ/QQ/QQi polynomial with a rigorous ball certificate.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The exact polynomial.
        argument: The bounded argument enclosure.
        context: Explicit output precision.
        budget: The shared work allowance.

    Returns:
        The enclosing value and rational endpoint trace.

    Raises:
        PolyError: DOMAIN for invalid values, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
    """
    return _evaluate_program[D](_dense_program(_exact_dense_boxes(poly, budget), numeric_bounds[D](argument, budget), budget), context, budget)


def check_enclose_dense[D: NumericDomain, C: CoefficientRing](poly: DensePoly[C], argument: D.Scalar, certificate: EnclosureEvaluation[D], mut budget: Budget) raises PolyError -> Bool:
    """Verify exact dense evaluation bounds without floating-point arithmetic or search.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The exact polynomial.
        argument: The original bounded argument.
        certificate: The proposed evaluation witness.
        budget: The shared work allowance.

    Returns:
        Whether every local bound and the reported value are justified.

    Raises:
        PolyError: DOMAIN for invalid inputs, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
    """
    return _check_program[D](_dense_program(_exact_dense_boxes(poly, budget), numeric_bounds[D](argument, budget), budget), certificate, budget)


def enclose_numeric[D: NumericDomain, E: NumericDomain](poly: NumericDensePoly[E], argument: D.Scalar, context: D.Context, mut budget: Budget) raises PolyError -> EnclosureEvaluation[D]:
    """Enclose evaluation of every polynomial in a numerical coefficient family.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        E: The source coefficient domain.

    Args:
        poly: The numerical point or family polynomial.
        argument: The bounded argument enclosure.
        context: Explicit output precision.
        budget: The shared work allowance.

    Returns:
        A family-wide enclosure and its bound trace.

    Raises:
        PolyError: DOMAIN for invalid values, UNSUPPORTED for other output domains, LIMIT or CAPACITY on work limits.
    """
    return _evaluate_program[D](_dense_program(_numeric_dense_boxes(poly, budget), numeric_bounds[D](argument, budget), budget), context, budget)


def check_enclose_numeric[D: NumericDomain, E: NumericDomain](poly: NumericDensePoly[E], argument: D.Scalar, certificate: EnclosureEvaluation[D], mut budget: Budget) raises PolyError -> Bool:
    """Check a numerical family evaluation using exact endpoint arithmetic.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        E: The source coefficient domain.

    Args:
        poly: The original numerical polynomial.
        argument: The original bounded argument.
        certificate: The proposed evaluation witness.
        budget: The shared work allowance.

    Returns:
        Whether the witness encloses all permitted evaluations.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _check_program[D](_dense_program(_numeric_dense_boxes(poly, budget), numeric_bounds[D](argument, budget), budget), certificate, budget)


def enclose_sparse[D: NumericDomain, C: CoefficientRing](poly: SparsePoly[C], assignments: List[Tuple[Int64, D.Scalar]], context: D.Context, mut budget: Budget) raises PolyError -> EnclosureEvaluation[D]:
    """Evaluate an exact multivariate polynomial with a checkable enclosure trace.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The exact sparse polynomial.
        assignments: Distinct assignments covering every occurring variable.
        context: Explicit output precision.
        budget: The shared work allowance.

    Returns:
        The enclosing value and local bound witnesses.

    Raises:
        PolyError: DOMAIN for invalid assignments, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
    """
    return _evaluate_program[D](_exact_sparse_program[D](poly, assignments, budget), context, budget)


def check_enclose_sparse[D: NumericDomain, C: CoefficientRing](poly: SparsePoly[C], assignments: List[Tuple[Int64, D.Scalar]], certificate: EnclosureEvaluation[D], mut budget: Budget) raises PolyError -> Bool:
    """Verify sparse exact-polynomial evaluation bounds independently.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The original exact polynomial.
        assignments: The original argument enclosures.
        certificate: The proposed evaluation witness.
        budget: The shared work allowance.

    Returns:
        Whether the complete bound trace is valid.

    Raises:
        PolyError: DOMAIN for invalid assignments, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
    """
    return _check_program[D](_exact_sparse_program[D](poly, assignments, budget), certificate, budget)


def enclose_numeric_sparse[D: NumericDomain, E: NumericDomain](poly: NumericSparsePoly[E], assignments: List[Tuple[Int64, D.Scalar]], context: D.Context, mut budget: Budget) raises PolyError -> EnclosureEvaluation[D]:
    """Enclose sparse numerical-family evaluation with cached variable powers.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        E: The source numerical domain.

    Args:
        poly: The numerical sparse polynomial.
        assignments: Distinct assignments covering every occurring variable.
        context: Explicit output precision.
        budget: The shared work allowance.

    Returns:
        A family-wide enclosure and its bound trace.

    Raises:
        PolyError: DOMAIN for invalid assignments, LIMIT or CAPACITY on work limits.
    """
    return _evaluate_program[D](_numeric_sparse_program[D](poly, assignments, budget), context, budget)


def check_enclose_numeric_sparse[D: NumericDomain, E: NumericDomain](poly: NumericSparsePoly[E], assignments: List[Tuple[Int64, D.Scalar]], certificate: EnclosureEvaluation[D], mut budget: Budget) raises PolyError -> Bool:
    """Verify sparse numerical-family evaluation using rational bounds.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        E: The source numerical domain.

    Args:
        poly: The original numerical polynomial.
        assignments: The original argument enclosures.
        certificate: The proposed evaluation witness.
        budget: The shared work allowance.

    Returns:
        Whether every local bound is justified.

    Raises:
        PolyError: DOMAIN for invalid assignments, LIMIT or CAPACITY on work limits.
    """
    return _check_program[D](_numeric_sparse_program[D](poly, assignments, budget), certificate, budget)


def adaptive_enclose_dense[D: NumericDomain, C: CoefficientRing](poly: DensePoly[C], argument: D.Scalar, options: NumericOptions, mut budget: Budget) raises PolyError -> AdaptiveEvaluation[D]:
    """Refine evaluation precision from the original exact coefficients.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The exact polynomial.
        argument: The original bounded argument, preserving its uncertainty.
        options: The precision ceiling and requested widths.
        budget: The shared work allowance.

    Returns:
        A certified enclosure with SUCCESS or an explicit precision-ceiling outcome.

    Raises:
        PolyError: DOMAIN for invalid options or inputs, LIMIT or CAPACITY on work limits.
    """
    return _adaptive_program[D](_dense_program(_exact_dense_boxes(poly, budget), numeric_bounds[D](argument, budget), budget), options, budget)


def adaptive_enclose_sparse[D: NumericDomain, C: CoefficientRing](poly: SparsePoly[C], assignments: List[Tuple[Int64, D.Scalar]], options: NumericOptions, mut budget: Budget) raises PolyError -> AdaptiveEvaluation[D]:
    """Refine sparse evaluation while retaining exact coefficients and input uncertainty.

    Parameters:
        D: BallDomain or ComplexBallDomain.
        C: ZZ, QQ, or QQi.

    Args:
        poly: The exact sparse polynomial.
        assignments: Distinct variable/enclosure pairs.
        options: The precision ceiling and requested widths.
        budget: The shared work allowance.

    Returns:
        The final checkable enclosure and its accuracy outcome.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return _adaptive_program[D](_exact_sparse_program[D](poly, assignments, budget), options, budget)


def check_numeric_accuracy[D: NumericDomain](result: AdaptiveEvaluation[D], options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    """Check an adaptive result's accuracy claim separately from its evaluation witness.

    Parameters:
        D: BallDomain or ComplexBallDomain.

    Args:
        result: The adaptive result whose evaluation witness must also be checked.
        options: The originally requested accuracy and precision limits.
        budget: The shared work allowance.

    Returns:
        Whether its outcome agrees with the enclosure width and precision schedule.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    options._validate()
    if len(result.evaluation._value) != 1:
        return False
    var bits = options.initial_precision
    while bits < result.evaluation.precision and bits < options.max_precision:
        budget.charge(1)
        bits = _next_precision(bits, options.max_precision)
    if bits != result.evaluation.precision:
        return False
    var accurate = _accurate(numeric_bounds[D](result.evaluation.value(budget), budget), options, budget)
    if result.outcome == NumericOutcome.SUCCESS:
        return accurate
    return result.outcome == NumericOutcome.PRECISION_CEILING and not accurate and bits == options.max_precision
