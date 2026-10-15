"""Incremental numerical witnesses with independently bounded rational checking."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .coefficients import CoefficientRing
from .dense import DensePoly
from .sparse import SparsePoly
from .numeric_dense import NumericDensePoly
from .numeric_sparse import NumericSparsePoly
from .numeric_domains import NumericDomain, _context_cost, _numeric_cost
from .numeric_bounds import NumericBox, numeric_bounds, _box_valid, _box_add, _box_mul, _evaluation_bounds
from ._numeric_program import _NumericProgram, _dense_program
from ._numeric_storage import _NumericStorage, _storage_plan
from .numeric_evaluate import _require_enclosure, _exact_dense_boxes, _numeric_dense_boxes, _exact_sparse_program, _numeric_sparse_program, _instruction_value, _enclosure_value, _enclosure_header, _enclosure_result


@fieldwise_init
struct StreamedEnclosure[D: NumericDomain](Copyable, Movable):
    """A completed scalar receipt requiring its separately transported witness.

    Parameters:
        D: BallDomain or ComplexBallDomain.
    """
    var precision: Int
    """The midpoint precision in bits."""
    var steps: Int
    """The number of endpoint boxes in the complete stream."""
    var _value: List[Self.D.Scalar]

    def value(self, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Copy the completed scalar enclosure.

        Args:
            budget: The shared allowance for copying the scalar.

        Returns:
            An independently owned value; its receipt alone is not a proof.

        Raises:
            PolyError: DOMAIN for malformed storage or LIMIT on exhaustion.
        """
        return _enclosure_value[Self.D](self._value, budget)


@fieldwise_init
struct EnclosureStream[D: NumericDomain](Movable):
    """Produce the existing instruction witness one endpoint box at a time.

    Parameters:
        D: BallDomain or ComplexBallDomain.
    """
    var _program: _NumericProgram
    var _storage: _NumericStorage
    var _values: List[Self.D.Scalar]
    var _context: List[Self.D.Context]
    var _position: Int

    @staticmethod
    def from_dense[C: CoefficientRing](poly: DensePoly[C], argument: Self.D.Scalar, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental producer for an exact dense polynomial.

        Parameters:
            C: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            argument: The bounded evaluation argument.
            context: The explicit output precision.
            budget: The shared preparation allowance.

        Returns:
            A producer owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _dense_program(_exact_dense_boxes(poly, budget), numeric_bounds[Self.D](argument, budget), budget)
        return Self._prepare(program^, context, budget)

    @staticmethod
    def from_sparse[C: CoefficientRing](poly: SparsePoly[C], assignments: List[Tuple[Int64, Self.D.Scalar]], context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental producer for an exact sparse polynomial.

        Parameters:
            C: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            assignments: Bounded values for every occurring variable.
            context: The explicit output precision.
            budget: The shared preparation allowance.

        Returns:
            A producer owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _exact_sparse_program[Self.D](poly, assignments, budget)
        return Self._prepare(program^, context, budget)

    @staticmethod
    def from_numeric[E: NumericDomain](poly: NumericDensePoly[E], argument: Self.D.Scalar, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental producer for a numerical dense polynomial.

        Parameters:
            E: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            argument: The bounded evaluation argument.
            context: The explicit output precision.
            budget: The shared preparation allowance.

        Returns:
            A producer owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _dense_program(_numeric_dense_boxes(poly, budget), numeric_bounds[Self.D](argument, budget), budget)
        return Self._prepare(program^, context, budget)

    @staticmethod
    def from_numeric_sparse[E: NumericDomain](poly: NumericSparsePoly[E], assignments: List[Tuple[Int64, Self.D.Scalar]], context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental producer for a numerical sparse polynomial.

        Parameters:
            E: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            assignments: Bounded values for every occurring variable.
            context: The explicit output precision.
            budget: The shared preparation allowance.

        Returns:
            A producer owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _numeric_sparse_program[Self.D](poly, assignments, budget)
        return Self._prepare(program^, context, budget)

    @staticmethod
    def _prepare(var program: _NumericProgram, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        _require_enclosure[Self.D]()
        _ = _context_cost[Self.D](context)
        var storage = _storage_plan[True](program, budget)
        _allocation[Self.D.Scalar](storage.count, budget)
        _allocation[Self.D.Context](1, budget)
        var values = List[Self.D.Scalar](capacity=storage.count)
        return Self(program^, storage^, values^, [context.copy()], 0)

    def step_count(self) -> Int:
        """Return the full witness length, including unconsumed steps.

        Returns:
            The number of endpoint boxes needed for completion.
        """
        return len(self._program.instructions)

    def next(mut self, mut budget: Budget) raises PolyError -> Optional[NumericBox]:
        """Produce the next bound while recycling dead scalar slots.

        Args:
            budget: The shared evaluation allowance.

        Returns:
            The next original witness box, or None after the final step.

        Raises:
            PolyError: DOMAIN for invalid arithmetic, LIMIT or CAPACITY on work limits; a failed step is not consumed.
        """
        if self._position == self.step_count():
            return None
        budget.charge(1)
        var slot = self._storage.slots[self._position]
        var box: NumericBox
        if len(self._program.aliases) and self._program.aliases[self._position] >= 0:
            box = _evaluation_bounds[Self.D](self._values[slot], budget)
        else:
            var value = _instruction_value[Self.D, True](self._program, self._storage, self._values, self._position, self._context[0], budget)
            box = _evaluation_bounds[Self.D](value, budget)
            if slot == len(self._values):
                self._values.append(value^)
            else:
                self._values[slot] = value^
        self._position += 1
        return Optional(box^)

    def finish(self, mut budget: Budget) raises PolyError -> StreamedEnclosure[Self.D]:
        """Create a receipt only after every endpoint has been produced.

        Args:
            budget: The allowance for copying the final scalar.

        Returns:
            A completed receipt whose witness is the consumed endpoint stream.

        Raises:
            PolyError: DOMAIN for an incomplete stream or LIMIT on exhaustion.
        """
        if self._position != self.step_count():
            raise PolyError(PolyErrorKind.DOMAIN, "incomplete enclosure stream")
        ref value = self._values[self._storage.slots[self._program.result]]
        _charge_size(_numeric_cost[Self.D](value), budget)
        _allocation[Self.D.Scalar](1, budget)
        return StreamedEnclosure[Self.D](_context_cost[Self.D](self._context[0]), self.step_count(), [value.copy()])


@fieldwise_init
struct EnclosureStreamChecker[D: NumericDomain](Movable):
    """Check streamed bounds independently while retaining only live witnesses.

    Parameters:
        D: BallDomain or ComplexBallDomain.
    """
    var _program: _NumericProgram
    var _storage: _NumericStorage
    var _bounds: List[NumericBox]
    var _position: Int
    var _failed: Bool

    @staticmethod
    def from_dense[C: CoefficientRing](poly: DensePoly[C], argument: Self.D.Scalar, mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental checker for an exact dense polynomial.

        Parameters:
            C: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            argument: The bounded evaluation argument.
            budget: The shared preparation allowance.

        Returns:
            A checker owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _dense_program(_exact_dense_boxes(poly, budget), numeric_bounds[Self.D](argument, budget), budget)
        return Self._prepare(program^, budget)

    @staticmethod
    def from_sparse[C: CoefficientRing](poly: SparsePoly[C], assignments: List[Tuple[Int64, Self.D.Scalar]], mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental checker for an exact sparse polynomial.

        Parameters:
            C: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            assignments: Bounded values for every occurring variable.
            budget: The shared preparation allowance.

        Returns:
            A checker owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _exact_sparse_program[Self.D](poly, assignments, budget)
        return Self._prepare(program^, budget)

    @staticmethod
    def from_numeric[E: NumericDomain](poly: NumericDensePoly[E], argument: Self.D.Scalar, mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental checker for an numerical dense polynomial.

        Parameters:
            E: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            argument: The bounded evaluation argument.
            budget: The shared preparation allowance.

        Returns:
            A checker owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _dense_program(_numeric_dense_boxes(poly, budget), numeric_bounds[Self.D](argument, budget), budget)
        return Self._prepare(program^, budget)

    @staticmethod
    def from_numeric_sparse[E: NumericDomain](poly: NumericSparsePoly[E], assignments: List[Tuple[Int64, Self.D.Scalar]], mut budget: Budget) raises PolyError -> Self:
        """Prepare an incremental checker for an numerical sparse polynomial.

        Parameters:
            E: The source coefficient domain.

        Args:
            poly: The original polynomial, borrowed during preparation.
            assignments: Bounded values for every occurring variable.
            budget: The shared preparation allowance.

        Returns:
            A checker owning its prepared input graph.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits.
        """
        _require_enclosure[Self.D]()
        var program = _numeric_sparse_program[Self.D](poly, assignments, budget)
        return Self._prepare(program^, budget)

    @staticmethod
    def _prepare(var program: _NumericProgram, mut budget: Budget) raises PolyError -> Self:
        _require_enclosure[Self.D]()
        var storage = _storage_plan(program, budget)
        _allocation[NumericBox](storage.count, budget)
        var bounds = List[NumericBox](capacity=storage.count)
        return Self(program^, storage^, bounds^, 0, False)

    def step_count(self) -> Int:
        """Return the number of bounds required by the original input graph.

        Returns:
            The full witness length.
        """
        return len(self._program.instructions)

    def accept(mut self, var bounds: NumericBox, mut budget: Budget) raises PolyError -> Bool:
        """Verify one next bound with exact rational arithmetic.

        Args:
            bounds: The next owned endpoint box in instruction order.
            budget: The shared checking allowance.

        Returns:
            True for a valid next box; False permanently invalidates this checker.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits; a failed computation does not advance the checker.
        """
        budget.charge(1)
        if self._failed or self._position == self.step_count() or not _box_valid(bounds, budget):
            self._failed = True
            return False
        ref instruction = self._program.instructions[self._position]
        var expected: NumericBox
        if instruction.operation == 0:
            expected = self._program.inputs[instruction.left].copy()
        elif instruction.operation == 1:
            expected = _box_add(self._bounds[self._storage.slots[instruction.left]], self._bounds[self._storage.slots[instruction.right]], budget)
        else:
            expected = _box_mul(self._bounds[self._storage.slots[instruction.left]], self._bounds[self._storage.slots[instruction.right]], budget)
        if not bounds.contains(expected, budget):
            self._failed = True
            return False
        var slot = self._storage.slots[self._position]
        if slot == len(self._bounds):
            self._bounds.append(bounds^)
        else:
            self._bounds[slot] = bounds^
        self._position += 1
        return True

    def finish(self, receipt: StreamedEnclosure[Self.D], mut budget: Budget) raises PolyError -> Bool:
        """Verify stream completeness, reported precision and final enclosure.

        Args:
            receipt: The completed producer receipt to verify.
            budget: The allowance for validating the final scalar.

        Returns:
            True only for a complete accepted stream and a consistent receipt.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        budget.charge(1)
        if self._failed or self._position != self.step_count() or receipt.steps != self.step_count() or not _enclosure_header[Self.D](receipt.precision, receipt._value):
            return False
        return _enclosure_result[Self.D](receipt._value, self._bounds[self._storage.slots[self._program.result]], budget)
