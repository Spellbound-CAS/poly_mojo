"""Independent checks for fixed-precision and adaptive enclosure evaluation."""
from std.memory import ArcPointer
from ._archive_meta import _Archive, _QuotientInput
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink, QQ
from std.builtin.rebind import downcast
from apn_mojo import BallContext
from polynomials import Budget, PolyError, PolyErrorKind, CoefficientRing, NumericDomain, BallDomain, ComplexBallDomain, DensePoly, SparsePoly, NumericDensePoly, NumericSparsePoly, EnclosureEvaluation, AdaptiveEvaluation, NumericOptions, NumericOutcome
from polynomials import enclose_dense, enclose_sparse, enclose_numeric, enclose_numeric_sparse, check_enclose_dense, check_enclose_sparse, check_enclose_numeric, check_enclose_numeric_sparse, adaptive_enclose_dense, adaptive_enclose_sparse, check_numeric_accuracy
from ._common import _allocation, _charge_size, _grow, _sum
from . import _numeric as _num
from . import _algebraic_numeric as _algnum
from ._algebraic_value import Real, Complex, _copy_witness
from polynomials import NoContext, numeric_bounds, NumericBox, RealAlgebraic, StreamedEnclosure, MonomialOrder

comptime _Core[C: AnyType, numerical: Bool, sparse: Bool] = (NumericSparsePoly[downcast[C, NumericDomain]] if sparse else NumericDensePoly[downcast[C, NumericDomain]]) if numerical else (SparsePoly[downcast[C, CoefficientRing]] if sparse else DensePoly[downcast[C, CoefficientRing]])
comptime _Domain[C: AnyType, X: AnyType] = ComplexBallDomain if _num._IsComplex[C] or _num._IsComplex[X] else BallDomain


@fieldwise_init
struct _Point[D: NumericDomain](Copyable, Movable):
    var axis: Int64
    var source: Complex
    var witness: _algnum._Witness


def _argument[D: NumericDomain, T: Copyable & Deinitable](source: T, axis: Int64, context: D.Context, mut witnesses: List[_Point[D]], mut budget: Budget) raises PolyError -> D.Scalar:
    comptime if T == Real or T == Complex:
        var value: Complex
        comptime if T == Real:
            value = Complex(rebind[Real](source).copy(), Real.zero(NoContext()))
        else:
            value = rebind[Complex](source).copy()
        _charge_size(value.size(), budget)
        var converted = _algnum._enclose[Complex, D](value, context, budget)
        _grow(witnesses, budget)
        witnesses.append(_Point[D](axis, value^, _algnum._Witness(_copy_witness(converted[1].real, budget), _copy_witness(converted[1].imag, budget))))
        _ = numeric_bounds[D](converted[0], budget)
        return converted[0].copy()
    else:
        return _num._value[D](source, context, budget)


@fieldwise_init
struct _Data[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](Copyable, Movable):
    var polynomial: _Core[Self.C, Self.numerical, Self.sparse]
    var points: List[Tuple[Int64, Self.D.Scalar]]
    var proof: AdaptiveEvaluation[Self.D]
    var options: Optional[NumericOptions]
    var algebraic: Optional[_algnum._Polynomial[Self.D, Self.sparse]]
    var algebraic_points: List[_Point[Self.D]]


@fieldwise_init
struct Enclosure[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](_Archive, _QuotientInput, Copyable, Movable):
    """An enclosing value with its original inputs and independently checkable trace.

    Parameters:
        C: The inferred source coefficient ring or numerical domain.
        D: The inferred real or complex ball domain.
        numerical: Whether the source coefficients are numerical.
        sparse: Whether the source is sparse.
    """
    var _data: ArcPointer[_Data[Self.C, Self.D, Self.numerical, Self.sparse]]

    def value(self, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Copy the enclosing scalar.

        Args:
            budget: The shared work allowance.

        Returns:
            The ball enclosing every evaluation allowed by the stored inputs.

        Raises:
            PolyError: LIMIT/CAPACITY on resource limits.
        """
        return self._data[].proof.evaluation.value(budget)

    def value(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Self.D.Scalar:
        """Copy the enclosing scalar with a work limit.

        Args:
            max_steps: The complete call's work allowance.

        Returns:
            An independent ball value.

        Raises:
            PolyError: LIMIT/CAPACITY on resource limits.
        """
        var budget = Budget(max_steps)
        return self.value(budget)

    def precision(self) -> Int:
        """Return the output midpoint precision in bits.

        Returns:
            The precision used for the last enclosure.
        """
        return self._data[].proof.evaluation.precision

    def outcome(self) -> Optional[NumericOutcome]:
        """Return the adaptive accuracy outcome, when accuracy was requested.

        Returns:
            SUCCESS or PRECISION_CEILING for adaptive calls; None for fixed precision.
        """
        if self._data[].options:
            return Optional(self._data[].proof.outcome)
        return Optional[NumericOutcome]()

    def accurate(self) -> Bool:
        """Report whether an adaptive accuracy request succeeded.

        Returns:
            True only for an adaptive SUCCESS; call check to verify the evidence.
        """
        if not self._data[].options:
            return False
        return self._data[].proof.outcome == NumericOutcome.SUCCESS

    def options(self) -> Optional[NumericOptions]:
        """Return the original accuracy request when present.

        Returns:
            The retained numerical options, or None for fixed precision.
        """
        return self._data[].options.copy()

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof.evaluation))] EnclosureEvaluation[Self.D]:
        """Borrow the evaluation trace.

        For algebraic coefficients, keep the complete result for conversion
        witnesses checked by evaluation.check; this accessor exposes only the trace.

        Returns:
            Immutable output and rational bounds for every arithmetic step.
        """
        return self._data[].proof.evaluation


    def _quotient_precision(self) -> Int:
        return self.precision()

    def _quotient_bounds(self, mut budget: Budget) raises PolyError -> NumericBox:
        return numeric_bounds[Self.D](self.value(budget), budget)

    def _same_point(self, other: Self, mut budget: Budget) raises PolyError -> Bool:
        ref left = self._data[]
        ref right = other._data[]
        if len(left.points) != len(right.points) or len(left.algebraic_points) != len(right.algebraic_points):
            return False
        for i in range(len(left.points)):
            if left.points[i][0] != right.points[i][0]: return False
            var x = numeric_bounds[Self.D](left.points[i][1], budget)
            var y = numeric_bounds[Self.D](right.points[i][1], budget)
            if not x.contains(y, budget) or not y.contains(x, budget): return False
        for i in range(len(left.algebraic_points)):
            if left.algebraic_points[i].axis != right.algebraic_points[i].axis: return False
            if not left.algebraic_points[i].source.sub(right.algebraic_points[i].source, budget).is_zero(): return False
        return True

    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io as codec
        comptime R = QQ if Self.numerical else downcast[Self.C, CoefficientRing]
        codec._depth(depth, budget)
        ref data = self._data[]
        codec._write_node[_Core[Self.C, Self.numerical, Self.sparse], R, Self.D, Int](data.polynomial, writer, sink, budget, depth + 1)
        writer.line(String(len(data.points)), sink, budget)
        for point in data.points:
            codec._write_node[Int64, R, Self.D, Int](point[0], writer, sink, budget, depth + 1)
            codec._write_node[Self.D.Scalar, R, Self.D, Int](point[1], writer, sink, budget, depth + 1)
        codec._write_node[AdaptiveEvaluation[Self.D], R, Self.D, Int](data.proof, writer, sink, budget, depth + 1)
        writer.line(String(Int(Bool(data.options))), sink, budget)
        if data.options:
            codec._write_node[NumericOptions, R, Self.D, Int](data.options.value(), writer, sink, budget, depth + 1)
        writer.line(String(Int(Bool(data.algebraic))), sink, budget)
        if data.algebraic:
            codec._write_node[_algnum._Polynomial[Self.D, Self.sparse], R, Self.D, Int](data.algebraic.value(), writer, sink, budget, depth + 1)
        writer.line(String(len(data.algebraic_points)), sink, budget)
        for point in data.algebraic_points:
            codec._write_node[Int64, R, Self.D, Int](point.axis, writer, sink, budget, depth + 1)
            codec._write_node[Complex, R, Self.D, Int](point.source, writer, sink, budget, depth + 1)
            codec._write_node[_algnum._Witness, R, Self.D, Int](point.witness, writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io as codec
        comptime R = QQ if Self.numerical else downcast[Self.C, CoefficientRing]
        codec._depth(depth, budget)
        var p = codec._read_node[_Core[Self.C, Self.numerical, Self.sparse], R, Self.D, Int](reader, budget, depth + 1)
        var points = List[Tuple[Int64, Self.D.Scalar]]()
        for _ in range(reader.count(budget)):
            var axis = codec._read_node[Int64, R, Self.D, Int](reader, budget, depth + 1)
            var value = codec._read_node[Self.D.Scalar, R, Self.D, Int](reader, budget, depth + 1)
            _grow(points, budget)
            points.append((axis, value^))
        var proof = codec._read_node[AdaptiveEvaluation[Self.D], R, Self.D, Int](reader, budget, depth + 1)
        var options = Optional[NumericOptions]()
        if reader.number(1, budget):
            options = codec._read_node[NumericOptions, R, Self.D, Int](reader, budget, depth + 1)
        var algebraic = Optional[_algnum._Polynomial[Self.D, Self.sparse]]()
        if reader.number(1, budget):
            algebraic = codec._read_node[_algnum._Polynomial[Self.D, Self.sparse], R, Self.D, Int](reader, budget, depth + 1)
        var algebraic_points = List[_Point[Self.D]]()
        for _ in range(reader.count(budget)):
            var axis = codec._read_node[Int64, R, Self.D, Int](reader, budget, depth + 1)
            var source = codec._read_node[Complex, R, Self.D, Int](reader, budget, depth + 1)
            var witness = codec._read_node[_algnum._Witness, R, Self.D, Int](reader, budget, depth + 1)
            _grow(algebraic_points, budget)
            algebraic_points.append(_Point[Self.D](axis, source^, witness^))
        _allocation[_Data[Self.C, Self.D, Self.numerical, Self.sparse]](1, budget)
        return Self(ArcPointer(_Data[Self.C, Self.D, Self.numerical, Self.sparse](p^, points^, proof^, options^, algebraic^, algebraic_points^)))


def _produce[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](p: _Core[C, numerical, sparse], points: List[Tuple[Int64, D.Scalar]], context: BallContext, mut budget: Budget) raises PolyError -> EnclosureEvaluation[D]:
    comptime assert D == BallDomain or D == ComplexBallDomain, "certified evaluation needs a ball domain"
    _num._context[D](rebind[D.Context](context), budget)
    comptime if numerical:
        comptime if sparse:
            return enclose_numeric_sparse[D, downcast[C, NumericDomain]](rebind[NumericSparsePoly[downcast[C, NumericDomain]]](p), points, rebind[D.Context](context), budget)
        else:
            return enclose_numeric[D, downcast[C, NumericDomain]](rebind[NumericDensePoly[downcast[C, NumericDomain]]](p), points[0][1], rebind[D.Context](context), budget)
    else:
        comptime if sparse:
            return enclose_sparse[D, downcast[C, CoefficientRing]](rebind[SparsePoly[downcast[C, CoefficientRing]]](p), points, rebind[D.Context](context), budget)
        else:
            return enclose_dense[D, downcast[C, CoefficientRing]](rebind[DensePoly[downcast[C, CoefficientRing]]](p), points[0][1], rebind[D.Context](context), budget)


def _options(options: NumericOptions) raises PolyError:
    if options.initial_precision < 2 or options.max_precision < options.initial_precision or options.absolute_bits < 0 or options.relative_bits < 0 or (options.absolute_bits == 0 and options.relative_bits == 0) or options.max_subdivisions < 0:
        raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical precision or accuracy options")


def _compute[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](var p: _Core[C, numerical, sparse], var points: List[Tuple[Int64, D.Scalar]], context: BallContext, options: Optional[NumericOptions], mut budget: Budget, var algebraic_points: List[_Point[D]] = List[_Point[D]]()) raises PolyError -> Enclosure[C, D, numerical, sparse]:
    comptime if not numerical and (C == Real or C == Complex):
        return rebind_var[Enclosure[C, D, numerical, sparse]](_algebraic_compute[C, D, sparse](rebind_var[_Core[C, False, sparse]](p^), points^, context, options, budget, algebraic_points^))
    var proof: AdaptiveEvaluation[D]
    if options:
        var request = options.value()
        _options(request)
        comptime if not numerical:
            comptime if sparse:
                proof = adaptive_enclose_sparse[D, downcast[C, CoefficientRing]](rebind[SparsePoly[downcast[C, CoefficientRing]]](p), points, request, budget)
            else:
                proof = adaptive_enclose_dense[D, downcast[C, CoefficientRing]](rebind[DensePoly[downcast[C, CoefficientRing]]](p), points[0][1], request, budget)
        else:
            var bits = request.initial_precision
            while True:
                budget.charge(1)
                var result = _produce[C, D, numerical, sparse](p, points, _context(bits, request.max_precision), budget)
                proof = AdaptiveEvaluation[D](NumericOutcome.SUCCESS, result^)
                if check_numeric_accuracy(proof, request, budget):
                    break
                if bits == request.max_precision:
                    proof.outcome = NumericOutcome.PRECISION_CEILING
                    break
                bits = request.max_precision if bits > request.max_precision // 2 else bits * 2
    else:
        proof = AdaptiveEvaluation[D](NumericOutcome.SUCCESS, _produce[C, D, numerical, sparse](p, points, context, budget))
    _allocation[_Data[C, D, numerical, sparse]](1, budget)
    return Enclosure[C, D, numerical, sparse](ArcPointer(_Data[C, D, numerical, sparse](p^, points^, proof^, options.copy(), None, algebraic_points^)))


def check[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](result: Enclosure[C, D, numerical, sparse], mut budget: Budget) raises PolyError -> Bool:
    """Verify every enclosure step and any retained accuracy claim.

    An inconclusive adaptive result can pass this check. Inspect accurate() or outcome() to learn whether the requested width was reached.

    Parameters:
        C: Inferred from the input types.
        D: The numerical coefficient domain.
        numerical: Inferred from the input types.
        sparse: Inferred from the input types.

    Args:
        result: The result and its retained evidence.
        budget: The shared work allowance.

    Returns:
        Whether the evidence is valid.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    if len(data.algebraic_points) > len(data.points):
        return False
    for i in range(len(data.algebraic_points)):
        budget.charge(1)
        ref point = data.algebraic_points[i]
        for j in range(i):
            budget.charge(1)
            if data.algebraic_points[j].axis == point.axis:
                return False
        var found = False
        for input in data.points:
            budget.charge(1)
            if input[0] == point.axis:
                if found or not _algnum._check[Complex, D](point.source, input[1], point.witness, budget):
                    return False
                found = True
        if not found:
            return False
    comptime if not sparse:
        if len(data.points) != 1 or data.points[0][0] != 0:
            return False
    var valid: Bool
    comptime if not numerical and (C == Real or C == Complex):
        if not data.algebraic or not _algnum._check_polynomial[downcast[C, CoefficientRing], D, sparse](rebind[_Core[C, False, sparse]](data.polynomial), data.algebraic.value(), budget):
            return False
        comptime if sparse:
            valid = check_enclose_numeric_sparse[D, D](rebind[NumericSparsePoly[D]](data.algebraic.value().value), data.points, data.proof.evaluation, budget)
        else:
            valid = check_enclose_numeric[D, D](rebind[NumericDensePoly[D]](data.algebraic.value().value), data.points[0][1], data.proof.evaluation, budget)
    elif numerical:
        if data.algebraic:
            return False
        comptime if sparse:
            valid = check_enclose_numeric_sparse[D, downcast[C, NumericDomain]](rebind[NumericSparsePoly[downcast[C, NumericDomain]]](data.polynomial), data.points, data.proof.evaluation, budget)
        else:
            valid = check_enclose_numeric[D, downcast[C, NumericDomain]](rebind[NumericDensePoly[downcast[C, NumericDomain]]](data.polynomial), data.points[0][1], data.proof.evaluation, budget)
    else:
        if data.algebraic:
            return False
        comptime if sparse:
            valid = check_enclose_sparse[D, downcast[C, CoefficientRing]](rebind[SparsePoly[downcast[C, CoefficientRing]]](data.polynomial), data.points, data.proof.evaluation, budget)
        else:
            valid = check_enclose_dense[D, downcast[C, CoefficientRing]](rebind[DensePoly[downcast[C, CoefficientRing]]](data.polynomial), data.points[0][1], data.proof.evaluation, budget)
    if not valid:
        return False
    return check_numeric_accuracy(data.proof, data.options.value(), budget) if data.options else True


def check[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool](result: Enclosure[C, D, numerical, sparse], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify every enclosure step and any retained accuracy claim.

    An inconclusive adaptive result can pass this check. Inspect accurate() or outcome() to learn whether the requested width was reached.

    Parameters:
        C: Inferred from the input types.
        D: The numerical coefficient domain.
        numerical: Inferred from the input types.
        sparse: Inferred from the input types.

    Args:
        result: The result and its retained evidence.
        max_steps: The complete call's work allowance.

    Returns:
        Whether the evidence is valid.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def _context(bits: Int, ceiling: Int) raises PolyError -> BallContext:
    try:
        return BallContext(bits, max_precision=ceiling)
    except error:
        raise PolyError(PolyErrorKind.DOMAIN, String(error))


from polynomials import EnclosureStream, EnclosureStreamChecker

comptime _RawStream[D: NumericDomain, checking: Bool] = EnclosureStreamChecker[D] if checking else EnclosureStream[D]


def _raw_stream[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool, checking: Bool](p: _Core[C, numerical, sparse], points: List[Tuple[Int64, D.Scalar]], context: BallContext, mut budget: Budget) raises PolyError -> _RawStream[D, checking]:
    comptime if numerical:
        comptime if sparse:
            comptime if checking:
                return rebind_var[_RawStream[D, checking]](EnclosureStreamChecker[D].from_numeric_sparse(rebind[NumericSparsePoly[downcast[C, NumericDomain]]](p), points, budget))
            else:
                return rebind_var[_RawStream[D, checking]](EnclosureStream[D].from_numeric_sparse(rebind[NumericSparsePoly[downcast[C, NumericDomain]]](p), points, rebind[D.Context](context), budget))
        else:
            comptime if checking:
                return rebind_var[_RawStream[D, checking]](EnclosureStreamChecker[D].from_numeric(rebind[NumericDensePoly[downcast[C, NumericDomain]]](p), points[0][1], budget))
            else:
                return rebind_var[_RawStream[D, checking]](EnclosureStream[D].from_numeric(rebind[NumericDensePoly[downcast[C, NumericDomain]]](p), points[0][1], rebind[D.Context](context), budget))
    else:
        comptime if sparse:
            comptime if checking:
                return rebind_var[_RawStream[D, checking]](EnclosureStreamChecker[D].from_sparse(rebind[SparsePoly[downcast[C, CoefficientRing]]](p), points, budget))
            else:
                return rebind_var[_RawStream[D, checking]](EnclosureStream[D].from_sparse(rebind[SparsePoly[downcast[C, CoefficientRing]]](p), points, rebind[D.Context](context), budget))
        else:
            comptime if checking:
                return rebind_var[_RawStream[D, checking]](EnclosureStreamChecker[D].from_dense(rebind[DensePoly[downcast[C, CoefficientRing]]](p), points[0][1], budget))
            else:
                return rebind_var[_RawStream[D, checking]](EnclosureStream[D].from_dense(rebind[DensePoly[downcast[C, CoefficientRing]]](p), points[0][1], rebind[D.Context](context), budget))

comptime _IsAlgebraic[T: AnyType] = T == Real or T == Complex
comptime _AlgebraicMember[previous: Bool, T: Movable] = previous or _IsAlgebraic[T]
comptime _AlgebraicTuple[T: type_of(Tuple)] = T.Ts.reduce[False, _AlgebraicMember]
comptime _Stream[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool, checking: Bool, algebraic_points: Bool] = AlgebraicEnclosureStream[C, D, numerical, sparse, checking] if _IsAlgebraic[C] or algebraic_points else _RawStream[D, checking]


@fieldwise_init
struct RefinedInput[D: NumericDomain](Copyable, Movable):
    """An outward numerical input and exact witnesses for its real and imaginary parts.

    Parameters:
        D: The inferred real or complex ball domain.
    """
    var value: Self.D.Scalar
    """The outward enclosure used by the arithmetic stream."""
    var real: RealAlgebraic
    """The refined real root, bound to the original input by the checker."""
    var imag: RealAlgebraic
    """The refined imaginary root; zero for real inputs."""


@fieldwise_init
struct EnclosureStep[D: NumericDomain](Copyable, Movable):
    """One coefficient refinement, point refinement, or arithmetic bound.

    Forward each step unchanged to the matching checker. Coefficients precede
    algebraic coordinates, which precede arithmetic. Only the active payload
    may be present.

    Parameters:
        D: The inferred real or complex ball domain.
    """
    var kind: Int
    """Zero for a coefficient, one for a coordinate, two for arithmetic."""
    var index: Int
    """Canonical coefficient index, original coordinate axis, or arithmetic step index."""
    var refinement: Optional[RefinedInput[Self.D]]
    """The input evidence, present only for coefficient and coordinate steps."""
    var bounds: Optional[NumericBox]
    """The exact arithmetic rectangle, present only for arithmetic steps."""


@fieldwise_init
struct _StreamPoint(Copyable, Movable):
    var axis: Int
    var source: Complex


def _stream_argument[D: NumericDomain, T: Copyable & Deinitable](source: T, axis: Int, context: D.Context, mut points: List[_StreamPoint], mut budget: Budget) raises PolyError -> D.Scalar:
    comptime if _IsAlgebraic[T]:
        var value: Complex
        comptime if T == Real:
            _charge_size(rebind[Real](source).size(), budget)
            value = Complex(rebind[Real](source).copy(), Real.zero(NoContext()))
        else:
            _charge_size(rebind[Complex](source).size(), budget)
            value = rebind[Complex](source).copy()
        var converted = _algnum._coarse[Complex, D](value, context, budget)
        _ = numeric_bounds[D](converted[0], budget)
        _grow(points, budget)
        points.append(_StreamPoint(axis, value^))
        return converted[0].copy()
    else:
        return _num._value[D](source, context, budget)


def _refined_input[C: CoefficientRing, D: NumericDomain](source: C, context: D.Context, mut budget: Budget) raises PolyError -> RefinedInput[D]:
    var pair = _algnum._enclose[C, D](source, context, budget)
    _ = numeric_bounds[D](pair[0], budget)
    return RefinedInput[D](pair[0].copy(), _copy_witness(pair[1].real, budget), _copy_witness(pair[1].imag, budget))


def _check_refinement[C: CoefficientRing, D: NumericDomain](source: C, input: RefinedInput[D], mut budget: Budget) raises PolyError -> Bool:
    try:
        return _algnum._check_parts[C, D](source, input.value, input.real, input.imag, budget)
    except error:
        if error.kind == PolyErrorKind.DOMAIN:
            return False
        raise error^


@fieldwise_init
struct AlgebraicEnclosureStream[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool, checking: Bool](Movable):
    """Refine algebraic inputs and check their evidence before streaming arithmetic.

    Construct this through dense.enclose_stream or sparse.enclose_stream, and
    construct its matching checker through enclose_checker. It keeps input
    storage and the live arithmetic state, but no history of emitted steps.

    Parameters:
        C: The original coefficient ring or numerical domain.
        D: The inferred real or complex ball domain.
        numerical: Whether the original polynomial has numerical coefficients.
        sparse: Whether the original polynomial uses sparse storage.
        checking: Whether this instance accepts evidence instead of producing it.
    """
    var _polynomial: _Core[Self.C, Self.numerical, Self.sparse]
    var _points: List[Tuple[Int64, Self.D.Scalar]]
    var _algebraic_points: List[_StreamPoint]
    var _converted: Optional[_algnum._Numeric[Self.D, Self.sparse]]
    var _engine: _RawStream[Self.D, Self.checking]
    var _context: BallContext
    var _coefficients: Int
    var _inputs: Int
    var _input_position: Int
    var _arithmetic_position: Int
    var _steps: Int
    var _ready: Bool
    var _failed: Bool

    def step_count(self) -> Int:
        """Return the complete stream length, including input evidence.

        Returns:
            The number of refinement and arithmetic events required for completion.
        """
        return self._steps

    def _prepare(mut self, mut budget: Budget) raises PolyError:
        if self._ready:
            return
        comptime if _IsAlgebraic[Self.C]:
            self._engine = _raw_stream[Self.D, Self.D, True, Self.sparse, Self.checking](self._converted.value(), self._points, self._context, budget)
        else:
            self._engine = _raw_stream[Self.C, Self.D, Self.numerical, Self.sparse, Self.checking](self._polynomial, self._points, self._context, budget)
        self._ready = True

    def _input(mut self, var value: Self.D.Scalar):
        # All fallible validation and copying precede this state transition.
        if self._input_position < self._coefficients:
            comptime if _IsAlgebraic[Self.C]:
                comptime if Self.sparse:
                    rebind[NumericSparsePoly[Self.D]](self._converted.value()).terms[self._input_position].coefficient = value^
                else:
                    rebind[NumericDensePoly[Self.D]](self._converted.value()).coefficients[self._input_position] = value^
        else:
            var axis = self._algebraic_points[self._input_position - self._coefficients].axis
            self._points[axis][1] = value^
        self._input_position += 1

    def next(mut self, mut budget: Budget) raises PolyError -> Optional[EnclosureStep[Self.D]]:
        """Refine one input or produce the next arithmetic bound.

        Args:
            budget: The shared allowance for refinement, storage, and arithmetic.

        Returns:
            The next event, or None after all events have been emitted.

        Raises:
            PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on resource
                exhaustion. A failed operation does not consume its event.
        """
        comptime assert not Self.checking, "next requires an enclosure producer"
        budget.charge(1)
        if self._input_position < self._inputs:
            var input: RefinedInput[Self.D]
            var kind = 1
            var index: Int
            if self._input_position < self._coefficients:
                kind = 0
                index = self._input_position
                comptime if _IsAlgebraic[Self.C]:
                    comptime if Self.sparse:
                        input = _refined_input[downcast[Self.C, CoefficientRing], Self.D](rebind[SparsePoly[downcast[Self.C, CoefficientRing]]](self._polynomial).terms[index].coefficient, rebind[Self.D.Context](self._context), budget)
                    else:
                        input = _refined_input[downcast[Self.C, CoefficientRing], Self.D](rebind[DensePoly[downcast[Self.C, CoefficientRing]]](self._polynomial).coefficients[index], rebind[Self.D.Context](self._context), budget)
                else:
                    raise PolyError(PolyErrorKind.DOMAIN, "unexpected algebraic coefficient")
            else:
                ref point = self._algebraic_points[self._input_position - self._coefficients]
                index = point.axis
                input = _refined_input[Complex, Self.D](point.source, rebind[Self.D.Context](self._context), budget)
            # Copy just the scalar into prepared arithmetic input storage.
            _ = numeric_bounds[Self.D](input.value, budget)
            var stored = input.value.copy()
            if kind == 0:
                comptime if _IsAlgebraic[Self.C]:
                    comptime if Self.sparse:
                        rebind[NumericSparsePoly[Self.D]](self._converted.value()).terms[index].coefficient = stored^
                    else:
                        rebind[NumericDensePoly[Self.D]](self._converted.value()).coefficients[index] = stored^
            else:
                self._points[index][1] = stored^
            self._input_position += 1
            return EnclosureStep[Self.D](kind, index, input^, None)
        self._prepare(budget)
        var bounds = rebind[EnclosureStream[Self.D]](self._engine).next(budget)
        if not bounds:
            return None
        var index = self._arithmetic_position
        self._arithmetic_position += 1
        return EnclosureStep[Self.D](2, index, None, bounds.take())

    def accept(mut self, var step: EnclosureStep[Self.D], mut budget: Budget) raises PolyError -> Bool:
        """Check the next root witness or arithmetic bound independently.

        Root identity, interval containment, and outward bounds are checked
        without repeating the producer's root refinement.

        Args:
            step: The next owned event in the producer's order.
            budget: The shared checking allowance.

        Returns:
            True on acceptance. False permanently invalidates this checker.

        Raises:
            PolyError: LIMIT or CAPACITY on exhaustion; the event is not consumed.
        """
        comptime assert Self.checking, "accept requires an enclosure checker"
        budget.charge(1)
        if self._failed:
            return False
        if self._input_position < self._inputs:
            var kind = 0 if self._input_position < self._coefficients else 1
            var index = self._input_position if kind == 0 else self._algebraic_points[self._input_position - self._coefficients].axis
            if step.kind != kind or step.index != index or not step.refinement or step.bounds:
                self._failed = True
                return False
            ref input = step.refinement.value()
            var valid = False
            if kind == 0:
                comptime if _IsAlgebraic[Self.C]:
                    comptime if Self.sparse:
                        valid = _check_refinement[downcast[Self.C, CoefficientRing], Self.D](rebind[SparsePoly[downcast[Self.C, CoefficientRing]]](self._polynomial).terms[index].coefficient, input, budget)
                    else:
                        valid = _check_refinement[downcast[Self.C, CoefficientRing], Self.D](rebind[DensePoly[downcast[Self.C, CoefficientRing]]](self._polynomial).coefficients[index], input, budget)
            else:
                valid = _check_refinement[Complex, Self.D](self._algebraic_points[self._input_position - self._coefficients].source, input, budget)
            if not valid:
                self._failed = True
                return False
            _ = numeric_bounds[Self.D](input.value, budget)
            self._input(input.value.copy())
            return True
        if step.kind != 2 or step.index != self._arithmetic_position or step.refinement or not step.bounds:
            self._failed = True
            return False
        self._prepare(budget)
        if not rebind[EnclosureStreamChecker[Self.D]](self._engine).accept(step.bounds.take(), budget):
            self._failed = True
            return False
        self._arithmetic_position += 1
        return True

    def finish(self, mut budget: Budget) raises PolyError -> StreamedEnclosure[Self.D]:
        """Return a receipt after every refinement and arithmetic event.

        Args:
            budget: The allowance for copying the final enclosure.

        Returns:
            A scalar receipt whose steps include refinement events.

        Raises:
            PolyError: DOMAIN for an incomplete stream, LIMIT or CAPACITY on exhaustion.
        """
        comptime assert not Self.checking, "producer finish requires no receipt"
        if not self._ready or self._input_position != self._inputs:
            raise PolyError(PolyErrorKind.DOMAIN, "incomplete algebraic enclosure stream")
        var receipt = rebind[EnclosureStream[Self.D]](self._engine).finish(budget)
        receipt.steps = self._steps
        return receipt^

    def finish(self, receipt: StreamedEnclosure[Self.D], mut budget: Budget) raises PolyError -> Bool:
        """Verify completion and the final scalar against accepted input evidence.

        Args:
            receipt: The producer's final scalar and total event count.
            budget: The shared allowance for final verification.

        Returns:
            Whether all expected events and the receipt have been checked.

        Raises:
            PolyError: LIMIT or CAPACITY on exhaustion.
        """
        comptime assert Self.checking, "checker finish requires a receipt"
        budget.charge(1)
        if self._failed or not self._ready or self._input_position != self._inputs or receipt.steps != self._steps or receipt.precision != self._context.precision().value() or len(receipt._value) != 1:
            return False
        var value = receipt.value(budget)
        _allocation[Self.D.Scalar](1, budget)
        var arithmetic = StreamedEnclosure[Self.D](receipt.precision, receipt.steps - self._inputs, [value^])
        return rebind[EnclosureStreamChecker[Self.D]](self._engine).finish(arithmetic, budget)


def _stream[C: AnyType, D: NumericDomain, numerical: Bool, sparse: Bool, checking: Bool, algebraic_points: Bool](var p: _Core[C, numerical, sparse], var points: List[Tuple[Int64, D.Scalar]], context: BallContext, mut budget: Budget, var roots: List[_StreamPoint]) raises PolyError -> _Stream[C, D, numerical, sparse, checking, algebraic_points]:
    comptime if _IsAlgebraic[C] or algebraic_points:
        var converted = Optional[_algnum._Numeric[D, sparse]]()
        var count = 0
        var engine: _RawStream[D, checking]
        comptime if _IsAlgebraic[C]:
            # Validate sources and prepare outward placeholders. No root is
            # refined until its event is requested, and no witness list survives.
            var coarse = _algnum._convert[downcast[C, CoefficientRing], D, sparse](rebind[_Core[C, False, sparse]](p), rebind[D.Context](context), budget, False)
            comptime if sparse:
                var value = NumericSparsePoly[D](MonomialOrder())
                swap(value, rebind[NumericSparsePoly[D]](coarse.value))
                converted = rebind_var[_algnum._Numeric[D, sparse]](value^)
            else:
                var value = NumericDensePoly[D]()
                swap(value, rebind[NumericDensePoly[D]](coarse.value))
                converted = rebind_var[_algnum._Numeric[D, sparse]](value^)
            comptime if sparse:
                count = len(rebind[SparsePoly[downcast[C, CoefficientRing]]](p).terms)
            else:
                count = len(rebind[DensePoly[downcast[C, CoefficientRing]]](p).coefficients)
            engine = _raw_stream[D, D, True, sparse, checking](converted.value(), points, context, budget)
        else:
            engine = _raw_stream[C, D, numerical, sparse, checking](p, points, context, budget)
        var inputs = _sum(count, len(roots))
        var arithmetic_steps: Int
        comptime if checking:
            arithmetic_steps = rebind[EnclosureStreamChecker[D]](engine).step_count()
        else:
            arithmetic_steps = rebind[EnclosureStream[D]](engine).step_count()
        var steps = _sum(inputs, arithmetic_steps)
        return rebind_var[_Stream[C, D, numerical, sparse, checking, algebraic_points]](AlgebraicEnclosureStream[C, D, numerical, sparse, checking](p^, points^, roots^, converted^, engine^, context, count, inputs, 0, 0, steps, inputs == 0, False))
    else:
        return rebind_var[_Stream[C, D, numerical, sparse, checking, algebraic_points]](_raw_stream[C, D, numerical, sparse, checking](p, points, context, budget))


def _algebraic_compute[C: AnyType, D: NumericDomain, sparse: Bool](var p: _Core[C, False, sparse], var points: List[Tuple[Int64, D.Scalar]], context: BallContext, options: Optional[NumericOptions], mut budget: Budget, var algebraic_points: List[_Point[D]] = List[_Point[D]]()) raises PolyError -> Enclosure[C, D, False, sparse]:
    if options:
        _options(options.value())
    _num._context[D](rebind[D.Context](context), budget)
    var bits = options.value().initial_precision if options else context.precision().value()
    var ceiling = options.value().max_precision if options else bits
    var proof: AdaptiveEvaluation[D]
    var retained: Optional[_algnum._Polynomial[D, sparse]]
    while True:
        budget.charge(1)
        var working = _context(bits, ceiling)
        var converted = _algnum._convert[downcast[C, CoefficientRing], D, sparse](p, rebind[D.Context](working), budget)
        var result = _produce[D, D, True, sparse](converted.value, points, working, budget)
        proof = AdaptiveEvaluation[D](NumericOutcome.SUCCESS, result^)
        retained = converted^
        if not options or check_numeric_accuracy(proof, options.value(), budget):
            break
        if bits == ceiling:
            proof.outcome = NumericOutcome.PRECISION_CEILING
            break
        bits = ceiling if bits > ceiling // 2 else bits * 2
    _allocation[_Data[C, D, False, sparse]](1, budget)
    return Enclosure[C, D, False, sparse](ArcPointer(_Data[C, D, False, sparse](p^, points^, proof^, options.copy(), retained^, algebraic_points^)))
