"""Algebraic coefficients bound to independently checked numerical root covers."""
from std.memory import ArcPointer
from polynomials import Budget, PolyError, PolyErrorKind, NumericOptions, NumericOutcome, BallDomain, ComplexBallDomain, RealRootResult, ComplexRootResult, RealRootSink, ComplexRootSink, RealRootStreamChecker, ComplexRootStreamChecker, RootStreamReceipt, RealRootCell, ComplexRootRegion, SturmSequence, PolynomialTextSink, PolynomialTextSource
from polynomials import real_roots_numeric, complex_roots_numeric, real_roots_numeric_to, complex_roots_numeric_to, check_real_roots_numeric, check_complex_roots_numeric
from . import dense, _algebraic_numeric as an, evaluation
from ._common import _Ring, _allocation
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer

comptime _Domain[z: Bool] = ComplexBallDomain if z else BallDomain
comptime _Proof[z: Bool] = ComplexRootResult if z else RealRootResult


def _convert[k: Int, z: Bool](source: dense.Array[k], options: NumericOptions, mut budget: Budget) raises PolyError -> an._Polynomial[_Domain[z], False]:
    options._validate()
    comptime assert k == 4 or (k == 5 and z), "algebraic real covers require AA coefficients"
    # A fixed, reproducible coefficient family lets producers and incremental
    # checkers share a mathematical claim without buffering or replaying streams.
    var result = an._convert[_Ring[k], _Domain[z], False](source._core[], rebind[_Domain[z].Context](evaluation._context(options.max_precision, options.max_precision)), budget)
    if not an._check_polynomial[_Ring[k], _Domain[z], False](source._core[], result, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic coefficient enclosure")
    return result^


@fieldwise_init
struct AlgebraicCover[k: Int, z: Bool](_Archive, Copyable, Movable):
    """A numerical cover bound to exact algebraic coefficients and conversion witnesses.

    Parameters:
        k: The inferred exact algebraic coefficient kind.
        z: Whether the cover contains complex rectangles.
    """
    var _source: dense.Array[Self.k]
    var _converted: ArcPointer[an._Polynomial[_Domain[Self.z], False]]
    var _options: NumericOptions
    var _proof: ArcPointer[_Proof[Self.z]]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._proof[]))] _Proof[Self.z]:
        """Borrow the cover certificate, including unresolved regions.

        Returns:
            The immutable numerical root evidence.
        """
        return self._proof[]

    def options(self) -> NumericOptions:
        """Copy the requested numerical guarantees.

        Returns:
            The retained precision and root-cover options.
        """
        return self._options

    def complete(self) -> Bool:
        """Report whether all roots are accounted for.

        Returns:
            The retained claim.
        """
        comptime if Self.z:
            return rebind[ComplexRootResult](self._proof[]).complete
        else:
            return rebind[RealRootResult](self._proof[]).complete

    def accurate(self) -> Bool:
        """Report whether regions meet the requested width.

        Returns:
            The retained claim.
        """
        comptime if Self.z:
            return rebind[ComplexRootResult](self._proof[]).accurate
        else:
            return rebind[RealRootResult](self._proof[]).accurate

    def isolated(self) -> Bool:
        """Report whether every counted region has count one.

        Returns:
            The retained claim.
        """
        comptime if Self.z:
            return rebind[ComplexRootResult](self._proof[]).isolated
        else:
            return rebind[RealRootResult](self._proof[]).isolated

    def outcome(self) -> NumericOutcome:
        """Report whether the requested guarantees succeeded.

        Returns:
            The retained numerical outcome.
        """
        comptime if Self.z:
            return rebind[ComplexRootResult](self._proof[]).outcome
        else:
            return rebind[RealRootResult](self._proof[]).outcome

    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return _check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io as codec
        codec._depth(depth, budget)
        codec._write_node[dense.Array[Self.k], _Ring[Self.k], _Domain[Self.z], Int](self._source, writer, sink, budget, depth + 1)
        codec._write_node[an._Polynomial[_Domain[Self.z], False], _Ring[Self.k], _Domain[Self.z], Int](self._converted[], writer, sink, budget, depth + 1)
        codec._write_node[NumericOptions, _Ring[Self.k], _Domain[Self.z], Int](self._options, writer, sink, budget, depth + 1)
        codec._write_node[_Proof[Self.z], _Ring[Self.k], _Domain[Self.z], Int](self._proof[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io as codec
        codec._depth(depth, budget)
        var p = codec._read_node[dense.Array[Self.k], _Ring[Self.k], _Domain[Self.z], Int](reader, budget, depth + 1)
        var converted = codec._read_node[an._Polynomial[_Domain[Self.z], False], _Ring[Self.k], _Domain[Self.z], Int](reader, budget, depth + 1)
        var options = codec._read_node[NumericOptions, _Ring[Self.k], _Domain[Self.z], Int](reader, budget, depth + 1)
        var proof = codec._read_node[_Proof[Self.z], _Ring[Self.k], _Domain[Self.z], Int](reader, budget, depth + 1)
        _allocation[an._Polynomial[_Domain[Self.z], False]](1, budget)
        _allocation[_Proof[Self.z]](1, budget)
        return Self(p^, ArcPointer(converted^), options, ArcPointer(proof^))


def _make[k: Int, z: Bool](source: dense.Array[k], options: NumericOptions, mut budget: Budget, complete_only: Bool) raises PolyError -> AlgebraicCover[k, z]:
    var converted = _convert[k, z](source, options, budget)
    var proof: _Proof[z]
    comptime if z:
        proof = rebind_var[_Proof[z]](complex_roots_numeric(converted.value, options, budget, complete_only=complete_only))
    else:
        proof = rebind_var[_Proof[z]](real_roots_numeric(converted.value, options, budget, complete_only=complete_only))
    _allocation[an._Polynomial[_Domain[z], False]](1, budget)
    _allocation[_Proof[z]](1, budget)
    return AlgebraicCover[k, z](source.copy(), ArcPointer(converted^), options, ArcPointer(proof^))


def _check[k: Int, z: Bool](result: AlgebraicCover[k, z], mut budget: Budget) raises PolyError -> Bool:
    if not an._check_polynomial[_Ring[k], _Domain[z], False](result._source._core[], result._converted[], budget):
        return False
    comptime if z:
        return check_complex_roots_numeric(result._converted[].value, rebind[ComplexRootResult](result._proof[]), result._options, budget)
    else:
        return check_real_roots_numeric(result._converted[].value, rebind[RealRootResult](result._proof[]), result._options, budget)


@fieldwise_init
struct AlgebraicRealRootChecker[k: Int](RealRootSink, Movable):
    """Check a numerical root stream against exact algebraic coefficients.

    Parameters:
        k: The inferred exact algebraic coefficient kind.
    """
    var _source: dense.Array[Self.k]
    var _converted: an._Polynomial[_Domain[False], False]
    var _checker: RealRootStreamChecker

    def begin(mut self, sequence: Optional[SturmSequence], mut budget: Budget) raises PolyError:
        """Accept the shared real-root stream header.

        Args:
            sequence: The optional rational Sturm witness.
            budget: The shared checking allowance.

        Raises:
            PolyError: DOMAIN on invalid headers, or LIMIT/CAPACITY.
        """
        self._checker.begin(sequence, budget)

    def accept(mut self, var witness: RealRootCell, mut budget: Budget) raises PolyError -> Bool:
        """Check and consume the next local root witness.

        Args:
            witness: The next owned numerical root certificate.
            budget: The shared checking allowance.

        Returns:
            Whether the witness is valid; rejection invalidates the checker.

        Raises:
            PolyError: DOMAIN on invalid evidence, or LIMIT/CAPACITY on resource limits.
        """
        return self._checker.accept(witness^, budget)

    def write(mut self, var witness: RealRootCell, mut budget: Budget) raises PolyError:
        """Check and consume the next local root witness.

        Args:
            witness: The next owned numerical root certificate.
            budget: The shared checking allowance.

        Raises:
            PolyError: DOMAIN on invalid evidence, or LIMIT/CAPACITY on resource limits.
        """
        self._checker.write(witness^, budget)

    def finish(self, receipt: RootStreamReceipt, mut budget: Budget) raises PolyError -> Bool:
        """Verify conversion, coverage, and the producer's completion claims.

        Args:
            receipt: The final receipt after all local witnesses.
            budget: The shared final checking allowance.

        Returns:
            Whether the complete stream certifies the retained algebraic input.

        Raises:
            PolyError: DOMAIN on invalid input, or LIMIT/CAPACITY on resource limits.
        """
        return an._check_polynomial[_Ring[Self.k], _Domain[False], False](self._source._core[], self._converted, budget) and self._checker.finish(receipt, budget)


def _real_checker[k: Int](source: dense.Array[k], options: NumericOptions, mut budget: Budget) raises PolyError -> AlgebraicRealRootChecker[k]:
    var converted = _convert[k, False](source, options, budget)
    var checker = RealRootStreamChecker.from_numeric(converted.value, options, budget)
    return AlgebraicRealRootChecker[k](source.copy(), converted^, checker^)


def _real_to[k: Int, S: RealRootSink](source: dense.Array[k], options: NumericOptions, mut sink: S, mut budget: Budget, complete_only: Bool) raises PolyError -> RootStreamReceipt:
    var converted = _convert[k, False](source, options, budget)
    return real_roots_numeric_to(converted.value, options, sink, budget, complete_only=complete_only)


@fieldwise_init
struct AlgebraicComplexRootChecker[k: Int](ComplexRootSink, Movable):
    """Check a numerical root stream against exact algebraic coefficients.

    Parameters:
        k: The inferred exact algebraic coefficient kind.
    """
    var _source: dense.Array[Self.k]
    var _converted: an._Polynomial[_Domain[True], False]
    var _checker: ComplexRootStreamChecker

    def accept(mut self, var witness: ComplexRootRegion, mut budget: Budget) raises PolyError -> Bool:
        """Check and consume the next local root witness.

        Args:
            witness: The next owned numerical root certificate.
            budget: The shared checking allowance.

        Returns:
            Whether the witness is valid; rejection invalidates the checker.

        Raises:
            PolyError: DOMAIN on invalid evidence, or LIMIT/CAPACITY on resource limits.
        """
        return self._checker.accept(witness^, budget)

    def write(mut self, var witness: ComplexRootRegion, mut budget: Budget) raises PolyError:
        """Check and consume the next local root witness.

        Args:
            witness: The next owned numerical root certificate.
            budget: The shared checking allowance.

        Raises:
            PolyError: DOMAIN on invalid evidence, or LIMIT/CAPACITY on resource limits.
        """
        self._checker.write(witness^, budget)

    def finish(self, receipt: RootStreamReceipt, mut budget: Budget) raises PolyError -> Bool:
        """Verify conversion, coverage, and the producer's completion claims.

        Args:
            receipt: The final receipt after all local witnesses.
            budget: The shared final checking allowance.

        Returns:
            Whether the complete stream certifies the retained algebraic input.

        Raises:
            PolyError: DOMAIN on invalid input, or LIMIT/CAPACITY on resource limits.
        """
        return an._check_polynomial[_Ring[Self.k], _Domain[True], False](self._source._core[], self._converted, budget) and self._checker.finish(receipt, budget)


def _complex_checker[k: Int](source: dense.Array[k], options: NumericOptions, mut budget: Budget) raises PolyError -> AlgebraicComplexRootChecker[k]:
    var converted = _convert[k, True](source, options, budget)
    var checker = ComplexRootStreamChecker.from_numeric(converted.value, options, budget)
    return AlgebraicComplexRootChecker[k](source.copy(), converted^, checker^)


def _complex_to[k: Int, S: ComplexRootSink](source: dense.Array[k], options: NumericOptions, mut sink: S, mut budget: Budget, complete_only: Bool) raises PolyError -> RootStreamReceipt:
    var converted = _convert[k, True](source, options, budget)
    return complex_roots_numeric_to(converted.value, options, sink, budget, complete_only=complete_only)
