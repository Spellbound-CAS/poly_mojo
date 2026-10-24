"""Root results, independent verification, and functional cover refinement."""
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink, BallDomain, QQ
from ._common import _grow
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from apn_mojo import Rational
from polynomials import Budget, PolyError, NoContext, DensePoly, RealAlgebraic, AlgebraicRootsCertificate, algebraic_roots, check_algebraic_roots, real_rational
from polynomials import DescartesCountCertificate, descartes_count as _descartes_count, check_descartes_count
from polynomials import NumericOptions, NumericOutcome, RealRootResult, ComplexRootResult, DescartesRootResult
from polynomials import check_real_roots_exact, check_complex_roots_exact, refine_real_roots_exact, refine_complex_roots_exact
from ._common import _Ring, _allocation, _charge_size, Values
from .algebraic import Real, _plain, _normal, _copy_witness, _same_witness
from polynomials import check_real_algebraic


@fieldwise_init
struct _ExactData(Copyable, Movable):
    var coefficients: List[RealAlgebraic]
    var proof: AlgebraicRootsCertificate


@fieldwise_init
struct ExactRoots( _Archive, Copyable, Movable, Iterable, Sized):
    """All distinct real roots in increasing order, with completeness evidence.

    Indexing and iteration return exact Real values. A repeated root appears
    once. Copies share both the values and the immutable certificate.
    """
    var _data: ArcPointer[_ExactData]
    var _values: Values[Real, 4]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _ExactIterator[ImmOrigin(iterable_origin)]
    """A borrowing iterator over exact root values.

    Parameters:
        iterable_mut: Whether the source origin permits mutation.
        iterable_origin: The source lifetime.
    """

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over exact values in increasing order.

        Returns:
            An iterator sharing the root witnesses.
        """
        return _ExactIterator[ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def __len__(self) -> Int:
        """Return the number of distinct real roots.

        Returns:
            The exact root count, excluding multiplicity.
        """
        return len(self._values)

    def __getitem__(self, index: Int) -> Real:
        """Return the exact root at a valid zero-based index.

        Args:
            index: The root's position in increasing order.

        Returns:
            An exact value sharing its stored witness.
        """
        return self._values[index]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] AlgebraicRootsCertificate:
        """Borrow the complete root certificate.

        Returns:
            Immutable completeness, selection, and evaluation evidence.
        """
        return self._data[].proof



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_ExactData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)
        writer.line(String(len(self._values)), sink, budget)
        for field655 in self._values:
            _result_io._write_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](field655, writer, sink, budget, depth + 1)
        writer.line(String(self._values.modulus()), sink, budget)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field657 = _result_io._read_node[_ExactData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_ExactData](1, budget)
        var field656 = ArcPointer(field657^)
        var field660 = reader.count(budget)
        var field659 = List[Real]()
        for _ in range(field660):
            var field661 = _result_io._read_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            _grow(field659, budget)
            field659.append(field661^)
        var field658 = _common.Values[Real, 4](field659^, reader.number(UInt64.MAX, budget), budget)
        return Self(field656^, field658^)

@fieldwise_init
struct _ExactIterator[origin: ImmOrigin](Iterator):
    comptime Element = Real
    var source: Pointer[ExactRoots, Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Real:
        if self.position == len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position]
        self.position += 1
        return result^


def _exact(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> ExactRoots:
    _allocation[RealAlgebraic](len(p.coefficients), budget)
    var coefficients = List[RealAlgebraic](capacity=len(p.coefficients))
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        coefficients.append(real_rational(c, budget))
    var proof = algebraic_roots(coefficients, budget)
    _allocation[Real](len(proof.roots), budget)
    var values = List[Real](capacity=len(proof.roots))
    # Candidates are every real root of the annihilator on its monic
    # irreducible factors, in increasing order: a root's index among its
    # polynomial's real roots is the number of earlier candidates sharing it.
    var seen = List[Int]()
    var counts = List[Int]()
    var k = 0
    for j in range(len(proof.candidates)):
        budget.charge(1)
        ref p = proof.candidates[j].polynomial
        var slot = -1
        for t in range(len(seen)):
            budget.charge(1)
            ref q = proof.candidates[seen[t]].polynomial
            if q.degree() == p.degree() and q == p:
                slot = t
                break
        if slot < 0:
            slot = len(seen)
            seen.append(j)
            counts.append(0)
        if k < len(proof.accepted) and proof.accepted[k] == j:
            values.append(_normal(_copy_witness(proof.roots[k], budget), counts[slot], budget))
            k += 1
        counts[slot] += 1
    _allocation[_ExactData](1, budget)
    return ExactRoots(ArcPointer(_ExactData(coefficients^, proof^)), Values[Real, 4](values^, 0, budget))


@fieldwise_init
struct _CountData(Copyable, Movable):
    var polynomial: DensePoly[QQ]
    var proof: DescartesCountCertificate


@fieldwise_init
struct RootCount( _Archive, Copyable, Movable):
    """A distinct-root count and its independently checkable Descartes evidence."""
    var _data: ArcPointer[_CountData]

    def count(self) -> Int:
        """Return the number of distinct roots in the requested open interval.

        Returns:
            The certified distinct-root count.
        """
        return self._data[].proof.count

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] DescartesCountCertificate:
        """Borrow the root-count certificate.

        Returns:
            Immutable endpoints, count, and the complete Descartes isolation behind it.
        """
        return self._data[].proof



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_CountData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field663 = _result_io._read_node[_CountData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_CountData](1, budget)
        var field662 = ArcPointer(field663^)
        return Self(field662^)

def _count(var p: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> RootCount:
    var proof = _descartes_count(p, low, high, budget)
    _allocation[_CountData](1, budget)
    return RootCount(ArcPointer(_CountData(p^, proof^)))


def check(result: ExactRoots, mut budget: Budget) raises PolyError -> Bool:
    """Verify that the result contains exactly all distinct real roots.

    Args:
        result: An exact-root result.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Whether the complete root certificate is valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    ref proof = result._data[].proof
    if len(result) != len(proof.roots) or not check_algebraic_roots(result._data[].coefficients, proof, budget):
        return False
    for i in range(len(result)):
        budget.charge(1)
        var value = result[i]
        if not _same_witness(value.witness(), proof.roots[i], budget) or not check_real_algebraic(value.witness(), budget):
            return False
    return True


def check(result: ExactRoots, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify that the result contains exactly all distinct real roots.

    Args:
        result: An exact-root result.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the complete root certificate is valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def check(result: RootCount, mut budget: Budget) raises PolyError -> Bool:
    """Verify an interval root count.

    Args:
        result: A retained root count.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Whether the Descartes evidence establishes the count.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return check_descartes_count(result._data[].polynomial, result._data[].proof, budget)


def check(result: RootCount, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify an interval root count.

    Args:
        result: A retained root count.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the Descartes evidence establishes the count.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


comptime _REAL_COVER_FORMAT = "format descartes-1"
"""The archive marker of a real cover's evidence format."""


@fieldwise_init
struct _RealCoverData[k: Int](Copyable, Movable):
    var polynomial: DensePoly[_Ring[Self.k]]
    var options: NumericOptions
    var proof: DescartesRootResult


@fieldwise_init
struct RealCover[k: Int]( _Archive, Copyable, Movable):
    """A real root cover retaining its polynomial, options, and certificate.

    Parameters:
        k: The inferred exact coefficient kind.
    """
    var _data: ArcPointer[_RealCoverData[Self.k]]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] DescartesRootResult:
        """Borrow the full Descartes evidence, including unresolved leaves.

        Returns:
            The immutable core real-root result and certificate.
        """
        return self._data[].proof

    def options(self) -> NumericOptions:
        """Return the precision and accuracy settings that define the claim.

        Returns:
            An independent copy of the retained options.
        """
        return self._data[].options

    def complete(self) -> Bool:
        """Report whether all roots are accounted for.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.complete

    def accurate(self) -> Bool:
        """Report whether counted regions meet the requested width.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.accurate

    def isolated(self) -> Bool:
        """Report whether each counted region has count one.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.isolated

    def outcome(self) -> NumericOutcome:
        """Report whether the requested claims succeeded or why they remain unresolved.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.outcome



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        # The evidence format: earlier archives held a Sturm cell cover.
        writer.line(_REAL_COVER_FORMAT, sink, budget)
        _result_io._write_node[_RealCoverData[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        if reader.line(budget) != _REAL_COVER_FORMAT:
            raise PolyError(PolyErrorKind.ENCODING, "real cover archive predates Descartes evidence; recompute the cover")
        var field669 = _result_io._read_node[_RealCoverData[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_RealCoverData[Self.k]](1, budget)
        var field668 = ArcPointer(field669^)
        return Self(field668^)

def check[k: Int](result: RealCover[k], mut budget: Budget) raises PolyError -> Bool:
    """Verify a real root cover against its retained input and options.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The real root cover.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Whether the cover is complete and its reported claims are supported. Incomplete covers return False.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return check_real_roots_exact(result._data[].polynomial, result._data[].proof, result._data[].options, budget)


def check[k: Int](result: RealCover[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify a real root cover against its retained input and options.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The real root cover.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the cover is complete and its reported claims are supported. Incomplete covers return False.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def refine[k: Int](result: RealCover[k], options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> RealCover[k]:
    """Refine a real cover while preserving the previous result.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The previous cover, including its original options.
        options: The new accuracy, precision, and subdivision settings.
        budget: The shared work allowance for preparation and computation.
        complete_only: Raise unless every requested property succeeds.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    # A cover's evidence was made here or checked when restored, so the
    # brackets are tightened without checking it again.
    var proof = refine_real_roots_exact(result._data[].polynomial, result._data[].proof, result._data[].options, options, budget, complete_only=complete_only, trusted=True)
    from ._interop import _dense_copy
    var p = _dense_copy(result._data[].polynomial, 0, budget)
    _allocation[_RealCoverData[k]](1, budget)
    return RealCover[k](ArcPointer(_RealCoverData[k](p^, options, proof^)))


def refine[k: Int](result: RealCover[k], options: NumericOptions, *, max_steps: UInt64 = UInt64.MAX, complete_only: Bool = False) raises PolyError -> RealCover[k]:
    """Refine a real cover while preserving the previous result.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The previous cover, including its original options.
        options: The new accuracy, precision, and subdivision settings.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        complete_only: Raise unless every requested property succeeds.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return refine(result, options, budget, complete_only=complete_only)


@fieldwise_init
struct _ComplexCoverData[k: Int](Copyable, Movable):
    var polynomial: DensePoly[_Ring[Self.k]]
    var options: NumericOptions
    var proof: ComplexRootResult


@fieldwise_init
struct ComplexCover[k: Int]( _Archive, Copyable, Movable):
    """A complex root cover retaining its polynomial, options, and certificate.

    Parameters:
        k: The inferred exact coefficient kind.
    """
    var _data: ArcPointer[_ComplexCoverData[Self.k]]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] ComplexRootResult:
        """Borrow the full cover evidence, including unresolved regions.

        Returns:
            The immutable core root-cover certificate.
        """
        return self._data[].proof

    def options(self) -> NumericOptions:
        """Return the precision and accuracy settings that define the claim.

        Returns:
            An independent copy of the retained options.
        """
        return self._data[].options

    def complete(self) -> Bool:
        """Report whether all roots are accounted for.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.complete

    def accurate(self) -> Bool:
        """Report whether counted regions meet the requested width.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.accurate

    def isolated(self) -> Bool:
        """Report whether each counted region has count one.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.isolated

    def outcome(self) -> NumericOutcome:
        """Report whether the requested claims succeeded or why they remain unresolved.

        Returns:
            The retained claim or numerical outcome.
        """
        return self._data[].proof.outcome



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_ComplexCoverData[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field671 = _result_io._read_node[_ComplexCoverData[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_ComplexCoverData[Self.k]](1, budget)
        var field670 = ArcPointer(field671^)
        return Self(field670^)

def check[k: Int](result: ComplexCover[k], mut budget: Budget) raises PolyError -> Bool:
    """Verify a complex root cover against its retained input and options.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The complex root cover.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Whether the cover is complete and its reported claims are supported. Incomplete covers return False.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return check_complex_roots_exact(result._data[].polynomial, result._data[].proof, result._data[].options, budget)


def check[k: Int](result: ComplexCover[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify a complex root cover against its retained input and options.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The complex root cover.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the cover is complete and its reported claims are supported. Incomplete covers return False.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def refine[k: Int](result: ComplexCover[k], options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> ComplexCover[k]:
    """Refine a complex cover while preserving the previous result.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The previous cover, including its original options.
        options: The new accuracy, precision, and subdivision settings.
        budget: The shared work allowance for preparation and computation.
        complete_only: Raise unless every requested property succeeds.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var proof = refine_complex_roots_exact(result._data[].polynomial, result._data[].proof, result._data[].options, options, budget, complete_only=complete_only)
    from ._interop import _dense_copy
    var p = _dense_copy(result._data[].polynomial, 0, budget)
    _allocation[_ComplexCoverData[k]](1, budget)
    return ComplexCover[k](ArcPointer(_ComplexCoverData[k](p^, options, proof^)))


def refine[k: Int](result: ComplexCover[k], options: NumericOptions, *, max_steps: UInt64 = UInt64.MAX, complete_only: Bool = False) raises PolyError -> ComplexCover[k]:
    """Refine a complex cover while preserving the previous result.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: The previous cover, including its original options.
        options: The new accuracy, precision, and subdivision settings.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        complete_only: Raise unless every requested property succeeds.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return refine(result, options, budget, complete_only=complete_only)


def intervals[k: Int](result: RealCover[k], mut budget: Budget) raises PolyError -> List[Tuple[Rational, Rational, Int]]:
    """Extract real intervals while keeping root counts explicit.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: A real cover; inspect complete and accurate before using it.
        budget: The shared work allowance for preparation and computation.

    Returns:
        (low, high, 1) for each distinct root in increasing order: an open interval containing it, or low == high for an exact rational root. Unresolved pieces appear as (low, high, -1).

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    ref certificate = result.evidence().certificate
    var values = List[Tuple[Rational, Rational, Int]]()
    _allocation[Tuple[Rational, Rational, Int]](len(certificate.leaves) + len(certificate.exact_roots), budget)
    var next_root = 0
    var next = 0
    for leaf in certificate.leaves:
        budget.charge(1)
        if next_root < len(certificate.exact_roots) and certificate.exact_roots[next_root] == leaf.low:
            ref root = certificate.exact_roots[next_root]
            _charge_size(2 * root.size(), budget)
            values.append((root.value.copy(), root.value.copy(), 1))
            next_root += 1
        if not leaf.resolved:
            _charge_size(leaf.low.size(), budget)
            _charge_size(leaf.high.size(), budget)
            values.append((leaf.low.value.copy(), leaf.high.value.copy(), -1))
            continue
        for _ in range(leaf.variations):
            ref bracket = certificate.brackets[next]
            _charge_size(bracket.low.size(), budget)
            _charge_size(bracket.high.size(), budget)
            values.append((bracket.low.value.copy(), bracket.high.value.copy(), 1))
            next += 1
    return values^


def intervals[k: Int](result: RealCover[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Rational, Rational, Int]]:
    """Extract real intervals while keeping root counts explicit.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: A real cover; inspect complete and accurate before using it.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        (low, high, 1) for each distinct root in increasing order: an open interval containing it, or low == high for an exact rational root. Unresolved pieces appear as (low, high, -1).

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return intervals(result, budget)


def rectangles[k: Int](result: ComplexCover[k], mut budget: Budget) raises PolyError -> List[Tuple[Rational, Rational, Rational, Rational, Int]]:
    """Extract complex rectangles with counts including multiplicity.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: A complex cover; inspect its completeness and accuracy.
        budget: The shared work allowance for preparation and computation.

    Returns:
        (real low, real high, imaginary low, imaginary high, count) tuples.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var values = List[Tuple[Rational, Rational, Rational, Rational, Int]]()
    _allocation[Tuple[Rational, Rational, Rational, Rational, Int]](len(result.evidence().regions), budget)
    values.reserve(len(result.evidence().regions))
    for region in result.evidence().regions:
        budget.charge(1)
        ref r = region.rectangle
        _charge_size(r.real_low.size(), budget)
        _charge_size(r.real_high.size(), budget)
        _charge_size(r.imag_low.size(), budget)
        _charge_size(r.imag_high.size(), budget)
        values.append((r.real_low.value.copy(), r.real_high.value.copy(), r.imag_low.value.copy(), r.imag_high.value.copy(), region.count))
    return values^


def rectangles[k: Int](result: ComplexCover[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Rational, Rational, Rational, Rational, Int]]:
    """Extract complex rectangles with counts including multiplicity.

    Parameters:
        k: The inferred coefficient kind.

    Args:
        result: A complex cover; inspect its completeness and accuracy.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        (real low, real high, imaginary low, imaginary high, count) tuples.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return rectangles(result, budget)


from apn_mojo import Complex, ComplexContext
from polynomials import NumericDomain, RootMultiplicityCertificate, real_root_multiplicities, check_real_root_multiplicities, RootDiagnostics, PolishedRootResult, NumericBox
from polynomials import root_diagnostics_exact, root_diagnostics_numeric, check_root_diagnostics_exact, check_root_diagnostics_numeric, polish_root_exact, polish_root_numeric, check_polished_root_exact, check_polished_root_numeric
from . import dense as _dense
from . import _numeric as _num
from ._common import _Kind, _Element
from . import _algebraic_numeric as _an
from polynomials import ComplexBallDomain, PolyErrorKind
from apn_mojo import BallContext


@fieldwise_init
struct Multiplicities[k: Int]( _Archive, Copyable, Movable):
    """Exact multiplicities aligned with the root intervals of a real root cover.

    Parameters:
        k: The inferred exact coefficient kind.
    """
    var _cover: RealCover[Self.k]
    var _proof: ArcPointer[RootMultiplicityCertificate]
    var _values: Values[Int, 0]

    def counts(self) -> Values[Int, 0]:
        """Return shared multiplicities in the same order as roots.intervals.

        Returns:
            Multiplicities aligned with the intervals returned by roots.intervals.
        """
        return self._values.copy()

    def evidence(self) -> ref[ImmOrigin(origin_of(self._proof[]))] RootMultiplicityCertificate:
        """Borrow the squarefree decomposition and the block of each root.

        Returns:
            An immutable reference to the retained core evidence.
        """
        return self._proof[]



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[RealCover[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._cover, writer, sink, budget, depth + 1)
        _result_io._write_node[RootMultiplicityCertificate, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._proof[], writer, sink, budget, depth + 1)
        writer.line(String(len(self._values)), sink, budget)
        for field676 in self._values:
            _result_io._write_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](field676, writer, sink, budget, depth + 1)
        writer.line(String(self._values.modulus()), sink, budget)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field677 = _result_io._read_node[RealCover[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field679 = _result_io._read_node[RootMultiplicityCertificate, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[RootMultiplicityCertificate](1, budget)
        var field678 = ArcPointer(field679^)
        var field682 = reader.count(budget)
        var field681 = List[Int]()
        for _ in range(field682):
            var field683 = _result_io._read_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            _grow(field681, budget)
            field681.append(field683)
        var field680 = _common.Values[Int, 0](field681^, reader.number(UInt64.MAX, budget), budget)
        return Self(field677^, field678^, field680^)

def multiplicities[k: Int](cover: RealCover[k], mut budget: Budget) raises PolyError -> Multiplicities[k]:
    """Attach exact multiplicities to a complete isolated real root cover.

    Parameters:
        k: The inferred exact coefficient kind.

    Args:
        cover: A complete real root cover isolating each distinct root.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        A result retaining the exact multiplicity counts and their witnesses.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var proof = real_root_multiplicities(cover._data[].polynomial, cover._data[].proof, cover._data[].options, budget)
    _allocation[Int](len(proof.multiplicities), budget)
    var counts = List[Int](capacity=len(proof.multiplicities))
    for i in range(len(proof.multiplicities)):
        budget.charge(1)
        counts.append(proof.multiplicities[i])
    _allocation[RootMultiplicityCertificate](1, budget)
    return Multiplicities[k](cover.copy(), ArcPointer(proof^), Values[Int, 0](counts^, 0, budget))


def multiplicities[k: Int](cover: RealCover[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Multiplicities[k]:
    """Return independently checkable multiplicities for the cover's distinct roots.

    Parameters:
        k: The inferred exact coefficient kind.

    Args:
        cover: A complete real root cover isolating each distinct root.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A result retaining the exact multiplicity counts and their witnesses.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return multiplicities(cover, budget)


def check[k: Int](result: Multiplicities[k], mut budget: Budget) raises PolyError -> Bool:
    """Verify the multiplicity evidence and its alignment with the returned counts.

    Parameters:
        k: The inferred exact coefficient kind.

    Args:
        result: The result whose retained evidence or values are being inspected.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Whether the independent evidence supports the stated result.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._cover._data[]
    if not check_real_root_multiplicities(data.polynomial, data.proof, data.options, result._proof[], budget):
        return False
    if len(result._values) != len(result._proof[].multiplicities):
        return False
    for i in range(len(result._values)):
        budget.charge(1)
        if result._values[i] != result._proof[].multiplicities[i]:
            return False
    return True


def check[k: Int](result: Multiplicities[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify exact root multiplicities with one shared work allowance.

    Parameters:
        k: The inferred exact coefficient kind.

    Args:
        result: The result whose retained evidence or values are being inspected.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the independent evidence supports the stated result.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


@fieldwise_init
struct Diagnostics[P: Iterable & Copyable & Deinitable]( _Archive, Copyable, Movable):
    """Residual and local sensitivity bounds retaining the polynomial and test point.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.
    """
    var _source: Self.P
    var _point: ArcPointer[NumericBox]
    var _proof: ArcPointer[RootDiagnostics]
    var _algebraic: Optional[_an._Polynomial[ComplexBallDomain, False]]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._proof[]))] RootDiagnostics:
        """Borrow residual, derivative, backward-error and sensitivity bounds.

        Returns:
            An immutable reference to the retained core evidence.
        """
        return self._proof[]



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = _Ring[_Kind[_Element[Self.P]]]
        comptime _ArchiveDomain = _num._StoredDomain[_Element[Self.P]]
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[Self.P, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._source, writer, sink, budget, depth + 1)
        _result_io._write_node[NumericBox, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._point[], writer, sink, budget, depth + 1)
        _result_io._write_node[RootDiagnostics, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._proof[], writer, sink, budget, depth + 1)
        writer.line(String(Int(Bool(self._algebraic))), sink, budget)
        if self._algebraic:
            _result_io._write_node[_an._Polynomial[ComplexBallDomain, False], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._algebraic.value(), writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[_Kind[_Element[Self.P]]]
        comptime _ArchiveDomain = _num._StoredDomain[_Element[Self.P]]
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field684 = _result_io._read_node[Self.P, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field686 = _result_io._read_node[NumericBox, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[NumericBox](1, budget)
        var field685 = ArcPointer(field686^)
        var field688 = _result_io._read_node[RootDiagnostics, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[RootDiagnostics](1, budget)
        var field687 = ArcPointer(field688^)
        var field689 = Optional[_an._Polynomial[ComplexBallDomain, False]]()
        if reader.number(1, budget):
            var field690 = _result_io._read_node[_an._Polynomial[ComplexBallDomain, False], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            field689 = field690^
        return Self(field684^, field685^, field687^, field689^)

@fieldwise_init
struct Polished[P: Iterable & Copyable & Deinitable]( _Archive, Copyable, Movable):
    """A Newton candidate with a one-root enclosure when independent validation succeeds.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.
    """
    var _source: Self.P
    var _options: NumericOptions
    var _proof: ArcPointer[PolishedRootResult]
    var _algebraic: Optional[_an._Polynomial[ComplexBallDomain, False]]

    def candidate(self) -> Complex:
        """Copy the advisory complex approximation.

        Returns:
            The advisory complex approximation.
        """
        return self._proof[].candidate.copy()

    def outcome(self) -> NumericOutcome:
        """Report the validation and requested-width outcome.

        Returns:
            The retained numerical outcome.
        """
        return self._proof[].outcome

    def evidence(self) -> ref[ImmOrigin(origin_of(self._proof[]))] PolishedRootResult:
        """Borrow the Newton candidate and optional contour certificate.

        Returns:
            An immutable reference to the retained core evidence.
        """
        return self._proof[]



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = _Ring[_Kind[_Element[Self.P]]]
        comptime _ArchiveDomain = _num._StoredDomain[_Element[Self.P]]
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[Self.P, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._source, writer, sink, budget, depth + 1)
        _result_io._write_node[NumericOptions, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._options, writer, sink, budget, depth + 1)
        _result_io._write_node[PolishedRootResult, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._proof[], writer, sink, budget, depth + 1)
        writer.line(String(Int(Bool(self._algebraic))), sink, budget)
        if self._algebraic:
            _result_io._write_node[_an._Polynomial[ComplexBallDomain, False], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._algebraic.value(), writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[_Kind[_Element[Self.P]]]
        comptime _ArchiveDomain = _num._StoredDomain[_Element[Self.P]]
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field691 = _result_io._read_node[Self.P, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field692 = _result_io._read_node[NumericOptions, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field694 = _result_io._read_node[PolishedRootResult, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[PolishedRootResult](1, budget)
        var field693 = ArcPointer(field694^)
        var field695 = Optional[_an._Polynomial[ComplexBallDomain, False]]()
        if reader.number(1, budget):
            var field696 = _result_io._read_node[_an._Polynomial[ComplexBallDomain, False], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            field695 = field696^
        return Self(field691^, field692^, field693^, field695^)

def diagnostics[I: Iterable, X: Copyable & Deinitable](values: I, point: X, mut budget: Budget, *, context: Optional[BallContext] = None) raises PolyError -> Diagnostics[_dense.Array[_Kind[_Element[I]]]] where _Kind[_Element[I]] >= 0:
    """Bound the residual, coefficient backward error and local sensitivity at a point.

    Parameters:
        I: The inferred input iterable type.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        point: The evaluation scalar or ordered coordinate tuple.
        budget: The shared work allowance, including preparation, computation, and result storage.
        context: Explicit ball precision for algebraic coefficients; unused for rational coefficients.

    Returns:
        Residual, derivative, backward-error, and local sensitivity evidence, with the source and point retained.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var p = _dense.asarray(values, budget)
    comptime assert _Kind[_Element[I]] >= 0 and _Kind[_Element[I]] != 3, "root diagnostics require characteristic-zero coefficients"
    var box = _num._bounds(point, budget)
    var enclosed = Optional[_an._Polynomial[ComplexBallDomain, False]]()
    var proof: RootDiagnostics
    comptime if _Kind[_Element[I]] >= 4:
        if not context:
            raise PolyError(PolyErrorKind.DOMAIN, "algebraic diagnostics require an explicit ball context")
        enclosed = _an._convert[_Ring[_Kind[_Element[I]]], ComplexBallDomain, False](p._core[], context.value(), budget)
        proof = root_diagnostics_numeric(enclosed.value().value, box, budget)
    else:
        proof = root_diagnostics_exact(p._core[], box, budget)
    _allocation[NumericBox](1, budget)
    _allocation[RootDiagnostics](1, budget)
    return Diagnostics[_dense.Array[_Kind[_Element[I]]]](p^, ArcPointer(box^), ArcPointer(proof^), enclosed^)


def diagnostics[I: Iterable, X: Copyable & Deinitable](values: I, point: X, *, context: Optional[BallContext] = None, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Diagnostics[_dense.Array[_Kind[_Element[I]]]] where _Kind[_Element[I]] >= 0:
    """Return independently checkable residual and local sensitivity bounds.

    Parameters:
        I: The inferred input iterable type.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        point: The evaluation scalar or ordered coordinate tuple.
        context: Explicit ball precision for algebraic coefficients; unused for rational coefficients.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Residual, derivative, backward-error, and local sensitivity evidence, with the source and point retained.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return diagnostics(values, point, budget, context=context)


def polish[I: Iterable, X: Copyable & Deinitable](values: I, seed: X, mut budget: Budget, *, context: ComplexContext, radius: Rational, iterations: Int = 8, options: NumericOptions = NumericOptions()) raises PolyError -> Polished[_dense.Array[_Kind[_Element[I]]]] where _Kind[_Element[I]] >= 0:
    """Polish a candidate with explicit precision and attempt a one-root contour proof.

    Parameters:
        I: The inferred input iterable type.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        seed: A finite scalar starting point for advisory Newton iteration.
        budget: The shared work allowance, including preparation, computation, and result storage.
        context: Explicit output precision, rounding settings, or ball working precision.
        radius: A positive exact radius for the final validation rectangle.
        iterations: The nonnegative Newton iteration limit; defaults to eight.
        options: Precision, accuracy, isolation, and subdivision requirements.

    Returns:
        An advisory Newton candidate and optional independent one-root certificate.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var p = _dense.asarray(values, budget)
    comptime assert _Kind[_Element[I]] >= 0 and _Kind[_Element[I]] != 3, "root diagnostics require characteristic-zero coefficients"
    var candidate = _num._value[ComplexDomain](seed, context, budget)
    var enclosed = Optional[_an._Polynomial[ComplexBallDomain, False]]()
    var proof: PolishedRootResult
    comptime if _Kind[_Element[I]] >= 4:
        var bits = max(context.real().format().precision(), context.imag().format().precision())
        from .evaluation import _context
        enclosed = _an._convert[_Ring[_Kind[_Element[I]]], ComplexBallDomain, False](p._core[], _context(bits, bits), budget)
        proof = polish_root_numeric(enclosed.value().value, candidate, QQ(radius), iterations, context, options, budget)
    else:
        proof = polish_root_exact(p._core[], candidate, QQ(radius), iterations, context, options, budget)
    _allocation[PolishedRootResult](1, budget)
    return Polished[_dense.Array[_Kind[_Element[I]]]](p^, options, ArcPointer(proof^), enclosed^)


def polish[I: Iterable, X: Copyable & Deinitable](values: I, seed: X, *, context: ComplexContext, radius: Rational, iterations: Int = 8, options: NumericOptions = NumericOptions(), max_steps: UInt64 = UInt64.MAX) raises PolyError -> Polished[_dense.Array[_Kind[_Element[I]]]] where _Kind[_Element[I]] >= 0:
    """Return a Newton candidate and its validation outcome, preserving the input.

    Parameters:
        I: The inferred input iterable type.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        seed: A finite scalar starting point for advisory Newton iteration.
        context: Explicit output precision, rounding settings, or ball working precision.
        radius: A positive exact radius for the final validation rectangle.
        iterations: The nonnegative Newton iteration limit; defaults to eight.
        options: Precision, accuracy, isolation, and subdivision requirements.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An advisory Newton candidate and optional independent one-root certificate.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return polish(values, seed, budget, context=context, radius=radius, iterations=iterations, options=options)


def diagnostics[D: NumericDomain, X: Copyable & Deinitable](values: _num.DenseArray[D], point: X, mut budget: Budget) raises PolyError -> Diagnostics[_num.DenseArray[D]]:
    """Bound the residual, coefficient backward error and local sensitivity at a point.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        point: The evaluation scalar or ordered coordinate tuple.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Residual, derivative, backward-error, and local sensitivity evidence, with the source and point retained.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var p = values.copy()
    var box = _num._bounds(point, budget)
    var proof = root_diagnostics_numeric(p._core[], box, budget)
    _allocation[NumericBox](1, budget)
    _allocation[RootDiagnostics](1, budget)
    return Diagnostics[_num.DenseArray[D]](p^, ArcPointer(box^), ArcPointer(proof^), None)


def diagnostics[D: NumericDomain, X: Copyable & Deinitable](values: _num.DenseArray[D], point: X, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Diagnostics[_num.DenseArray[D]]:
    """Return independently checkable residual and local sensitivity bounds.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        point: The evaluation scalar or ordered coordinate tuple.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Residual, derivative, backward-error, and local sensitivity evidence, with the source and point retained.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return diagnostics[D, X](values, point, budget)


def polish[D: NumericDomain, X: Copyable & Deinitable](values: _num.DenseArray[D], seed: X, mut budget: Budget, *, context: ComplexContext, radius: Rational, iterations: Int = 8, options: NumericOptions = NumericOptions()) raises PolyError -> Polished[_num.DenseArray[D]]:
    """Polish a candidate with explicit precision and attempt a one-root contour proof.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        seed: A finite scalar starting point for advisory Newton iteration.
        budget: The shared work allowance, including preparation, computation, and result storage.
        context: Explicit output precision, rounding settings, or ball working precision.
        radius: A positive exact radius for the final validation rectangle.
        iterations: The nonnegative Newton iteration limit; defaults to eight.
        options: Precision, accuracy, isolation, and subdivision requirements.

    Returns:
        An advisory Newton candidate and optional independent one-root certificate.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var p = values.copy()
    var candidate = _num._value[ComplexDomain](seed, context, budget)
    var proof = polish_root_numeric(p._core[], candidate, QQ(radius), iterations, context, options, budget)
    _allocation[PolishedRootResult](1, budget)
    return Polished[_num.DenseArray[D]](p^, options, ArcPointer(proof^), None)


def polish[D: NumericDomain, X: Copyable & Deinitable](values: _num.DenseArray[D], seed: X, *, context: ComplexContext, radius: Rational, iterations: Int = 8, options: NumericOptions = NumericOptions(), max_steps: UInt64 = UInt64.MAX) raises PolyError -> Polished[_num.DenseArray[D]]:
    """Return a Newton candidate and its validation outcome, preserving the input.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        X: The evaluation point or seed type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        seed: A finite scalar starting point for advisory Newton iteration.
        context: Explicit output precision, rounding settings, or ball working precision.
        radius: A positive exact radius for the final validation rectangle.
        iterations: The nonnegative Newton iteration limit; defaults to eight.
        options: Precision, accuracy, isolation, and subdivision requirements.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An advisory Newton candidate and optional independent one-root certificate.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return polish[D, X](values, seed, budget, context=context, radius=radius, iterations=iterations, options=options)

from polynomials import ComplexDomain


def check[P: Iterable & Copyable & Deinitable](result: Diagnostics[P], mut budget: Budget) raises PolyError -> Bool:
    """Verify the numerical evidence against the retained polynomial and settings.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Whether the independent evidence supports the stated result.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _Kind[_Element[P]] >= 4:
        if not result._algebraic:
            return False
        if not _an._check_polynomial[_Ring[_Kind[_Element[P]]], ComplexBallDomain, False](rebind[_dense.Array[_Kind[_Element[P]]]](result._source)._core[], result._algebraic.value(), budget):
            return False
        return check_root_diagnostics_numeric(result._algebraic.value().value, result._point[], result._proof[], budget)
    elif _Kind[_Element[P]] >= 0:
        return check_root_diagnostics_exact(rebind[_dense.Array[_Kind[_Element[P]]]](result._source)._core[], result._point[], result._proof[], budget)
    else:
        return check_root_diagnostics_numeric(rebind[_num.DenseArray[_num._StoredDomain[_Element[P]]]](result._source)._core[], result._point[], result._proof[], budget)


def check[P: Iterable & Copyable & Deinitable](result: Diagnostics[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Independently check the retained numerical evidence.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the independent evidence supports the stated result.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def check[P: Iterable & Copyable & Deinitable](result: Polished[P], mut budget: Budget) raises PolyError -> Bool:
    """Verify the numerical evidence against the retained polynomial and settings.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Whether the independent evidence supports the stated result.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _Kind[_Element[P]] >= 4:
        if not result._algebraic:
            return False
        if not _an._check_polynomial[_Ring[_Kind[_Element[P]]], ComplexBallDomain, False](rebind[_dense.Array[_Kind[_Element[P]]]](result._source)._core[], result._algebraic.value(), budget):
            return False
        return check_polished_root_numeric(result._algebraic.value().value, result._proof[], result._options, budget)
    elif _Kind[_Element[P]] >= 0:
        return check_polished_root_exact(rebind[_dense.Array[_Kind[_Element[P]]]](result._source)._core[], result._proof[], result._options, budget)
    else:
        return check_polished_root_numeric(rebind[_num.DenseArray[_num._StoredDomain[_Element[P]]]](result._source)._core[], result._proof[], result._options, budget)


def check[P: Iterable & Copyable & Deinitable](result: Polished[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Independently check the retained numerical evidence.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the independent evidence supports the stated result.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


from polynomials import real_roots_numeric, complex_roots_numeric, check_real_roots_numeric, check_complex_roots_numeric, refine_real_roots_numeric, refine_complex_roots_numeric

comptime _NumericProof[z: Bool] = ComplexRootResult if z else RealRootResult


@fieldwise_init
struct NumericCover[D: NumericDomain, complex: Bool]( _Archive, Copyable, Movable):
    """A validated numerical root cover, including unresolved cells and accuracy status.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        complex: Whether this cover uses complex rectangles rather than real intervals.
    """
    var _source: _num.DenseArray[Self.D]
    var _options: NumericOptions
    var _proof: ArcPointer[_NumericProof[Self.complex]]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._proof[]))] _NumericProof[Self.complex]:
        """Borrow the full cover, including exclusions and unresolved regions.

        Returns:
            An immutable reference to the retained core evidence.
        """
        return self._proof[]

    def options(self) -> NumericOptions:
        """Copy the requested precision, accuracy and subdivision settings.

        Returns:
            A copy of the retained numerical requirements.
        """
        return self._options.copy()

    def complete(self) -> Bool:
        """Report whether all roots are accounted for.

        Returns:
            Whether all roots are accounted for.
        """
        comptime if Self.complex:
            return rebind[ComplexRootResult](self._proof[]).complete
        else:
            return rebind[RealRootResult](self._proof[]).complete

    def accurate(self) -> Bool:
        """Report whether counted regions meet the requested width.

        Returns:
            Whether counted regions meet the requested width.
        """
        comptime if Self.complex:
            return rebind[ComplexRootResult](self._proof[]).accurate
        else:
            return rebind[RealRootResult](self._proof[]).accurate

    def isolated(self) -> Bool:
        """Report whether every counted region isolates one root.

        Returns:
            Whether every counted region contains one root.
        """
        comptime if Self.complex:
            return rebind[ComplexRootResult](self._proof[]).isolated
        else:
            return rebind[RealRootResult](self._proof[]).isolated

    def outcome(self) -> NumericOutcome:
        """Report the requested numerical guarantee or unresolved outcome.

        Returns:
            The retained numerical outcome.
        """
        comptime if Self.complex:
            return rebind[ComplexRootResult](self._proof[]).outcome
        else:
            return rebind[RealRootResult](self._proof[]).outcome



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = Self.D
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_num.DenseArray[Self.D], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._source, writer, sink, budget, depth + 1)
        _result_io._write_node[NumericOptions, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._options, writer, sink, budget, depth + 1)
        _result_io._write_node[_NumericProof[Self.complex], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._proof[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = Self.D
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field672 = _result_io._read_node[_num.DenseArray[Self.D], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field673 = _result_io._read_node[NumericOptions, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field675 = _result_io._read_node[_NumericProof[Self.complex], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_NumericProof[Self.complex]](1, budget)
        var field674 = ArcPointer(field675^)
        return Self(field672^, field673^, field674^)

def check[D: NumericDomain, z: Bool](result: NumericCover[D, z], mut budget: Budget) raises PolyError -> Bool:
    """Independently validate a numerical root cover against its retained coefficients.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Whether the evidence supports a complete cover. Accuracy or isolation
        may remain unresolved; inspect those flags separately.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    comptime if z:
        return check_complex_roots_numeric(result._source._core[], rebind[ComplexRootResult](result._proof[]), result._options, budget)
    else:
        return check_real_roots_numeric(result._source._core[], rebind[RealRootResult](result._proof[]), result._options, budget)


def check[D: NumericDomain, z: Bool](result: NumericCover[D, z], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Check all claims of a numerical root cover; incomplete covers return False.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the evidence supports a complete cover. Accuracy or isolation
        may remain unresolved; inspect those flags separately.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def refine[D: NumericDomain, z: Bool](result: NumericCover[D, z], options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> NumericCover[D, z]:
    """Refine a numerical root cover while preserving its coefficient family and previous result.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        options: Precision, accuracy, isolation, and subdivision requirements.
        budget: The shared work allowance, including preparation, computation, and result storage.
        complete_only: Raise unless all requested numerical guarantees are established.

    Returns:
        A new root cover preserving the previous result and coefficient family.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var proof: _NumericProof[z]
    comptime if z:
        proof = rebind_var[_NumericProof[z]](refine_complex_roots_numeric(result._source._core[], rebind[ComplexRootResult](result._proof[]), result._options, options, budget, complete_only=complete_only))
    else:
        proof = rebind_var[_NumericProof[z]](refine_real_roots_numeric(result._source._core[], rebind[RealRootResult](result._proof[]), result._options, options, budget, complete_only=complete_only))
    _allocation[_NumericProof[z]](1, budget)
    return NumericCover[D, z](result._source.copy(), options, ArcPointer(proof^))


def refine[D: NumericDomain, z: Bool](result: NumericCover[D, z], options: NumericOptions, *, max_steps: UInt64 = UInt64.MAX, complete_only: Bool = False) raises PolyError -> NumericCover[D, z]:
    """Refine numerical root enclosures under new accuracy and subdivision settings.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        options: Precision, accuracy, isolation, and subdivision requirements.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        complete_only: Raise unless all requested numerical guarantees are established.

    Returns:
        A new root cover preserving the previous result and coefficient family.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return refine(result, options, budget, complete_only=complete_only)


def intervals[D: NumericDomain](result: NumericCover[D, False], mut budget: Budget) raises PolyError -> List[Tuple[Rational, Rational, Int]]:
    """Extract real intervals while keeping root counts explicit.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        result: A real cover; inspect complete and accurate before using it.
        budget: The shared work allowance for preparation and computation.

    Returns:
        (low, high, distinct count) tuples, omitting zero-count cells. A count of -1 remains unresolved.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var values = List[Tuple[Rational, Rational, Int]]()
    _allocation[Tuple[Rational, Rational, Int]](len(result.evidence().cells), budget)
    values.reserve(len(result.evidence().cells))
    for cell in result.evidence().cells:
        budget.charge(1)
        if cell.count == 0:
            continue
        _charge_size(cell.interval.real_low.size(), budget)
        _charge_size(cell.interval.real_high.size(), budget)
        values.append((cell.interval.real_low.value.copy(), cell.interval.real_high.value.copy(), cell.count))
    return values^


def intervals[D: NumericDomain](result: NumericCover[D, False], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Rational, Rational, Int]]:
    """Extract real intervals while keeping root counts explicit.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        result: A real cover; inspect complete and accurate before using it.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        (low, high, distinct count) tuples, omitting zero-count cells. A count of -1 remains unresolved.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return intervals(result, budget)



def rectangles[D: NumericDomain](result: NumericCover[D, True], mut budget: Budget) raises PolyError -> List[Tuple[Rational, Rational, Rational, Rational, Int]]:
    """Extract complex rectangles with counts including multiplicity.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        result: A complex cover; inspect its completeness and accuracy.
        budget: The shared work allowance for preparation and computation.

    Returns:
        (real low, real high, imaginary low, imaginary high, count) tuples.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var values = List[Tuple[Rational, Rational, Rational, Rational, Int]]()
    _allocation[Tuple[Rational, Rational, Rational, Rational, Int]](len(result.evidence().regions), budget)
    values.reserve(len(result.evidence().regions))
    for region in result.evidence().regions:
        budget.charge(1)
        ref r = region.rectangle
        _charge_size(r.real_low.size(), budget)
        _charge_size(r.real_high.size(), budget)
        _charge_size(r.imag_low.size(), budget)
        _charge_size(r.imag_high.size(), budget)
        values.append((r.real_low.value.copy(), r.real_high.value.copy(), r.imag_low.value.copy(), r.imag_high.value.copy(), region.count))
    return values^


def rectangles[D: NumericDomain](result: NumericCover[D, True], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Rational, Rational, Rational, Rational, Int]]:
    """Extract complex rectangles with counts including multiplicity.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        result: A complex cover; inspect its completeness and accuracy.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        (real low, real high, imaginary low, imaginary high, count) tuples.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return rectangles(result, budget)


from std.iter import iter, next
from std.builtin.rebind import downcast
from . import algebraic as _algebraic
from ._common import _check_batch


def _exact_coefficients[I: Iterable](values: I, mut budget: Budget) raises PolyError -> ExactRoots:
    _check_batch(values)
    var coefficients = List[RealAlgebraic]()
    var cursor = iter(values)
    while True:
        budget.charge(1)
        var value: downcast[_Element[I], Copyable & Deinitable]
        try:
            value = rebind_var[downcast[_Element[I], Copyable & Deinitable]](next(cursor))
        except StopIteration:
            break
        var real = _algebraic.asreal(value, budget)
        _grow(coefficients, budget)
        coefficients.append(_copy_witness(real.witness(), budget))
    var proof = algebraic_roots(coefficients, budget)
    _allocation[Real](len(proof.roots), budget)
    var result = List[Real](capacity=len(proof.roots))
    for root in proof.roots:
        result.append(_plain(_copy_witness(root, budget), budget))
    _allocation[_ExactData](1, budget)
    return ExactRoots(ArcPointer(_ExactData(coefficients^, proof^)), Values[Real, 4](result^, 0, budget))


@fieldwise_init
struct SpecializedRoots[P: Copyable & Deinitable, k: Int = 1]( _Archive, Copyable, Movable, Iterable, Sized):
    """All distinct real roots after specialization, retaining both proofs.

    Parameters:
        P: The original sparse power shape.
        k: The specialization domain; QQ by default, or AA.
    """
    var _specialization: _Specialization[Self.P, Self.k]
    var _roots: ExactRoots
    var _axis: Int
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _SpecializedIterator[Self.P, Self.k, ImmOrigin(iterable_origin)]
    """Borrow the sorted exact roots.

    Parameters:
        iterable_mut: Whether the source origin is mutable.
        iterable_origin: The source origin.
    """

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over distinct roots in increasing order.

        Returns:
            An iterator borrowing the root values.
        """
        return _SpecializedIterator[Self.P, Self.k, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def __len__(self) -> Int:
        """Return the distinct-root count.

        Returns:
            The number of roots, excluding multiplicity.
        """
        return len(self._roots)

    def __getitem__(self, index: Int) -> Real:
        """Read a root at a valid zero-based index.

        Args:
            index: The root's position in increasing order.

        Returns:
            An exact real value sharing its witness.
        """
        return self._roots[index]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._roots._data[].proof))] AlgebraicRootsCertificate:
        """Borrow the complete root certificate.

        Returns:
            The elimination, completeness and filtering evidence.
        """
        return self._roots._data[].proof

    def specialization(self) -> _Specialization[Self.P, Self.k]:
        """Return the retained specialization with its coefficient evidence.

        Returns:
            A value sharing the immutable specialization.
        """
        return self._specialization.copy()



    @staticmethod
    def _archive_shape() -> Int:
        from .sparse import _width
        return _width[Self.P]()

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_Specialization[Self.P, Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._specialization, writer, sink, budget, depth + 1)
        _result_io._write_node[ExactRoots, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._roots, writer, sink, budget, depth + 1)
        _result_io._write_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._axis, writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field697 = _result_io._read_node[_Specialization[Self.P, Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field698 = _result_io._read_node[ExactRoots, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field699 = _result_io._read_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        return Self(field697^, field698^, field699)

@fieldwise_init
struct _SpecializedIterator[P: Copyable & Deinitable, k: Int, origin: ImmOrigin](Iterator):
    comptime Element = Real
    var source: Pointer[SpecializedRoots[Self.P, Self.k], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Real:
        if self.position == len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position]
        self.position += 1
        return result^


def check[P: Copyable & Deinitable, k: Int](result: SpecializedRoots[P, k], mut budget: Budget) raises PolyError -> Bool:
    """Verify specialization, root completeness, and their coefficient binding.

    Parameters:
        P: Inferred from the input types.
        k: The retained specialization domain.

    Args:
        result: The result and its retained evidence.
        budget: The shared work allowance.

    Returns:
        Whether the evidence is valid.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    from .sparse import _width
    if result._axis < 0 or result._axis >= _width[P]():
        return False
    var valid: Bool
    var values: Values[Real, 4]
    comptime if k == 4:
        ref specialized = rebind[_algebraic.AlgebraicSpecialization[P]](result._specialization)
        for assignment in specialized._point[]:
            budget.charge(1)
            if assignment[0] == Int64(result._axis):
                return False
        valid = _algebraic.check(specialized, budget)
        values = _algebraic.coefficients[_algebraic.AlgebraicSpecialization[P], P, Real](specialized, budget, axis=result._axis)
    else:
        ref specialized = rebind[_algebraic.Specialization[P]](result._specialization)
        for assignment in specialized._data[].point:
            budget.charge(1)
            if assignment[0] == Int64(result._axis):
                return False
        valid = _algebraic.check(specialized, budget)
        values = _algebraic.coefficients(specialized, budget, axis=result._axis)
    if not valid or not check(result._roots, budget):
        return False
    ref supplied = result._roots._data[].coefficients
    if len(values) != len(supplied):
        return False
    for i in range(len(values)):
        budget.charge(1)
        if not _same_witness(values[i].witness(), supplied[i], budget):
            return False
    return True


def check[P: Copyable & Deinitable, k: Int](result: SpecializedRoots[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify specialization, root completeness, and their coefficient binding.

    Parameters:
        P: Inferred from the input types.
        k: The retained specialization domain.

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


from . import _complex_roots as _complex
from ._complex_roots import ComplexRoots, ComplexRootsCertificate


def check(result: ComplexRoots, mut budget: Budget) raises PolyError -> Bool:
    """Verify exact complex roots, multiplicities, ordering, and completeness.

    Args:
        result: An exact complex-root result.
        budget: The shared verification work allowance.

    Returns:
        Whether valid component witnesses reconstruct the original polynomial.

    Raises:
        PolyError: DOMAIN for malformed values, or LIMIT/CAPACITY on resource limits.
    """
    return _complex._check(result._data[].polynomial, result.evidence(), budget)


def check(result: ComplexRoots, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify exact complex roots, multiplicities, ordering, and completeness.

    Args:
        result: An exact complex-root result.
        max_steps: The complete check's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether valid component witnesses reconstruct the original polynomial.

    Raises:
        PolyError: DOMAIN for malformed values, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def check_complex[I: Iterable](values: I, certificate: ComplexRootsCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check independently supplied exact complex-root evidence.

    Parameters:
        I: The inferred dense coefficient iterable type.

    Args:
        values: The original ascending exact coefficients or dense array.
        certificate: Component witnesses and claimed root multiplicities.
        budget: The shared preparation and verification allowance.

    Returns:
        Whether the certificate gives all roots with their exact multiplicities.
        A certificate for the zero polynomial is always rejected.

    Raises:
        PolyError: DOMAIN for malformed input, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[_Element[I]] >= 0 and _Kind[_Element[I]] != 3, "exact complex roots require characteristic-zero exact coefficients"
    var p = _dense._promote[5](_dense._prepare(values, budget), NoContext(), budget)
    return _complex._check(p._core[], certificate, budget)


def check_complex[I: Iterable](values: I, certificate: ComplexRootsCertificate, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Check independently supplied exact complex-root evidence.

    Parameters:
        I: The inferred dense coefficient iterable type.

    Args:
        values: The original ascending exact coefficients or dense array.
        certificate: Component witnesses and claimed root multiplicities.
        max_steps: The complete check's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the certificate gives all roots with their exact multiplicities.
        A certificate for the zero polynomial is always rejected.

    Raises:
        PolyError: DOMAIN for malformed input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check_complex(values, certificate, budget)


@fieldwise_init
struct AlgebraicCount( _Archive, Copyable, Movable):
    """An open-interval count retaining all real roots and exact endpoints."""
    var _roots: ExactRoots
    var _low: Real
    var _high: Real
    var _count: Int

    def count(self) -> Int:
        """Return the distinct-root count.

        Returns:
            The number of roots strictly between the retained endpoints.
        """
        return self._count

    def evidence(self) -> ref[ImmOrigin(origin_of(self._roots._data[].proof))] AlgebraicRootsCertificate:
        """Borrow the complete real-root evidence.

        Returns:
            The certificate used to establish the interval count.
        """
        return self._roots._data[].proof



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[ExactRoots, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._roots, writer, sink, budget, depth + 1)
        _result_io._write_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._low, writer, sink, budget, depth + 1)
        _result_io._write_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._high, writer, sink, budget, depth + 1)
        _result_io._write_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._count, writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field664 = _result_io._read_node[ExactRoots, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field665 = _result_io._read_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field666 = _result_io._read_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field667 = _result_io._read_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        return Self(field664^, field665^, field666^, field667)

def _interval_count(result: ExactRoots, low: Real, high: Real, mut budget: Budget) raises PolyError -> Int:
    from . import algebraic
    if algebraic.compare(low, high, budget) >= 0:
        raise PolyError(PolyErrorKind.DOMAIN, "root count requires ordered nonroot endpoints")
    var count = 0
    for value in result:
        var a = algebraic.compare(value, low, budget)
        var b = algebraic.compare(value, high, budget)
        if a == 0 or b == 0:
            raise PolyError(PolyErrorKind.DOMAIN, "root count requires ordered nonroot endpoints")
        if a > 0 and b < 0:
            count += 1
    return count


def check(result: AlgebraicCount, mut budget: Budget) raises PolyError -> Bool:
    """Verify completeness and the exact interval comparisons for a root count.

    Args:
        result: The retained interval count.
        budget: The shared work allowance.

    Returns:
        Whether the complete root evidence and count agree.

    Raises:
        PolyError: LIMIT or CAPACITY on resource limits.
    """
    from . import algebraic
    if not check(result._roots, budget) or not algebraic.check(result._low, budget) or not algebraic.check(result._high, budget):
        return False
    try:
        return result._count == _interval_count(result._roots, result._low, result._high, budget)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


def check(result: AlgebraicCount, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify an algebraic interval count with a bounded work allowance.

    Args:
        result: The retained interval count.
        max_steps: The complete call's work allowance.

    Returns:
        Whether the complete root evidence and count agree.

    Raises:
        PolyError: LIMIT or CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


from . import _common

from ._complex_roots import ExactComplexRootEvent, ExactComplexRootReceipt, ExactComplexRootSink, ExactComplexRootChecker

from . import _algebraic_covers as _ac
from ._algebraic_covers import AlgebraicCover, AlgebraicRealRootChecker, AlgebraicComplexRootChecker


def check[k: Int, z: Bool](result: AlgebraicCover[k, z], mut budget: Budget) raises PolyError -> Bool:
    """Independently validate a numerical root cover against its retained coefficients.

    Parameters:
        k: The inferred exact algebraic coefficient kind.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Whether the evidence supports a complete cover. Accuracy or isolation
        may remain unresolved; inspect those flags separately.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    return _ac._check(result, budget)


def check[k: Int, z: Bool](result: AlgebraicCover[k, z], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Check all claims of a numerical root cover; incomplete covers return False.

    Parameters:
        k: The inferred exact algebraic coefficient kind.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the evidence supports a complete cover. Accuracy or isolation
        may remain unresolved; inspect those flags separately.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def refine[k: Int, z: Bool](result: AlgebraicCover[k, z], options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> AlgebraicCover[k, z]:
    """Re-enclose the original algebraic coefficients and compute a fresh root cover.

    Parameters:
        k: The inferred exact algebraic coefficient kind.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        options: Precision, accuracy, isolation, and subdivision requirements.
        budget: The shared work allowance, including preparation, computation, and result storage.
        complete_only: Raise unless all requested numerical guarantees are established.

    Returns:
        A new cover of the original exact polynomial, with coefficient enclosures
        prepared at options.max_precision. The previous result is unchanged.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    if not _ac.an._check_polynomial[_Ring[k], _ac._Domain[z], False](result._source._core[], result._converted[], budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic coefficient enclosure")
    return _ac._make[k, z](result._source, options, budget, complete_only)


def refine[k: Int, z: Bool](result: AlgebraicCover[k, z], options: NumericOptions, *, max_steps: UInt64 = UInt64.MAX, complete_only: Bool = False) raises PolyError -> AlgebraicCover[k, z]:
    """Refine numerical root enclosures under new accuracy and subdivision settings.

    Parameters:
        k: The inferred exact algebraic coefficient kind.
        z: Whether the cover contains complex rectangles.

    Args:
        result: The result whose retained evidence or values are being inspected.
        options: Precision, accuracy, isolation, and subdivision requirements.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        complete_only: Raise unless all requested numerical guarantees are established.

    Returns:
        A new cover of the original exact polynomial, with coefficient enclosures
        prepared at options.max_precision. The previous result is unchanged.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return refine(result, options, budget, complete_only=complete_only)


def intervals[k: Int](result: AlgebraicCover[k, False], mut budget: Budget) raises PolyError -> List[Tuple[Rational, Rational, Int]]:
    """Extract real intervals while keeping root counts explicit.

    Parameters:
        k: The inferred algebraic coefficient domain.

    Args:
        result: A real cover; inspect complete and accurate before using it.
        budget: The shared work allowance for preparation and computation.

    Returns:
        (low, high, distinct count) tuples, omitting zero-count cells. A count of -1 remains unresolved.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var values = List[Tuple[Rational, Rational, Int]]()
    _allocation[Tuple[Rational, Rational, Int]](len(result.evidence().cells), budget)
    values.reserve(len(result.evidence().cells))
    for cell in result.evidence().cells:
        budget.charge(1)
        if cell.count == 0:
            continue
        _charge_size(cell.interval.real_low.size(), budget)
        _charge_size(cell.interval.real_high.size(), budget)
        values.append((cell.interval.real_low.value.copy(), cell.interval.real_high.value.copy(), cell.count))
    return values^


def intervals[k: Int](result: AlgebraicCover[k, False], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Rational, Rational, Int]]:
    """Extract real intervals while keeping root counts explicit.

    Parameters:
        k: The inferred algebraic coefficient domain.

    Args:
        result: A real cover; inspect complete and accurate before using it.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        (low, high, distinct count) tuples, omitting zero-count cells. A count of -1 remains unresolved.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return intervals(result, budget)



def rectangles[k: Int](result: AlgebraicCover[k, True], mut budget: Budget) raises PolyError -> List[Tuple[Rational, Rational, Rational, Rational, Int]]:
    """Extract complex rectangles with counts including multiplicity.

    Parameters:
        k: The inferred algebraic coefficient domain.

    Args:
        result: A complex cover; inspect its completeness and accuracy.
        budget: The shared work allowance for preparation and computation.

    Returns:
        (real low, real high, imaginary low, imaginary high, count) tuples.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var values = List[Tuple[Rational, Rational, Rational, Rational, Int]]()
    _allocation[Tuple[Rational, Rational, Rational, Rational, Int]](len(result.evidence().regions), budget)
    values.reserve(len(result.evidence().regions))
    for region in result.evidence().regions:
        budget.charge(1)
        ref r = region.rectangle
        _charge_size(r.real_low.size(), budget)
        _charge_size(r.real_high.size(), budget)
        _charge_size(r.imag_low.size(), budget)
        _charge_size(r.imag_high.size(), budget)
        values.append((r.real_low.value.copy(), r.real_high.value.copy(), r.imag_low.value.copy(), r.imag_high.value.copy(), region.count))
    return values^


def rectangles[k: Int](result: AlgebraicCover[k, True], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Rational, Rational, Rational, Rational, Int]]:
    """Extract complex rectangles with counts including multiplicity.

    Parameters:
        k: The inferred algebraic coefficient domain.

    Args:
        result: A complex cover; inspect its completeness and accuracy.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        (real low, real high, imaginary low, imaginary high, count) tuples.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return rectangles(result, budget)



comptime _Specialization[P: Copyable & Deinitable, k: Int] = _algebraic.AlgebraicSpecialization[P] if k == 4 else _algebraic.Specialization[P]

from . import _sturm
from ._sturm import SturmPlan, SturmCount


def check(plan: SturmPlan, mut budget: Budget) raises PolyError -> Bool:
    """Verify the Sturm recurrence and strict degree descent without polynomial division.

    Args:
        plan: The immutable prepared Sturm chain.
        budget: The shared work allowance.

    Returns:
        Whether all chain identities and their source binding are valid.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    return _sturm._check(plan, budget)


def check(plan: SturmPlan, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the Sturm recurrence and strict degree descent without polynomial division.

    Args:
        plan: The immutable prepared Sturm chain.
        max_steps: The complete call allowance.

    Returns:
        Whether all chain identities and their source binding are valid.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    var budget = Budget(max_steps)
    return check(plan, budget)


def check(result: SturmCount, mut budget: Budget) raises PolyError -> Bool:
    """Verify the retained Sturm chain, endpoints, and distinct-root count.

    Args:
        result: The retained interval count.
        budget: The shared work allowance.

    Returns:
        Whether the independent evidence establishes the count.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    return _sturm._check_count(result, budget)


def check(result: SturmCount, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the retained Sturm chain, endpoints, and distinct-root count.

    Args:
        result: The retained interval count.
        max_steps: The complete call allowance.

    Returns:
        Whether the independent evidence establishes the count.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def sign_changes[T: Copyable & Deinitable](plan: SturmPlan, point: T, mut budget: Budget) raises PolyError -> Int:
    """Count sign variations at an integer or rational point, omitting zeros.

    Parameters:
        T: The integer or rational point type.

    Args:
        plan: The immutable prepared Sturm chain.
        point: An integer or rational; roots of chain members are allowed.
        budget: The shared work allowance.

    Returns:
        The number of adjacent opposite nonzero signs.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    return _sturm.sign_changes(plan._data[].sequence, _sturm._point(point, budget), budget)


def sign_changes[T: Copyable & Deinitable](plan: SturmPlan, point: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Int:
    """Count sign variations at an integer or rational point, omitting zeros.

    Parameters:
        T: The integer or rational point type.

    Args:
        plan: The immutable prepared Sturm chain.
        point: An integer or rational; roots of chain members are allowed.
        max_steps: The complete call allowance.

    Returns:
        The number of adjacent opposite nonzero signs.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    var budget = Budget(max_steps)
    return sign_changes(plan, point, budget)


def count_roots[A: Copyable & Deinitable, B: Copyable & Deinitable](plan: SturmPlan, low: A, high: B, mut budget: Budget) raises PolyError -> SturmCount:
    """Count distinct roots in an open interval using a prepared Sturm chain.

    Parameters:
        A: The lower endpoint type.
        B: The upper endpoint type.

    Args:
        plan: The immutable prepared Sturm chain.
        low: The lower integer or rational endpoint; must not be a root.
        high: The strictly larger nonroot endpoint.
        budget: The shared work allowance.

    Returns:
        A checked-result object sharing the plan and retaining both endpoints.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    var a = _sturm._point(low, budget)
    var b = _sturm._point(high, budget)
    var count = _sturm._count(plan, a, b, budget)
    return SturmCount(plan.copy(), a^, b^, count)


def count_roots[A: Copyable & Deinitable, B: Copyable & Deinitable](plan: SturmPlan, low: A, high: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> SturmCount:
    """Count distinct roots in an open interval using a prepared Sturm chain.

    Parameters:
        A: The lower endpoint type.
        B: The upper endpoint type.

    Args:
        plan: The immutable prepared Sturm chain.
        low: The lower integer or rational endpoint; must not be a root.
        high: The strictly larger nonroot endpoint.
        max_steps: The complete call allowance.

    Returns:
        A checked-result object sharing the plan and retaining both endpoints.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    var budget = Budget(max_steps)
    return count_roots(plan, low, high, budget)


def polynomials(plan: SturmPlan, mut budget: Budget) raises PolyError -> List[_dense.Array[1]]:
    """Copy the prepared chain into ordinary rational coefficient arrays.

    Args:
        plan: The immutable prepared Sturm chain.
        budget: The shared work allowance.

    Returns:
        Independent array values in chain order, starting with the input polynomial.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    _allocation[_dense.Array[1]](len(plan), budget)
    var values = List[_dense.Array[1]](capacity=len(plan))
    for p in plan._data[].sequence.polynomials:
        values.append(_dense.from_core(p, budget))
    return values^


def polynomials(plan: SturmPlan, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[_dense.Array[1]]:
    """Copy the prepared chain into ordinary rational coefficient arrays.

    Args:
        plan: The immutable prepared Sturm chain.
        max_steps: The complete call allowance.

    Returns:
        Independent array values in chain order, starting with the input polynomial.

    Raises:
        PolyError: DOMAIN for invalid inputs, or LIMIT/CAPACITY.
    """
    var budget = Budget(max_steps)
    return polynomials(plan, budget)
