"""Retained specialization transcripts for real algebraic coefficients."""
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from polynomials import Budget, PolyError, PolyErrorKind, QQ, NoContext, SparsePoly, Term, Monomial, MonomialOrder, RealAlgebraic, AlgebraicSpecializationCertificate, BallDomain, PolynomialTextSource, PolynomialTextSink
from . import algebraic, sparse
from ._algebraic_value import Real, _copy_witness, _same_witness
from ._common import _allocation, _grow, _sum
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer


def _encode[P: Copyable & Deinitable](source: sparse.Array[P, 4], point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> Tuple[SparsePoly[QQ], List[Tuple[Int64, RealAlgebraic]]]:
    source._core[].check_invariants(budget)
    var assignments = List[Tuple[Int64, RealAlgebraic]]()
    var width = sparse._width[P]()
    for entry in point:
        budget.charge(1)
        if entry[0] < 0 or entry[0] >= Int64(width):
            raise PolyError(PolyErrorKind.DOMAIN, "specialization axis is outside the power shape")
        for previous in assignments:
            budget.charge(1)
            if previous[0] == entry[0]:
                raise PolyError(PolyErrorKind.DOMAIN, "specialization axes must be distinct")
        _grow(assignments, budget)
        assignments.append((entry[0], _copy_witness(entry[1], budget)))
    var terms = List[Term[QQ]]()
    var auxiliary = width
    for term in source._core[].terms:
        if not algebraic.check(term.coefficient, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic coefficient witness")
        _allocation[Tuple[Int64, UInt32]](_sum(len(term.monomial.pairs), 1), budget)
        var pairs = List[Tuple[Int64, UInt32]](capacity=_sum(len(term.monomial.pairs), 1))
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] < 0 or pair[0] >= Int64(width):
                raise PolyError(PolyErrorKind.DOMAIN, "unmapped specialization axis")
            pairs.append(pair)
        var coefficient = QQ.one(NoContext())
        ref witness = term.coefficient.witness()
        if witness.rational:
            coefficient = witness.low.copy()
        else:
            pairs.append((Int64(auxiliary), UInt32(1)))
            _grow(assignments, budget)
            assignments.append((Int64(auxiliary), _copy_witness(witness, budget)))
            auxiliary = _sum(auxiliary, 1)
        _grow(terms, budget)
        terms.append(Term[QQ](Monomial(pairs^, budget), coefficient^))
    return SparsePoly[QQ].from_terms(terms^, MonomialOrder(), NoContext(), budget), assignments^


@fieldwise_init
struct AlgebraicSpecialization[P: Copyable & Deinitable](_Archive, Copyable, Movable, Iterable, Sized):
    """Specialized AA coefficients retaining the source, assignments, and proofs.

    Parameters:
        P: The original sparse power shape.
    """
    var _source: sparse.Array[Self.P, 4]
    var _point: ArcPointer[List[Tuple[Int64, RealAlgebraic]]]
    var _lifted: algebraic.Specialization[Self.P]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _Iterator[Self.P, ImmOrigin(iterable_origin)]
    """Borrow powers and exact real coefficients.

    Parameters:
        iterable_mut: Whether the source origin permits mutation.
        iterable_origin: The source lifetime.
    """

    def __len__(self) -> Int:
        """Return the coefficient-group count, including exact zeros.

        Returns:
            The number of specialized coefficient groups.
        """
        return len(self._lifted)

    def __getitem__(self, index: Int) -> Tuple[Self.P, Real]:
        """Copy one power and coefficient pair.

        Args:
            index: A valid zero-based group index.

        Returns:
            The power and exact algebraic coefficient.
        """
        return self._lifted[index]

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over the specialized terms in canonical order.

        Returns:
            A borrowing iterator over coefficient groups.
        """
        return _Iterator[Self.P, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def domain(self) -> String:
        """Return the exact real algebraic coefficient domain.

        Returns:
            AA.
        """
        return "AA"

    def ndim(self) -> Int:
        """Return the retained polynomial axis count.

        Returns:
            The original power-shape width.
        """
        return sparse._width[Self.P]()

    def evidence(self) -> ref[ImmOrigin(origin_of(self._lifted._data[].proof))] AlgebraicSpecializationCertificate:
        """Borrow the coefficient evaluation transcript.

        Returns:
            The rational lifting and algebraic evaluation evidence.
        """
        return self._lifted.evidence()

    @staticmethod
    def _archive_shape() -> Int:
        return sparse._width[Self.P]()

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return _check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io as codec
        codec._depth(depth, budget)
        codec._write_node[sparse.Array[Self.P, 4], Real, BallDomain, Int](self._source, writer, sink, budget, depth + 1)
        writer.line(String(len(self._point[])), sink, budget)
        for entry in self._point[]:
            codec._write_node[Int64, Real, BallDomain, Int](entry[0], writer, sink, budget, depth + 1)
            codec._write_node[RealAlgebraic, Real, BallDomain, Int](entry[1], writer, sink, budget, depth + 1)
        codec._write_node[algebraic.Specialization[Self.P], Real, BallDomain, Int](self._lifted, writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io as codec
        codec._depth(depth, budget)
        var source = codec._read_node[sparse.Array[Self.P, 4], Real, BallDomain, Int](reader, budget, depth + 1)
        var count = reader.count(budget)
        var point = List[Tuple[Int64, RealAlgebraic]]()
        for _ in range(count):
            var axis = codec._read_node[Int64, Real, BallDomain, Int](reader, budget, depth + 1)
            var value = codec._read_node[RealAlgebraic, Real, BallDomain, Int](reader, budget, depth + 1)
            _grow(point, budget)
            point.append((axis, value^))
        var lifted = codec._read_node[algebraic.Specialization[Self.P], Real, BallDomain, Int](reader, budget, depth + 1)
        _allocation[List[Tuple[Int64, RealAlgebraic]]](1, budget)
        return Self(source^, ArcPointer(point^), lifted^)


@fieldwise_init
struct _Iterator[P: Copyable & Deinitable, origin: ImmOrigin](Iterator):
    comptime Element = Tuple[Self.P, Real]
    var source: Pointer[AlgebraicSpecialization[Self.P], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.Element:
        if self.position == len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position]
        self.position += 1
        return result^


def _make[P: Copyable & Deinitable](source: sparse.Array[P, 4], var point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> AlgebraicSpecialization[P]:
    var encoded = _encode(source, point, budget)
    var polynomial = SparsePoly[QQ].zero(MonomialOrder(), NoContext(), budget)
    var assignments = List[Tuple[Int64, RealAlgebraic]]()
    swap(polynomial, encoded[0])
    swap(assignments, encoded[1])
    var lifted = algebraic._specialize[P](polynomial^, assignments^, budget)
    _allocation[List[Tuple[Int64, RealAlgebraic]]](1, budget)
    return AlgebraicSpecialization[P](source.copy(), ArcPointer(point^), lifted^)


def _check[P: Copyable & Deinitable](value: AlgebraicSpecialization[P], mut budget: Budget) raises PolyError -> Bool:
    var encoded = _encode(value._source, value._point[], budget)
    if not encoded[0].sub(value._lifted._data[].polynomial, budget).is_zero():
        return False
    ref supplied = value._lifted._data[].point
    if len(encoded[1]) != len(supplied):
        return False
    for i in range(len(supplied)):
        budget.charge(1)
        if encoded[1][i][0] != supplied[i][0] or not _same_witness(encoded[1][i][1], supplied[i][1], budget):
            return False
    return algebraic.check(value._lifted, budget)
