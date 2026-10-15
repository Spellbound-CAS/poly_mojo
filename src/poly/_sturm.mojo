"""Reusable rational Sturm chains and checked interval counts."""
from std.memory import ArcPointer
from polynomials import Budget, PolyError, PolyErrorKind, QQ, NoContext, DensePoly, SturmSequence, sturm_sequence, sign_changes, evaluate_dense, derivative_dense, BallDomain, PolynomialTextSource, PolynomialTextSink
from ._common import _allocation, _Kind, _convert
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer


@fieldwise_init
struct _Data(Copyable, Movable):
    var source: DensePoly[QQ]
    var sequence: SturmSequence


@fieldwise_init
struct SturmPlan(_Archive, Copyable, Movable, Sized):
    """An immutable rational Sturm chain shared by repeated interval queries.

    Construct with dense.sturm_sequence or sparse.sturm_sequence. Integer
    coefficients promote to QQ. The zero polynomial has no finite root count.
    """
    var _data: ArcPointer[_Data]

    def __len__(self) -> Int:
        """Return the number of nonzero polynomials in the chain.

        Returns:
            The chain length, including the original polynomial.
        """
        return len(self._data[].sequence.polynomials)

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].sequence))] SturmSequence:
        """Borrow the recurrence and quotient witnesses.

        Returns:
            The immutable core Sturm certificate.
        """
        return self._data[].sequence

    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return _check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io as codec
        codec._depth(depth, budget)
        codec._write_node[DensePoly[QQ], QQ, BallDomain, Int](self._data[].source, writer, sink, budget, depth + 1)
        codec._write_node[SturmSequence, QQ, BallDomain, Int](self._data[].sequence, writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io as codec
        codec._depth(depth, budget)
        var source = codec._read_node[DensePoly[QQ], QQ, BallDomain, Int](reader, budget, depth + 1)
        var sequence = codec._read_node[SturmSequence, QQ, BallDomain, Int](reader, budget, depth + 1)
        _allocation[_Data](1, budget)
        return Self(ArcPointer(_Data(source^, sequence^)))


def _make(var source: DensePoly[QQ], mut budget: Budget) raises PolyError -> SturmPlan:
    var sequence = sturm_sequence(source, budget)
    _allocation[_Data](1, budget)
    return SturmPlan(ArcPointer(_Data(source^, sequence^)))


def _equal(a: DensePoly[QQ], b: DensePoly[QQ], mut budget: Budget) raises PolyError -> Bool:
    if a.variable != b.variable: return False
    return a.sub(b, budget).is_zero()


def _check(plan: SturmPlan, mut budget: Budget) raises PolyError -> Bool:
    ref p = plan._data[].source
    ref proof = plan._data[].sequence
    p.check_invariants(budget)
    if p.is_zero() or p.variable != 0 or not len(proof.polynomials) or len(proof.quotients) != len(proof.polynomials) - 1:
        return False
    for i in range(len(proof.polynomials)):
        budget.charge(1)
        ref q = proof.polynomials[i]
        q.check_invariants(budget)
        if q.is_zero() or q.variable != p.variable: return False
        if i and q.degree() >= proof.polynomials[i - 1].degree(): return False
    if not _equal(p, proof.polynomials[0], budget): return False
    if p.degree() == 0: return len(proof.polynomials) == 1
    if len(proof.polynomials) < 2 or not _equal(derivative_dense(p, budget), proof.polynomials[1], budget): return False
    for i in range(1, len(proof.polynomials)):
        ref q = proof.quotients[i - 1]
        q.check_invariants(budget)
        if q.variable != p.variable: return False
        var relation = q.mul(proof.polynomials[i], budget)
        if i + 1 < len(proof.polynomials):
            relation = relation.sub(proof.polynomials[i + 1], budget)
        if not _equal(proof.polynomials[i - 1], relation, budget): return False
    return True


def _point[T: Copyable & Deinitable](point: T, mut budget: Budget) raises PolyError -> QQ:
    comptime assert _Kind[T] == 0 or _Kind[T] == 1, "Sturm queries require integer or rational points"
    return QQ(_convert[1](point, NoContext(), budget))


def _count(plan: SturmPlan, low: QQ, high: QQ, mut budget: Budget) raises PolyError -> Int:
    if low.sub(high, budget).value.sign() >= 0:
        raise PolyError(PolyErrorKind.DOMAIN, "root count requires increasing endpoints")
    ref p = plan._data[].source
    if evaluate_dense(p, low, budget).is_zero() or evaluate_dense(p, high, budget).is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "root count requires nonroot endpoints")
    return sign_changes(plan._data[].sequence, low, budget) - sign_changes(plan._data[].sequence, high, budget)


@fieldwise_init
struct SturmCount(_Archive, Copyable, Movable):
    """An open-interval root count sharing its prepared Sturm evidence."""
    var _plan: SturmPlan
    var _low: QQ
    var _high: QQ
    var _count: Int

    def count(self) -> Int:
        """Return the claimed number of distinct roots.

        Returns:
            The root count, excluding multiplicity and interval endpoints.
        """
        return self._count

    def plan(self) -> SturmPlan:
        """Return the shared immutable chain used for this query.

        Returns:
            The prepared Sturm plan without copying polynomial buffers.
        """
        return self._plan.copy()

    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return _check_count(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io as codec
        codec._depth(depth, budget)
        codec._write_node[SturmPlan, QQ, BallDomain, Int](self._plan, writer, sink, budget, depth + 1)
        codec._write_node[QQ, QQ, BallDomain, Int](self._low, writer, sink, budget, depth + 1)
        codec._write_node[QQ, QQ, BallDomain, Int](self._high, writer, sink, budget, depth + 1)
        codec._write_node[Int, QQ, BallDomain, Int](self._count, writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io as codec
        codec._depth(depth, budget)
        var plan = codec._read_node[SturmPlan, QQ, BallDomain, Int](reader, budget, depth + 1)
        var low = codec._read_node[QQ, QQ, BallDomain, Int](reader, budget, depth + 1)
        var high = codec._read_node[QQ, QQ, BallDomain, Int](reader, budget, depth + 1)
        var count = codec._read_node[Int, QQ, BallDomain, Int](reader, budget, depth + 1)
        return Self(plan^, low^, high^, count)


def _check_count(result: SturmCount, mut budget: Budget) raises PolyError -> Bool:
    if not _check(result._plan, budget) or result._count < 0: return False
    if result._low.sub(result._high, budget).value.sign() >= 0: return False
    ref p = result._plan._data[].source
    if evaluate_dense(p, result._low, budget).is_zero() or evaluate_dense(p, result._high, budget).is_zero(): return False
    return result._count == _count(result._plan, result._low, result._high, budget)
