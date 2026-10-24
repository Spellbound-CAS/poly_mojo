"""Readable polynomial sign constraints, selective solving, and CAD streams."""
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink, BallDomain, QQ
from ._common import _grow
from std.iter import Iterable, iter, next, StopIteration
from std.memory import ArcPointer
from std.collections import Dict
from polynomials import Budget, PolyError, PolyErrorKind, SparsePoly, CADCertificate, CADCoveringsCertificate, CADSink, CADStreamReceipt, CADStreamChecker, CADProjectionCertificate, CADCell, CADFiber
from polynomials import CADFinalCertificate, CADFinalCell, cad_final, check_cad_final
from . import sparse, _cad_algebraic as _aa, _cad_fixed as _cf

comptime _Checker[k: Int] = _aa._Checker if k == 4 else CADStreamChecker
from .algebraic import Real, _plain, _copy_witness
from polynomials import real_root
import apn_mojo as apn
from ._common import _Element, _Nested, _Kind, _allocation, _Ring, _FieldKind


@fieldwise_init
struct Constraint[P: Copyable & Deinitable](Copyable, Movable):
    """A sparse exact-real polynomial and the signs allowed for its value.

    Construct constraints with negative, zero, positive, nonpositive,
    nonnegative, or nonzero. Every constraint compares the polynomial to zero.

    Parameters:
        P: The polynomial's sparse power shape.
    """
    var _polynomial: sparse.Array[Self.P, 1]
    var _mask: UInt8
    var _algebraic: Optional[sparse.Array[Self.P, 4]]


def negative[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be less than zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "CAD requires exact real coefficients"
    var p = sparse._prepare[I, P, C](values, budget)
    from polynomials import NoContext
    comptime if _Kind[C] == 4:
        from polynomials import MonomialOrder
        return Constraint[P](sparse.Array[P, 1](SparsePoly[QQ](MonomialOrder(), NoContext()), budget), UInt8(1), sparse._promote[4](p, NoContext(), budget))
    else:
        return Constraint[P](sparse._promote[1](p, NoContext(), budget), UInt8(1), None)


def negative[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be less than zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return negative[I, P, C](values, budget)


def zero[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be equal to zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "CAD requires exact real coefficients"
    var p = sparse._prepare[I, P, C](values, budget)
    from polynomials import NoContext
    comptime if _Kind[C] == 4:
        from polynomials import MonomialOrder
        return Constraint[P](sparse.Array[P, 1](SparsePoly[QQ](MonomialOrder(), NoContext()), budget), UInt8(2), sparse._promote[4](p, NoContext(), budget))
    else:
        return Constraint[P](sparse._promote[1](p, NoContext(), budget), UInt8(2), None)


def zero[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be equal to zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return zero[I, P, C](values, budget)


def positive[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be greater than zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "CAD requires exact real coefficients"
    var p = sparse._prepare[I, P, C](values, budget)
    from polynomials import NoContext
    comptime if _Kind[C] == 4:
        from polynomials import MonomialOrder
        return Constraint[P](sparse.Array[P, 1](SparsePoly[QQ](MonomialOrder(), NoContext()), budget), UInt8(4), sparse._promote[4](p, NoContext(), budget))
    else:
        return Constraint[P](sparse._promote[1](p, NoContext(), budget), UInt8(4), None)


def positive[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be greater than zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return positive[I, P, C](values, budget)


def nonpositive[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be less than or equal to zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "CAD requires exact real coefficients"
    var p = sparse._prepare[I, P, C](values, budget)
    from polynomials import NoContext
    comptime if _Kind[C] == 4:
        from polynomials import MonomialOrder
        return Constraint[P](sparse.Array[P, 1](SparsePoly[QQ](MonomialOrder(), NoContext()), budget), UInt8(3), sparse._promote[4](p, NoContext(), budget))
    else:
        return Constraint[P](sparse._promote[1](p, NoContext(), budget), UInt8(3), None)


def nonpositive[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be less than or equal to zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return nonpositive[I, P, C](values, budget)


def nonnegative[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be greater than or equal to zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "CAD requires exact real coefficients"
    var p = sparse._prepare[I, P, C](values, budget)
    from polynomials import NoContext
    comptime if _Kind[C] == 4:
        from polynomials import MonomialOrder
        return Constraint[P](sparse.Array[P, 1](SparsePoly[QQ](MonomialOrder(), NoContext()), budget), UInt8(6), sparse._promote[4](p, NoContext(), budget))
    else:
        return Constraint[P](sparse._promote[1](p, NoContext(), budget), UInt8(6), None)


def nonnegative[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be greater than or equal to zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return nonnegative[I, P, C](values, budget)


def nonzero[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be different from zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "CAD requires exact real coefficients"
    var p = sparse._prepare[I, P, C](values, budget)
    from polynomials import NoContext
    comptime if _Kind[C] == 4:
        from polynomials import MonomialOrder
        return Constraint[P](sparse.Array[P, 1](SparsePoly[QQ](MonomialOrder(), NoContext()), budget), UInt8(5), sparse._promote[4](p, NoContext(), budget))
    else:
        return Constraint[P](sparse._promote[1](p, NoContext(), budget), UInt8(5), None)


def nonzero[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Constraint[P] where _Element[I] == Tuple[P, C]:
    """Require a polynomial value to be different from zero.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        values: An integer, rational, or real-algebraic sparse term iterable or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A reusable constraint for solve; inputs remain unchanged.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return nonzero[I, P, C](values, budget)


def _inputs[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> _aa._Prepared where _Element[_Nested[I]] == Tuple[P, C]:
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "CAD requires exact real coefficients"
    comptime if _Kind[C] == 4:
        var prepared = sparse._system[4, I, P, C](values, budget)
        _allocation[SparsePoly[Real]](len(prepared[0]), budget)
        var originals = List[SparsePoly[Real]](capacity=len(prepared[0]))
        for p in prepared[0]:
            originals.append(p._core[].reorder(p._core[].order, budget))
        return _aa._encode(originals^, sparse._width[P](), budget)
    else:
        var prepared = sparse._system[1, I, P, C](values, budget)
        _allocation[SparsePoly[QQ]](len(prepared[0]), budget)
        var result = List[SparsePoly[QQ]](capacity=len(prepared[0]))
        for p in prepared[0]:
            result.append(p._core[].reorder(p._core[].order, budget))
        return _aa._Prepared(result^, None)


@fieldwise_init
struct _Constraints(Movable):
    var problem: _aa._Prepared
    var masks: List[UInt8]


def _constraints[I: Iterable, P: Copyable & Deinitable](constraints: I, mut budget: Budget) raises PolyError -> _Constraints where _Element[I] == Constraint[P]:
    from polynomials import NoContext
    var items = List[Constraint[P]]()
    var masks = List[UInt8]()
    var algebraic = False
    var cursor = iter(constraints)
    while True:
        budget.charge(1)
        var item: type_of(cursor).Element
        try:
            item = next(cursor)
        except StopIteration:
            break
        var constraint = rebind_var[Constraint[P]](item^)
        algebraic = algebraic or Bool(constraint._algebraic)
        _grow(masks, budget)
        masks.append(constraint._mask)
        _grow(items, budget)
        items.append(constraint^)
    if algebraic:
        _allocation[SparsePoly[Real]](len(items), budget)
        var originals = List[SparsePoly[Real]](capacity=len(items))
        for item in items:
            var p = item._algebraic.value().copy() if item._algebraic else sparse._promote[4](item._polynomial, NoContext(), budget)
            originals.append(p._core[].reorder(p._core[].order, budget))
        var prepared = _aa._encode(originals^, sparse._width[P](), budget)
        _aa._masks(masks, prepared.embedding, budget)
        return _Constraints(prepared^, masks^)
    else:
        _allocation[SparsePoly[QQ]](len(items), budget)
        var inputs = List[SparsePoly[QQ]](capacity=len(items))
        for item in items:
            ref p = item._polynomial._core[]
            inputs.append(p.reorder(p.order, budget))
        return _Constraints(_aa._Prepared(inputs^, None), masks^)


@fieldwise_init
struct _DecompositionData(Copyable, Movable):
    var inputs: List[SparsePoly[QQ]]
    var proof: CADFinalCertificate
    var embedding: Optional[ArcPointer[_aa._Embedding]]
    var fixed: List[Tuple[Int, Real]]


@fieldwise_init
struct Decomposition[P: Copyable & Deinitable]( _Archive, Copyable, Movable):
    """A sign-invariant decomposition at the supplied coefficients and fixed coordinates.

    Parameters:
        P: The input polynomials' power shape.
    """
    var _data: ArcPointer[_DecompositionData]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] CADFinalCertificate:
        """Borrow the decomposition and its proof on the selected coordinate slice.

        Returns:
            Immutable compact final evidence; use cad.check to verify its scope and cad.cells to request absolute geometry.
        """
        return self._data[].proof



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
        _result_io._write_node[_DecompositionData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field703 = _result_io._read_node[_DecompositionData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_DecompositionData](1, budget)
        var field702 = ArcPointer(field703^)
        return Self(field702^)

@fieldwise_init
struct _SolutionData(Copyable, Movable):
    var inputs: List[SparsePoly[QQ]]
    var masks: List[UInt8]
    var proof: CADCoveringsCertificate
    var embedding: Optional[ArcPointer[_aa._Embedding]]
    var fixed: List[Tuple[Int, Real]]


@fieldwise_init
struct Solution[P: Copyable & Deinitable]( _Archive, Copyable, Movable):
    """A witness or a proof that a conjunction has no solution on the selected slice.

    Parameters:
        P: The constraints' shared sparse power shape.
    """
    var _data: ArcPointer[_SolutionData]

    def satisfiable(self) -> Bool:
        """Report whether the search found a point satisfying every constraint.

        Returns:
            True for a witness, False for a covering proof of failure.
        """
        return self._data[].proof.satisfiable()

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] CADCoveringsCertificate:
        """Borrow the coverings certificate: a satisfying point or a covering tree.

        Returns:
            Immutable evidence for this conjunction, not a full decomposition.
        """
        return self._data[].proof



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
        _result_io._write_node[_SolutionData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field705 = _result_io._read_node[_SolutionData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_SolutionData](1, budget)
        var field704 = ArcPointer(field705^)
        return Self(field704^)

def decompose[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](polynomials: I, mut budget: Budget, *, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Decomposition[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Decompose real space into cells where every input has a constant sign.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        budget: The shared work allowance for preparation and computation.
        axes: All axes in base-to-top order; defaults to increasing axis order.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return decompose[I, P, C, List[Tuple[Int, Int]]](polynomials, budget, axes=axes, lazard=lazard, fixed=List[Tuple[Int, Int]]())


def decompose[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](polynomials: I, *, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Decomposition[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Decompose real space into cells where every input has a constant sign.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order; defaults to increasing axis order.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return decompose[I, P, C](polynomials, budget, axes=axes, lazard=lazard)


def solve[I: Iterable, P: Copyable & Deinitable](constraints: I, mut budget: Budget, *, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> Solution[P] where _Element[I] == Constraint[P]:
    """Find a point satisfying every sign constraint, or prove none exists.

    The constraints form a conjunction. Budget exhaustion raises and never
    masquerades as an unsatisfiable result.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        budget: The shared work allowance for preparation and computation.
        axes: All axes in base-to-top order; defaults to increasing axis order.
        lazard: Use Lazard projection; False selects McCallum projection.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return solve[I, P, List[Tuple[Int, Int]]](constraints, budget, axes=axes, lazard=lazard, fixed=List[Tuple[Int, Int]]())


def solve[I: Iterable, P: Copyable & Deinitable](constraints: I, *, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> Solution[P] where _Element[I] == Constraint[P]:
    """Find a point satisfying every sign constraint, or prove none exists.

    The constraints form a conjunction. Budget exhaustion raises and never
    masquerades as an unsatisfiable result.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order; defaults to increasing axis order.
        lazard: Use Lazard projection; False selects McCallum projection.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return solve[I, P](constraints, budget, axes=axes, lazard=lazard)


def check[P: Copyable & Deinitable](result: Decomposition[P], mut budget: Budget) raises PolyError -> Bool:
    """Verify the retained geometric claim independently.

    Parameters:
        P: The inferred power shape.

    Args:
        result: The original result with its input and certificate.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Whether the complete certificate supports the stated result.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    var variables = data.proof.variables(budget)
    if not _cf._check(data.inputs, data.embedding, data.fixed, variables, sparse._width[P](), budget):
        return False
    return check_cad_final(data.inputs, variables, List[UInt8](), data.proof, budget, fixed=_cf._prefix(data.embedding, data.fixed, budget))


def check[P: Copyable & Deinitable](result: Decomposition[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the retained geometric claim independently.

    Parameters:
        P: The inferred power shape.

    Args:
        result: The original result with its input and certificate.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the complete certificate supports the stated result.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def check[P: Copyable & Deinitable](result: Solution[P], mut budget: Budget) raises PolyError -> Bool:
    """Verify the retained geometric claim independently.

    Parameters:
        P: The inferred power shape.

    Args:
        result: The original result with its input and certificate.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Whether the complete certificate supports the stated result.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    if not _cf._check(data.inputs, data.embedding, data.fixed, data.proof.variables, sparse._width[P](), budget) or not _cf._check_masks(data.inputs, data.masks, data.embedding, data.fixed, budget):
        return False
    return _cf._fiber(data.embedding, data.fixed, budget).check(data.inputs, data.masks, data.proof, budget)


def check[P: Copyable & Deinitable](result: Solution[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the retained geometric claim independently.

    Parameters:
        P: The inferred power shape.

    Args:
        result: The original result with its input and certificate.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the complete certificate supports the stated result.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


comptime _OUTPUT_BITS = 830
"""Irrational output coordinates are refined to width at most 2^-_OUTPUT_BITS."""


def _grid_sign(g: List[apn.Integer], j: apn.Integer, k: Int) raises -> Int:
    """The exact sign of G(j / 2^k) for integer G, from 2^(k deg G) G(j / 2^k) by Horner."""
    var n = len(g) - 1
    var acc = g[n].copy()
    for i in range(n - 1, -1, -1):
        acc = acc * j + (g[i] << (k * (n - i)))
    return acc.sign()


def _fixed_values(g: List[apn.Integer], t: apn.Integer, precision: Int) raises -> Tuple[apn.Integer, apn.Integer]:
    """2^precision G(t / 2^precision) and the same for G', floor-rounded."""
    var n = len(g) - 1
    var value = g[n] << precision
    var slope = (g[n] * apn.Integer(n)) << precision if n else apn.Integer(0)
    for i in range(n - 1, -1, -1):
        value = ((value * t) >> precision) + (g[i] << precision)
        if i:
            slope = ((slope * t) >> precision) + ((g[i] * apn.Integer(i)) << precision)
    return (value^, slope^)


def _dyadic_cell(sample: RealAlgebraic, mut work: Int) raises -> Optional[Tuple[apn.Rational, apn.Rational]]:
    """The interval `_output` refines an irrational sample to, without its bisection, or None.

    Refinement to the output width bisects the stored isolating interval
    [L, H]: its validation splits and its narrowing both halve dyadic cells
    of [L, H], and they stop at the first depth k with width (H - L) / 2^k at
    most half the output width. The result is therefore the depth-k cell
    that holds the root strictly inside, whenever no midpoint is the root.
    Here the root is located by safeguarded Newton on G(t) = c f(L + (H - L) t)
    in fixed point, and a cell is accepted only when G has opposite nonzero
    signs at its ends, which also excludes every coarser midpoint.
    """
    var low = sample.low.value.copy()
    var span = sample.high.value - low
    # Refinement returns an interval narrower than the output width unchanged.
    if span < apn.Rational(apn.Integer(1), apn.Integer(1) << _OUTPUT_BITS):
        return None
    # The first depth k with span / 2^k <= 2^-(_OUTPUT_BITS + 1).
    var ratio = span * apn.Rational(apn.Integer(1) << (_OUTPUT_BITS + 1))
    var k = max(0, ratio.numerator().magnitude_bit_length() - ratio.denominator().magnitude_bit_length() - 1)
    while ratio.numerator() > ratio.denominator() << k:
        k += 1
    while k > 0 and ratio.numerator() <= ratio.denominator() << (k - 1):
        k -= 1
    # G(t) = f(L + span t) by Horner over polynomials in t, then integral.
    ref coefficients = sample.polynomial.coefficients
    var n = len(coefficients) - 1
    if n < 1:
        return None
    work += n * n + 1
    var shifted: List[apn.Rational] = [coefficients[n].value.copy()]
    for i in range(n - 1, -1, -1):
        var next = List[apn.Rational](length=len(shifted) + 1, fill=apn.Rational(0))
        for d in range(len(shifted)):
            next[d] = next[d] + shifted[d] * low
            next[d + 1] = next[d + 1] + shifted[d] * span
        next[0] = next[0] + coefficients[i].value
        shifted = next^
    var common = apn.Integer(1)
    for c in shifted:
        common = apn.lcm(common, c.denominator())
    var g = List[apn.Integer](capacity=len(shifted))
    for c in shifted:
        g.append((c * apn.Rational(common)).to_integer_exact())
    # Refinement first validates the interval: Descartes' bound of exactly one
    # leaves it whole. The variations of (1 + x)^n G(1 / (1 + x)), the bound on
    # (0, 1), must be exactly one here; a nonsquarefree f only adds fallbacks.
    var reflected = List[apn.Integer](capacity=n + 1)
    for i in range(n, -1, -1):
        reflected.append(g[i].copy())
    for i in range(n):
        for d in range(n - 1, i - 1, -1):
            reflected[d] = reflected[d] + reflected[d + 1]
    var variations = 0
    var previous = 0
    for c in reflected:
        var sign = c.sign()
        if sign != 0:
            if previous != 0 and sign != previous:
                variations += 1
            previous = sign
    if variations != 1:
        return None
    # The root lies in (0, 1), where G changes sign.
    var below = g[0].sign()
    var whole = apn.Integer(0)
    for c in g:
        whole = whole + c
    if below == 0 or whole.sign() == 0 or whole.sign() == below:
        return None
    var precision = k + 64
    var lo = apn.Integer(0)
    var hi = apn.Integer(1) << precision
    var t = hi >> 1
    var tolerance = apn.Integer(1) << max(0, precision - k - 16)
    var converged = False
    for _ in range(4 * precision + 64):
        work += n + 1
        var values = _fixed_values(g, t, precision)
        var s = values[0].sign()
        if s == 0:
            converged = True
            break
        if s == below:
            lo = t.copy()
        else:
            hi = t.copy()
        var candidate = (lo + hi) >> 1
        if values[1].sign() != 0:
            var newton = t - (values[0] << precision) // values[1]
            if newton > lo and newton < hi:
                candidate = newton^
        var step = candidate - t
        if step.sign() < 0:
            step = apn.Integer(0) - step
        t = candidate^
        if step < tolerance or hi - lo < tolerance:
            converged = True
            break
    if not converged:
        return None
    var cells = apn.Integer(1) << k
    var j = t >> (precision - k)
    # Rounding can leave the estimate one cell off; the exact signs decide.
    for offset in [0, -1, 1]:
        var cell = j + apn.Integer(offset)
        if cell.sign() < 0 or not (cell < cells):
            continue
        work += 2 * n + 2
        var left = _grid_sign(g, cell, k)
        var right = _grid_sign(g, cell + apn.Integer(1), k)
        if left != 0 and right != 0 and left != right:
            var step = span / apn.Rational(cells)
            var a = low + step * apn.Rational(cell)
            var b = a + step
            return (a^, b^)
    return None


def _dyadic_output(sample: RealAlgebraic, mut budget: Budget) raises PolyError -> Optional[Real]:
    """`_output` of an irrational sample from `_dyadic_cell`, or None to refine it."""
    var work = 0
    var cell: Optional[Tuple[apn.Rational, apn.Rational]]
    try:
        cell = _dyadic_cell(sample, work)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    budget.charge(UInt64(work))
    if not cell:
        return None
    ref ends = cell.value()
    return _plain(real_root(sample.polynomial, QQ(ends[0].copy()), QQ(ends[1].copy()), budget), budget)


def _output(sample: RealAlgebraic, mut budget: Budget) raises PolyError -> Real:
    # A plain copy of a stored sample at the output precision; lifting keeps the stored interval.
    if not sample.rational and _OUTPUT_BITS != 0:
        var direct = _dyadic_output(sample, budget)
        if direct:
            return direct.take()
    var value = _plain(_copy_witness(sample, budget), budget)
    if sample.rational or _OUTPUT_BITS == 0:
        return value^
    var width: apn.Rational
    try:
        width = apn.Rational(apn.Integer(1), apn.Integer(1) << apn.Integer(_OUTPUT_BITS))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return _alg.refine(value, width, budget)


def _sample_key(sample: RealAlgebraic, mut budget: Budget) raises PolyError -> String:
    """The defining polynomial and isolating interval, all that `_output` reads of an irrational sample."""
    var key = String()
    for i in range(len(sample.polynomial.coefficients)):
        key.write(sample.polynomial.coefficients[i].value, ",")
    key.write(":", sample.low.value, ":", sample.high.value)
    budget.charge(UInt64(key.byte_length() // 64 + 1))
    return key^


def _point[P: Copyable & Deinitable](stored_cells: List[CADFinalCell], variables: List[Int64], index: Int, mut budget: Budget) raises PolyError -> List[Real]:
    var outputs = List[Optional[Real]](length=len(stored_cells), fill=None)
    var refined = Dict[String, Int]()
    return _point_cached[P](stored_cells, variables, index, outputs, refined, budget)


def _point_cached[P: Copyable & Deinitable](stored_cells: List[CADFinalCell], variables: List[Int64], index: Int, mut outputs: List[Optional[Real]], mut refined: Dict[String, Int], mut budget: Budget) raises PolyError -> List[Real]:
    var width = sparse._width[P]()
    _allocation[Int](width, budget)
    var indices = List[Int](length=width, fill=-1)
    var cursor = index
    while cursor >= 0 and stored_cells[cursor].level >= 0:
        budget.charge(1)
        ref cell = stored_cells[cursor]
        var axis = variables[cell.level]
        if axis >= 0 and axis < Int64(width):
            indices[Int(axis)] = cursor
        cursor = cell.parent
    _allocation[Real](width, budget)
    var result = List[Real](capacity=width)
    for entry in indices:
        budget.charge(1)
        if entry < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "sample path does not include every axis")
        if not outputs[entry]:
            ref sample = stored_cells[entry].sample
            if sample.rational:
                outputs[entry] = _output(sample, budget)
            else:
                # Cells in different stacks often share a root; its output is
                # a function of the key alone, so it is refined once.
                var key = _sample_key(sample, budget)
                var earlier = refined.get(key)
                if earlier:
                    outputs[entry] = outputs[earlier.value()].value().copy()
                else:
                    outputs[entry] = _output(sample, budget)
                    refined[key^] = entry
        result.append(outputs[entry].value().copy())
    return result^


def sample[P: Copyable & Deinitable](result: Solution[P], mut budget: Budget) raises PolyError -> Optional[List[Real]]:
    """Extract a satisfying point in polynomial axis order.

    Parameters:
        P: The inferred power shape.

    Args:
        result: A selective sign-solving result.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Exact real coordinates, or None when the conjunction is unsatisfiable.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    if not result.satisfiable():
        budget.charge(1)
        return None
    ref proof = result._data[].proof
    var width = sparse._width[P]()
    _allocation[Int](width, budget)
    var indices = List[Int](length=width, fill=-1)
    for k in range(len(proof.variables)):
        budget.charge(1)
        var axis = proof.variables[k]
        if axis >= 0 and axis < Int64(width) and k < len(proof.witness):
            indices[Int(axis)] = k
    _allocation[Real](width, budget)
    var point = List[Real](capacity=width)
    for entry in indices:
        budget.charge(1)
        if entry < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "sample does not include every axis")
        point.append(_output(proof.witness[entry], budget))
    return Optional(point^)


def sample[P: Copyable & Deinitable](result: Solution[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Optional[List[Real]]:
    """Extract a satisfying point in polynomial axis order.

    Parameters:
        P: The inferred power shape.

    Args:
        result: A selective sign-solving result.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Exact real coordinates, or None when the conjunction is unsatisfiable.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return sample(result, budget)


@fieldwise_init
struct Cell(Copyable, Movable):
    """A terminal cell's dimension, exact sample point, and input signs."""
    var index: Int
    """The corresponding cell index in the full decomposition certificate."""
    var dimension: Int
    """The number of sector choices along this cell's path."""
    var sample: List[Real]
    """An exact sample point in polynomial axis order."""
    var signs: List[Int]
    """Minus one, zero, or one in original input polynomial order."""


def cells[P: Copyable & Deinitable](result: Decomposition[P], mut budget: Budget) raises PolyError -> List[Cell]:
    """Describe every terminal cell using its dimension, sample, and signs.

    Parameters:
        P: The inferred power shape.

    Args:
        result: A full decomposition.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Cells in certificate order. Samples name cells; the full certificate proves their geometry.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return _describe_cells[P](result._data[].proof, result._data[].embedding, result._data[].fixed, budget)


def cells[P: Copyable & Deinitable](result: Decomposition[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Cell]:
    """Describe every terminal cell using its dimension, sample, and signs.

    Parameters:
        P: The inferred power shape.

    Args:
        result: A full decomposition.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Cells in certificate order. Samples name cells; the full certificate proves their geometry.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return cells(result, budget)


def decompose_to[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, S: CADSink](polynomials: I, mut sink: S, mut budget: Budget, *, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> CADStreamReceipt where _Element[_Nested[I]] == Tuple[P, C]:
    """Stream projection and cell witnesses without retaining their full history.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        S: The inferred input type.

    Args:
        polynomials: The integer, rational, or real-algebraic sparse polynomial inputs.
        sink: The consumer of projection and depth-first witness events.
        budget: The shared work allowance for preparation and computation.
        axes: All axes in base-to-top order.
        lazard: Use Lazard projection by default.

    Returns:
        A completion receipt. Events written before an error do not prove a complete decomposition.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return decompose_to[I, P, C, S, List[Tuple[Int, Int]]](polynomials, sink, budget, axes=axes, lazard=lazard, fixed=List[Tuple[Int, Int]]())


def decompose_to[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, S: CADSink](polynomials: I, mut sink: S, *, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> CADStreamReceipt where _Element[_Nested[I]] == Tuple[P, C]:
    """Stream projection and cell witnesses without retaining their full history.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        S: The inferred input type.

    Args:
        polynomials: The integer, rational, or real-algebraic sparse polynomial inputs.
        sink: The consumer of projection and depth-first witness events.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order.
        lazard: Use Lazard projection by default.

    Returns:
        A completion receipt. Events written before an error do not prove a complete decomposition.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return decompose_to[I, P, C, S](polynomials, sink, budget, axes=axes, lazard=lazard)


def checker[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](polynomials: I, mut budget: Budget) raises PolyError -> _Checker[_Kind[C]] where _Element[_Nested[I]] == Tuple[P, C]:
    """Prepare independent incremental checking of a decomposition stream.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        polynomials: The original input polynomials.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A checker accepting streamed events; finish(receipt, budget) verifies complete traversal.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _inputs[I, P, C](polynomials, budget)
    var core = _aa._fiber(prepared.embedding, budget).checker(prepared.inputs, budget)
    comptime if _Kind[C] == 4:
        return rebind_var[_Checker[_Kind[C]]](_aa._Checker(core^, prepared.embedding.copy(), sparse._width[P](), False))
    else:
        return rebind_var[_Checker[_Kind[C]]](core^)


def checker[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](polynomials: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Checker[_Kind[C]] where _Element[_Nested[I]] == Tuple[P, C]:
    """Prepare independent incremental checking of a decomposition stream.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.

    Args:
        polynomials: The original input polynomials.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A checker accepting streamed events; finish(receipt, budget) verifies complete traversal.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return checker[I, P, C](polynomials, budget)


from std.builtin.rebind import downcast
from polynomials import RealAlgebraic, CADCellCertificate, CADCoveringCertificate, LazardEvaluationCertificate
from polynomials import cad_projection, check_cad_projection, cad_cell, check_cad_cell, lazard_evaluate as _lazard_evaluate, check_lazard_evaluation
from . import algebraic as _alg
from ._interop import _sparse_copy, _variables


@fieldwise_init
struct _ProjectionData(Copyable, Movable):
    var inputs: List[SparsePoly[QQ]]
    var variables: List[Int64]
    var proof: CADProjectionCertificate
    var embedding: Optional[ArcPointer[_aa._Embedding]]


@fieldwise_init
struct _CellAtData(Copyable, Movable):
    var inputs: List[SparsePoly[QQ]]
    var variables: List[Int64]
    var point: List[RealAlgebraic]
    var proof: CADCellCertificate
    var embedding: Optional[ArcPointer[_aa._Embedding]]


@fieldwise_init
struct _LazardData(Copyable, Movable):
    var polynomial: SparsePoly[QQ]
    var point: List[Tuple[Int64, RealAlgebraic]]
    var axis: Int64
    var proof: LazardEvaluationCertificate
    var embedding: Optional[ArcPointer[_aa._Embedding]]


@fieldwise_init
struct _CoveringData(Copyable, Movable):
    var inputs: List[SparsePoly[QQ]]
    var masks: List[UInt8]
    var proof: CADFinalCertificate
    var embedding: Optional[ArcPointer[_aa._Embedding]]
    var fixed: List[Tuple[Int, Real]]


@fieldwise_init
struct Projection[P: Copyable & Deinitable, k: Int = 1]( _Archive, Copyable, Movable):
    """A retained CAD projection closure in an explicit axis order.

    Parameters:
        P: The original sparse power shape.
        k: The inferred real coefficient field; normally left implicit.
    """
    var _data: ArcPointer[_ProjectionData]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] CADProjectionCertificate:
        """Borrow the retained core evidence.

        Returns:
            An immutable certificate bound to the original inputs.
        """
        return self._data[].proof



    @staticmethod
    def _archive_shape() -> Int:
        from .sparse import _width
        return _width[Self.P]()

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_ProjectionData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field707 = _result_io._read_node[_ProjectionData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_ProjectionData](1, budget)
        var field706 = ArcPointer(field707^)
        return Self(field706^)

@fieldwise_init
struct CellAt[P: Copyable & Deinitable]( _Archive, Copyable, Movable):
    """A certified sign-invariant cell containing the supplied exact point.

    Parameters:
        P: The original sparse power shape.
    """
    var _data: ArcPointer[_CellAtData]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] CADCellCertificate:
        """Borrow the retained core evidence.

        Returns:
            An immutable certificate bound to the original inputs.
        """
        return self._data[].proof



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
        _result_io._write_node[_CellAtData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field709 = _result_io._read_node[_CellAtData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_CellAtData](1, budget)
        var field708 = ArcPointer(field709^)
        return Self(field708^)

@fieldwise_init
struct Lazard[P: Copyable & Deinitable]( _Archive, Copyable, Movable):
    """A lifting polynomial with certified vanishing orders.

    Parameters:
        P: The original sparse power shape.
    """
    var _data: ArcPointer[_LazardData]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] LazardEvaluationCertificate:
        """Borrow the retained core evidence.

        Returns:
            An immutable certificate bound to the original inputs.
        """
        return self._data[].proof



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
        _result_io._write_node[_LazardData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field711 = _result_io._read_node[_LazardData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_LazardData](1, budget)
        var field710 = ArcPointer(field711^)
        return Self(field710^)

@fieldwise_init
struct Covering[P: Copyable & Deinitable]( _Archive, Copyable, Movable):
    """A full-decomposition feasibility witness or conflict covering.

    Parameters:
        P: The original sparse power shape.
    """
    var _data: ArcPointer[_CoveringData]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] CADFinalCertificate:
        """Borrow the retained core evidence.

        Returns:
            An immutable certificate bound to the original inputs.
        """
        return self._data[].proof

    def satisfiable(self) -> Bool:
        """Report whether the covering contains a satisfying witness cell.

        Returns:
            True for a witness; verify the claim with cad.check.
        """
        return self._data[].proof.satisfiable()



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
        _result_io._write_node[_CoveringData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field713 = _result_io._read_node[_CoveringData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_CoveringData](1, budget)
        var field712 = ArcPointer(field713^)
        return Self(field712^)

def project[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](polynomials: I, mut budget: Budget, *, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> Projection[P, _FieldKind[_Kind[C]]] where _Element[_Nested[I]] == Tuple[P, C]:
    """Compute a CAD projection closure without constructing lifted cells.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        budget: The shared work allowance, including preparation and output storage.
        axes: All axes in base-to-top cylindrical order; empty uses increasing order.
        lazard: Use Lazard projection; False selects McCallum projection.

    Returns:
        A projection result with independently checkable closure and basis evidence.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _inputs[I, P, C](polynomials, budget)
    var inputs = List[SparsePoly[QQ]]()
    swap(inputs, prepared.inputs)
    var original_axes = sparse._axes[P](axes, budget)
    var variables = _aa._axes(original_axes, prepared.embedding, budget)
    var proof = cad_projection(inputs, variables, budget, lazard=lazard)
    _allocation[_ProjectionData](1, budget)
    return Projection[P, _FieldKind[_Kind[C]]](ArcPointer(_ProjectionData(inputs^, variables^, proof^, prepared.embedding.copy())))


def project[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](polynomials: I, *, axes: List[Int] = List[Int](), lazard: Bool = True, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Projection[P, _FieldKind[_Kind[C]]] where _Element[_Nested[I]] == Tuple[P, C]:
    """Compute a CAD projection closure without constructing lifted cells.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        axes: All axes in base-to-top cylindrical order; empty uses increasing order.
        lazard: Use Lazard projection; False selects McCallum projection.
        max_steps: The complete call’s work allowance.

    Returns:
        A projection result with independently checkable closure and basis evidence.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return project[I, P, C](polynomials, budget, axes=axes, lazard=lazard)


def _query_point[P: Copyable & Deinitable, A: Iterable](point: A, variables: List[Int64], mut budget: Budget) raises PolyError -> List[RealAlgebraic]:
    var original = List[RealAlgebraic]()
    var cursor = iter(point)
    while True:
        budget.charge(1)
        var item: downcast[_Element[A], Copyable & Deinitable]
        try:
            item = rebind_var[downcast[_Element[A], Copyable & Deinitable]](next(cursor))
        except StopIteration:
            break
        if len(original) == sparse._width[P]():
            raise PolyError(PolyErrorKind.DOMAIN, "CAD point has too many coordinates")
        var value = _alg.asreal(item, budget)
        _grow(original, budget)
        original.append(_copy_witness(value.witness(), budget))
    if len(original) != sparse._width[P]():
        raise PolyError(PolyErrorKind.DOMAIN, "CAD point has the wrong number of coordinates")
    _allocation[RealAlgebraic](len(original), budget)
    var ordered = List[RealAlgebraic](capacity=len(original))
    for axis in variables:
        ordered.append(_copy_witness(original[Int(axis)], budget))
    return ordered^


def cell_at[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, A: Iterable](polynomials: I, point: A, mut budget: Budget, *, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> CellAt[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Construct only the CAD cell path containing an exact query point.

    Coordinates use original polynomial axis order, even when axes changes the CAD traversal order.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.
        A: Inferred from the input types unless explicitly supplied.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        point: The scalar argument or coordinates in original axis order.
        budget: The shared work allowance, including preparation and output storage.
        axes: All axes in base-to-top cylindrical order; empty uses increasing order.
        lazard: Use Lazard projection; False selects McCallum projection.

    Returns:
        A cell result retaining query-point containment and constant-sign evidence.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _inputs[I, P, C](polynomials, budget)
    var inputs = List[SparsePoly[QQ]]()
    swap(inputs, prepared.inputs)
    var original_axes = sparse._axes[P](axes, budget)
    var variables = _aa._axes(original_axes, prepared.embedding, budget)
    var query = _aa._query(_query_point[P](point, original_axes, budget), prepared.embedding, budget)
    var proof = cad_cell(inputs, variables, query, budget, lazard=lazard)
    _allocation[_CellAtData](1, budget)
    return CellAt[P](ArcPointer(_CellAtData(inputs^, variables^, query^, proof^, prepared.embedding.copy())))


def cell_at[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, A: Iterable](polynomials: I, point: A, *, axes: List[Int] = List[Int](), lazard: Bool = True, max_steps: UInt64 = UInt64.MAX) raises PolyError -> CellAt[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Construct only the CAD cell path containing an exact query point.

    Coordinates use original polynomial axis order, even when axes changes the CAD traversal order.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.
        A: Inferred from the input types unless explicitly supplied.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        point: The scalar argument or coordinates in original axis order.
        axes: All axes in base-to-top cylindrical order; empty uses increasing order.
        lazard: Use Lazard projection; False selects McCallum projection.
        max_steps: The complete call’s work allowance.

    Returns:
        A cell result retaining query-point containment and constant-sign evidence.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return cell_at[I, P, C, A](polynomials, point, budget, axes=axes, lazard=lazard)


def cell_at[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, *Ts: Copyable & Deinitable](polynomials: I, point: Tuple[*Ts], mut budget: Budget, *, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> CellAt[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Construct the CAD cell containing a coordinate tuple.

    Tuple positions use original polynomial axis order.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.
        Ts: Inferred from the input types unless explicitly supplied.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        point: The scalar argument or coordinates in original axis order.
        budget: The shared work allowance, including preparation and output storage.
        axes: All axes in base-to-top cylindrical order; empty uses increasing order.
        lazard: Use Lazard projection; False selects McCallum projection.

    Returns:
        A retained cell containing the exact query point.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert Tuple[*Ts].Ts.length == sparse._width[P](), "CAD point has the wrong number of coordinates"
    _allocation[Real](Tuple[*Ts].Ts.length, budget)
    var coordinates = List[Real](capacity=Tuple[*Ts].Ts.length)
    comptime for i in range(Tuple[*Ts].Ts.length):
        coordinates.append(_alg.asreal(point[i], budget))
    return cell_at[I, P, C, List[Real]](polynomials, coordinates, budget, axes=axes, lazard=lazard)


def cell_at[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, *Ts: Copyable & Deinitable](polynomials: I, point: Tuple[*Ts], *, axes: List[Int] = List[Int](), lazard: Bool = True, max_steps: UInt64 = UInt64.MAX) raises PolyError -> CellAt[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Construct the CAD cell containing a coordinate tuple.

    Tuple positions use original polynomial axis order.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.
        Ts: Inferred from the input types unless explicitly supplied.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        point: The scalar argument or coordinates in original axis order.
        axes: All axes in base-to-top cylindrical order; empty uses increasing order.
        lazard: Use Lazard projection; False selects McCallum projection.
        max_steps: The complete call’s work allowance.

    Returns:
        A retained cell containing the exact query point.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return cell_at[I, P, C, *Ts](polynomials, point, budget, axes=axes, lazard=lazard)


def lazard_evaluate[I: Iterable, A: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, X: Copyable & Deinitable](values: I, assignments: A, mut budget: Budget, *, axis: Int = 0) raises PolyError -> Lazard[P] where _Element[I] == Tuple[P, C] and _Element[A] == Tuple[Int, X]:
    """Remove vanishing factors before ordered algebraic substitutions.

    Assignments follow base-to-top cylindrical order and must leave the lifting axis unassigned.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        A: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.
        X: Inferred from the input types unless explicitly supplied.

    Args:
        values: Input coefficients, sparse terms, or a prepared array.
        assignments: Distinct axis and exact-value pairs in the required order.
        budget: The shared work allowance, including preparation and output storage.
        axis: The zero-based polynomial axis.

    Returns:
        Certified vanishing orders and a nonzero lifting polynomial up to positive factorials.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[C] == 0 or _Kind[C] == 1 or _Kind[C] == 4, "Lazard evaluation requires exact real coefficients"
    sparse._axis[P](axis)
    var p = sparse._prepare[I, P, C](values, budget)
    _allocation[type_of(p)](1, budget)
    var family = List[type_of(p)](length=1, fill=p)
    var prepared = _inputs[List[type_of(p)], P, _Element[type_of(p)].Ts[1]](family, budget)
    var point = sparse._algebraic_assignments[P, A, X](assignments, budget)
    var extended = _aa._assignments(point, prepared.embedding, budget)
    var polynomial = prepared.inputs[0].reorder(prepared.inputs[0].order, budget)
    var proof = _lazard_evaluate(polynomial, Int64(axis), extended, budget)
    _allocation[_LazardData](1, budget)
    return Lazard[P](ArcPointer(_LazardData(polynomial^, extended^, Int64(axis), proof^, prepared.embedding.copy())))


def lazard_evaluate[I: Iterable, A: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, X: Copyable & Deinitable](values: I, assignments: A, *, axis: Int = 0, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Lazard[P] where _Element[I] == Tuple[P, C] and _Element[A] == Tuple[Int, X]:
    """Remove vanishing factors before ordered algebraic substitutions.

    Assignments follow base-to-top cylindrical order and must leave the lifting axis unassigned.

    Parameters:
        I: Inferred from the input types unless explicitly supplied.
        A: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        C: Inferred from the input types unless explicitly supplied.
        X: Inferred from the input types unless explicitly supplied.

    Args:
        values: Input coefficients, sparse terms, or a prepared array.
        assignments: Distinct axis and exact-value pairs in the required order.
        axis: The zero-based polynomial axis.
        max_steps: The complete call’s work allowance.

    Returns:
        Certified vanishing orders and a nonzero lifting polynomial up to positive factorials.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return lazard_evaluate[I, A, P, C, X](values, assignments, budget, axis=axis)


def check[P: Copyable & Deinitable, k: Int](result: Projection[P, k], mut budget: Budget) raises PolyError -> Bool:
    """Verify the retained projection evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.
        k: The inferred real coefficient field.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    if not _aa._check(data.inputs, data.embedding, sparse._width[P](), budget) or not _aa._check_axes(data.variables, data.embedding, sparse._width[P](), budget):
        return False
    if len(data.variables) != len(data.proof.variables):
        return False
    for i in range(len(data.variables)):
        budget.charge(1)
        if data.variables[i] != data.proof.variables[i]:
            return False
    return check_cad_projection(data.inputs, data.proof, budget)


def check[P: Copyable & Deinitable, k: Int](result: Projection[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the retained projection evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.
        k: The inferred real coefficient field.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def check[P: Copyable & Deinitable](result: CellAt[P], mut budget: Budget) raises PolyError -> Bool:
    """Verify the retained cellat evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    if not _aa._check(data.inputs, data.embedding, sparse._width[P](), budget) or not _aa._check_axes(data.variables, data.embedding, sparse._width[P](), budget):
        return False
    if len(data.variables) != len(data.proof.projection.variables):
        return False
    for i in range(len(data.variables)):
        budget.charge(1)
        if data.variables[i] != data.proof.projection.variables[i]:
            return False
    if not _aa._check_query(data.point, data.embedding, budget):
        return False
    return check_cad_cell(data.inputs, data.point, data.proof, budget)


def check[P: Copyable & Deinitable](result: CellAt[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the retained cellat evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def check[P: Copyable & Deinitable](result: Lazard[P], mut budget: Budget) raises PolyError -> Bool:
    """Verify the retained lazard evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    if data.axis < 0 or data.axis >= Int64(sparse._width[P]()):
        return False
    if not _aa._check_lazard(data.polynomial, data.point, data.embedding, sparse._width[P](), budget):
        return False
    if data.embedding:
        var count = len(data.embedding.value()[].roots)
        if len(data.proof.orders) < count:
            return False
        for i in range(count):
            budget.charge(1)
            if data.proof.orders[i] != 0:
                return False
    return check_lazard_evaluation(data.polynomial, data.axis, data.point, data.proof, budget)


def check[P: Copyable & Deinitable](result: Lazard[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the retained lazard evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def check[P: Copyable & Deinitable](result: Covering[P], mut budget: Budget) raises PolyError -> Bool:
    """Verify the retained covering evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    var variables = data.proof.variables(budget)
    if not _cf._check(data.inputs, data.embedding, data.fixed, variables, sparse._width[P](), budget) or not _cf._check_masks(data.inputs, data.masks, data.embedding, data.fixed, budget):
        return False
    return check_cad_final(data.inputs, variables, data.masks, data.proof, budget, fixed=_cf._prefix(data.embedding, data.fixed, budget), query=True)


def check[P: Copyable & Deinitable](result: Covering[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the retained covering evidence against its original inputs.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Whether all claimed relations and required bindings are valid.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def projection_polynomials[P: Copyable & Deinitable, k: Int](result: Projection[P, k], mut budget: Budget, *, axis: Int = -1) raises PolyError -> List[sparse.Array[P, k]]:
    """Copy projection basis polynomials for one axis or the complete closure.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.
        k: The inferred real coefficient field.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.
        axis: The original polynomial axis; -1 collects every level in CAD order.

    Returns:
        Independent sparse arrays retaining all original axes.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    if axis < -1 or axis >= sparse._width[P]():
        raise PolyError(PolyErrorKind.DOMAIN, "projection axis outside the power shape")
    var values = List[sparse.Array[P, k]]()
    ref proof = result._data[].proof
    for i in range(len(proof.levels)):
        budget.charge(1)
        if proof.variables[i] >= Int64(sparse._width[P]()) or (axis != -1 and proof.variables[i] != Int64(axis)):
            continue
        for p in proof.levels[i]:
            _grow(values, budget)
            comptime if k == 4:
                var specialized = _aa._specialize(p, result._data[].embedding, budget)
                if not specialized.is_zero():
                    values.append(sparse.Array[P, k](rebind_var[SparsePoly[_Ring[k]]](specialized^), budget))
            else:
                values.append(sparse.Array[P, k](rebind_var[SparsePoly[_Ring[k]]](_sparse_copy(p, _variables(sparse._width[P](), List[Int64](), budget), False, budget)), budget))
    return values^


def projection_polynomials[P: Copyable & Deinitable, k: Int](result: Projection[P, k], *, axis: Int = -1, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[sparse.Array[P, k]]:
    """Copy projection basis polynomials for one axis or the complete closure.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.
        k: The inferred real coefficient field.

    Args:
        result: The result retaining its original inputs and evidence.
        axis: The original polynomial axis; -1 collects every level in CAD order.
        max_steps: The complete call’s work allowance.

    Returns:
        Independent sparse arrays retaining all original axes.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return projection_polynomials(result, budget, axis=axis)


def sample[P: Copyable & Deinitable](result: CellAt[P], mut budget: Budget) raises PolyError -> List[Real]:
    """Copy the retained query point in original polynomial axis order.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Exact coordinates of the point whose cell was requested.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    _allocation[Real](len(data.point), budget)
    var output = List[Real](capacity=len(data.point))
    for axis in range(sparse._width[P]()):
        for i in range(len(data.variables)):
            budget.charge(1)
            if data.variables[i] == Int64(axis):
                output.append(_plain(_copy_witness(data.point[i], budget), budget))
    return output^


def sample[P: Copyable & Deinitable](result: CellAt[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Real]:
    """Copy the retained query point in original polynomial axis order.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Exact coordinates of the point whose cell was requested.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return sample(result, budget)


def signs[P: Copyable & Deinitable](result: CellAt[P], mut budget: Budget) raises PolyError -> List[Int]:
    """Return the constant input-polynomial signs on the selected cell.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Minus one, zero or one for each input polynomial in input order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref proof = result._data[].proof
    ref embedding = result._data[].embedding
    var count = len(embedding.value()[].originals) if embedding else len(proof.signs)
    if count > len(proof.signs):
        raise PolyError(PolyErrorKind.DOMAIN, "CAD cell has too few input signs")
    _allocation[Int](count, budget)
    var values = List[Int](capacity=count)
    for i in range(count):
        budget.charge(1)
        values.append(proof.signs[i].relation)
    return values^


def signs[P: Copyable & Deinitable](result: CellAt[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Int]:
    """Return the constant input-polynomial signs on the selected cell.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Minus one, zero or one for each input polynomial in input order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return signs(result, budget)


def vanishing_orders[P: Copyable & Deinitable](result: Lazard[P], mut budget: Budget) raises PolyError -> List[Tuple[Int, Int]]:
    """Copy the minimal vanishing order for each assigned variable.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Axis and vanishing-order pairs in the original assignment order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    if len(data.point) != len(data.proof.orders):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid Lazard order count")
    _allocation[Tuple[Int, Int]](len(data.point), budget)
    var values = List[Tuple[Int, Int]](capacity=len(data.point))
    for i in range(len(data.point)):
        budget.charge(1)
        if data.point[i][0] < Int64(sparse._width[P]()):
            values.append((Int(data.point[i][0]), data.proof.orders[i]))
    return values^


def vanishing_orders[P: Copyable & Deinitable](result: Lazard[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Int, Int]]:
    """Copy the minimal vanishing order for each assigned variable.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Axis and vanishing-order pairs in the original assignment order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return vanishing_orders(result, budget)


def lifting_coefficients[P: Copyable & Deinitable](result: Lazard[P], mut budget: Budget) raises PolyError -> List[Real]:
    """Copy ascending exact-real coefficients of the Lazard lifting polynomial.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Exact coefficients usable directly with dense.real_roots; verify the Lazard result with cad.check.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref data = result._data[]
    if not len(data.proof.evaluations):
        raise PolyError(PolyErrorKind.DOMAIN, "missing Lazard specialization")
    ref proof = data.proof.evaluations[len(data.proof.evaluations) - 1]
    if len(proof.monomials) != len(proof.evaluations):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid specialization coefficient count")
    var degree = -1
    for powers in proof.monomials:
        budget.charge(1)
        var exponent = 0
        for pair in powers.pairs:
            budget.charge(1)
            if pair[0] != data.axis:
                raise PolyError(PolyErrorKind.DOMAIN, "lifting polynomial has another free axis")
            exponent = Int(pair[1])
        degree = max(degree, exponent)
    _allocation[Real](degree + 1, budget)
    var zero = _alg.asreal(0, budget)
    var values = List[Real](length=degree + 1, fill=zero)
    for i in range(len(proof.monomials)):
        budget.charge(1)
        var exponent = 0
        if len(proof.monomials[i].pairs):
            exponent = Int(proof.monomials[i].pairs[0][1])
        values[exponent] = _plain(_copy_witness(proof.evaluations[i].value, budget), budget)
    return values^


def lifting_coefficients[P: Copyable & Deinitable](result: Lazard[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Real]:
    """Copy ascending exact-real coefficients of the Lazard lifting polynomial.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Exact coefficients usable directly with dense.real_roots; verify the Lazard result with cad.check.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return lifting_coefficients(result, budget)


def cover[I: Iterable, P: Copyable & Deinitable](constraints: I, mut budget: Budget, *, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Covering[P] where _Element[I] == Constraint[P]:
    """Decide sign constraints using a complete CAD decomposition.

    The result retains the full decomposition and a witness or conflict covering. Use solve when a selective search is sufficient.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        budget: The shared work allowance for preparation and computation.
        axes: All axes in base-to-top order; defaults to increasing axis order.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    return cover[I, P, List[Tuple[Int, Int]]](constraints, budget, axes=axes, lazard=lazard, fixed=List[Tuple[Int, Int]]())


def cover[I: Iterable, P: Copyable & Deinitable](constraints: I, *, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Covering[P] where _Element[I] == Constraint[P]:
    """Decide sign constraints using a complete CAD decomposition.

    The result retains the full decomposition and a witness or conflict covering. Use solve when a selective search is sufficient.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order; defaults to increasing axis order.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return cover[I, P](constraints, budget, axes=axes, lazard=lazard)


def sample[P: Copyable & Deinitable](result: Covering[P], mut budget: Budget) raises PolyError -> Optional[List[Real]]:
    """Copy a satisfying sample from a full-decomposition covering.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Exact coordinates in original axis order, or None for an infeasible conjunction.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    ref proof = result._data[].proof
    var point = proof.sample(budget)
    if not point:
        return None
    var variables = proof.variables(budget)
    _allocation[Int](sparse._width[P](), budget)
    var indices = List[Int](length=sparse._width[P](), fill=-1)
    for k in range(len(variables)):
        budget.charge(1)
        if variables[k] >= 0 and variables[k] < Int64(len(indices)):
            indices[Int(variables[k])] = k
    _allocation[Real](len(indices), budget)
    var output = List[Real](capacity=len(indices))
    for index in indices:
        if index < 0 or index >= len(point.value()):
            raise PolyError(PolyErrorKind.DOMAIN, "incomplete covering sample")
        output.append(_output(point.value()[index], budget))
    return Optional(output^)


def sample[P: Copyable & Deinitable](result: Covering[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Optional[List[Real]]:
    """Copy a satisfying sample from a full-decomposition covering.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Exact coordinates in original axis order, or None for an infeasible conjunction.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return sample(result, budget)


def _describe_cells[P: Copyable & Deinitable](proof: CADFinalCertificate, embedding: Optional[ArcPointer[_aa._Embedding]], fixed: List[Tuple[Int, Real]], mut budget: Budget) raises PolyError -> List[Cell]:
    var values = List[Cell]()
    var stored = proof.cells(budget)
    var variables = proof.variables(budget)
    # Each distinct stored sample is refined to the output precision once.
    _allocation[Optional[Real]](len(stored), budget)
    var outputs = List[Optional[Real]](length=len(stored), fill=None)
    var refined = Dict[String, Int]()
    var level = len(variables) - 1
    for i in range(len(stored)):
        budget.charge(1)
        ref cell = stored[i]
        if cell.level != level:
            continue
        if embedding:
            ref bound = embedding.value()[]
            if len(cell.signs) != len(bound.originals) + 3 * len(bound.roots) + len(fixed):
                raise PolyError(PolyErrorKind.DOMAIN, "CAD cell has the wrong sign count")
            var on_fiber = True
            for k in range(len(bound.roots)):
                budget.charge(3)
                var offset = len(bound.originals) + 3 * k
                on_fiber = on_fiber and cell.signs[offset] == 0 and cell.signs[offset + 1] == 1 and cell.signs[offset + 2] == 1
            if not on_fiber:
                continue
        var dimension = 0
        var cursor = i
        while cursor >= 0 and stored[cursor].level >= 0:
            budget.charge(1)
            if stored[cursor].position % 2 == 0 and variables[stored[cursor].level] < Int64(sparse._width[P]()):
                dimension += 1
            cursor = stored[cursor].parent
        var count = len(embedding.value()[].originals) if embedding else len(cell.signs) - len(fixed)
        _allocation[Int](count, budget)
        var signs = List[Int](capacity=count)
        for j in range(count):
            budget.charge(1)
            signs.append(cell.signs[j])
        _grow(values, budget)
        values.append(Cell(i, dimension, _point_cached[P](stored, variables, i, outputs, refined, budget), signs^))
    return values^


def cells[P: Copyable & Deinitable](result: Covering[P], mut budget: Budget) raises PolyError -> List[Cell]:
    """Copy every terminal cell from a full-decomposition covering.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Cells in certificate order, with dimensions, original-axis samples and input signs.

    Raises:
        PolyError: DOMAIN for incompatible input, ENCODING for invalid transport records where applicable, or LIMIT/CAPACITY on resource limits.
    """
    return _describe_cells[P](result._data[].proof, result._data[].embedding, result._data[].fixed, budget)


def cells[P: Copyable & Deinitable](result: Covering[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Cell]:
    """Copy every terminal cell from a full-decomposition covering.

    Parameters:
        P: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Cells in certificate order, with dimensions, original-axis samples and input signs.

    Raises:
        PolyError: DOMAIN for incompatible input, ENCODING for invalid transport records where applicable, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return cells(result, budget)



def decompose[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, F: AnyType](polynomials: I, mut budget: Budget, *, fixed: F, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Decomposition[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Decompose the requested coordinate slice into sign-invariant cells.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        budget: The shared work allowance for preparation and computation.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        axes: All axes in base-to-top order, fixed axes first; empty puts increasing fixed axes first.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _inputs[I, P, C](polynomials, budget)
    var bindings = _cf._collect(fixed, sparse._width[P](), budget)
    var original_axes = _cf._order[P](bindings, axes, budget)
    var inputs = List[SparsePoly[QQ]]()
    swap(inputs, prepared.inputs)
    _cf._augment(inputs, bindings, budget)
    var variables = _aa._axes(original_axes, prepared.embedding, budget)
    if not _cf._check(inputs, prepared.embedding, bindings, variables, sparse._width[P](), budget):
        raise PolyError(PolyErrorKind.DOMAIN, "CAD input binding failed verification")
    var proof = cad_final(inputs, variables, budget, fixed=_cf._prefix(prepared.embedding, bindings, budget), lazard=lazard)
    _allocation[_DecompositionData](1, budget)
    return Decomposition[P](ArcPointer(_DecompositionData(inputs^, proof^, prepared.embedding.copy(), bindings^)))


def decompose[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, F: AnyType](polynomials: I, *, fixed: F, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Decomposition[P] where _Element[_Nested[I]] == Tuple[P, C]:
    """Decompose the requested coordinate slice into sign-invariant cells.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        polynomials: A finite iterable of integer, rational, or real-algebraic sparse polynomials.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order, fixed axes first; empty puts increasing fixed axes first.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return decompose[I, P, C, F](polynomials, budget, axes=axes, lazard=lazard, fixed=fixed)


def solve[I: Iterable, P: Copyable & Deinitable, F: AnyType](constraints: I, mut budget: Budget, *, fixed: F, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> Solution[P] where _Element[I] == Constraint[P]:
    """Find a satisfying point on the requested slice, or prove none exists there.

    The constraints form a conjunction. Budget exhaustion raises and never
    masquerades as an unsatisfiable result.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        budget: The shared work allowance for preparation and computation.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        axes: All axes in base-to-top order, fixed axes first; empty puts increasing fixed axes first.
        lazard: Use Lazard projection; False selects McCallum projection.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _constraints[I, P](constraints, budget)
    var bindings = _cf._collect(fixed, sparse._width[P](), budget)
    var original_axes = _cf._order[P](bindings, axes, budget)
    var inputs = List[SparsePoly[QQ]]()
    var masks = List[UInt8]()
    swap(inputs, prepared.problem.inputs)
    swap(masks, prepared.masks)
    _cf._augment(inputs, bindings, budget)
    _cf._masks(masks, bindings, budget)
    var variables = _aa._axes(original_axes, prepared.problem.embedding, budget)
    var proof = _cf._fiber(prepared.problem.embedding, bindings, budget).search(inputs, masks, variables, budget, lazard=lazard)
    _allocation[_SolutionData](1, budget)
    return Solution[P](ArcPointer(_SolutionData(inputs^, masks^, proof^, prepared.problem.embedding.copy(), bindings^)))


def solve[I: Iterable, P: Copyable & Deinitable, F: AnyType](constraints: I, *, fixed: F, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> Solution[P] where _Element[I] == Constraint[P]:
    """Find a satisfying point on the requested slice, or prove none exists there.

    The constraints form a conjunction. Budget exhaustion raises and never
    masquerades as an unsatisfiable result.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order, fixed axes first; empty puts increasing fixed axes first.
        lazard: Use Lazard projection; False selects McCallum projection.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return solve[I, P, F](constraints, budget, axes=axes, lazard=lazard, fixed=fixed)


def cover[I: Iterable, P: Copyable & Deinitable, F: AnyType](constraints: I, mut budget: Budget, *, fixed: F, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Covering[P] where _Element[I] == Constraint[P]:
    """Decide sign constraints using a complete decomposition of the requested slice.

    The result retains the full decomposition and a witness or conflict covering. Use solve when a selective search is sufficient.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        budget: The shared work allowance for preparation and computation.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        axes: All axes in base-to-top order, fixed axes first; empty puts increasing fixed axes first.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _constraints[I, P](constraints, budget)
    var bindings = _cf._collect(fixed, sparse._width[P](), budget)
    var original_axes = _cf._order[P](bindings, axes, budget)
    var inputs = List[SparsePoly[QQ]]()
    var masks = List[UInt8]()
    swap(inputs, prepared.problem.inputs)
    swap(masks, prepared.masks)
    _cf._augment(inputs, bindings, budget)
    _cf._masks(masks, bindings, budget)
    var variables = _aa._axes(original_axes, prepared.problem.embedding, budget)
    if not _cf._check(inputs, prepared.problem.embedding, bindings, variables, sparse._width[P](), budget) or not _cf._check_masks(inputs, masks, prepared.problem.embedding, bindings, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "CAD covering input binding failed verification")
    var proof = cad_final(inputs, variables, budget, fixed=_cf._prefix(prepared.problem.embedding, bindings, budget), masks=masks, query=True, lazard=lazard)
    _allocation[_CoveringData](1, budget)
    return Covering[P](ArcPointer(_CoveringData(inputs^, masks^, proof^, prepared.problem.embedding.copy(), bindings^)))


def cover[I: Iterable, P: Copyable & Deinitable, F: AnyType](constraints: I, *, fixed: F, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Optional[Bool] = None) raises PolyError -> Covering[P] where _Element[I] == Constraint[P]:
    """Decide sign constraints using a complete decomposition of the requested slice.

    The result retains the full decomposition and a witness or conflict covering. Use solve when a selective search is sufficient.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        constraints: A finite iterable of constraints sharing one power shape.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order, fixed axes first; empty puts increasing fixed axes first.
        lazard: None selects Revision 15 with interval lifting and complete final verification; True requests plain Lazard, False uses the checked legacy McCallum backend.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return cover[I, P, F](constraints, budget, axes=axes, lazard=lazard, fixed=fixed)


def decompose_to[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, S: CADSink, F: AnyType](polynomials: I, mut sink: S, mut budget: Budget, *, fixed: F, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> CADStreamReceipt where _Element[_Nested[I]] == Tuple[P, C]:
    """Stream a decomposition of the requested slice without retaining its full history.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        S: The inferred input type.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        polynomials: The integer, rational, or real-algebraic sparse polynomial inputs.
        sink: The consumer of projection and depth-first witness events.
        budget: The shared work allowance for preparation and computation.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        axes: All axes in base-to-top order, fixed axes first; empty chooses their order.
        lazard: Use Lazard projection by default.

    Returns:
        A completion receipt. Events written before an error do not prove a complete decomposition.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _inputs[I, P, C](polynomials, budget)
    var bindings = _cf._collect(fixed, sparse._width[P](), budget)
    var original_axes = _cf._order[P](bindings, axes, budget)
    var inputs = List[SparsePoly[QQ]]()
    swap(inputs, prepared.inputs)
    _cf._augment(inputs, bindings, budget)
    var variables = _aa._axes(original_axes, prepared.embedding, budget)
    return _cf._fiber(prepared.embedding, bindings, budget).decompose_to(inputs, variables, sink, budget, lazard=lazard)


def decompose_to[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, S: CADSink, F: AnyType](polynomials: I, mut sink: S, *, fixed: F, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int](), lazard: Bool = True) raises PolyError -> CADStreamReceipt where _Element[_Nested[I]] == Tuple[P, C]:
    """Stream a decomposition of the requested slice without retaining its full history.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        S: The inferred input type.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        polynomials: The integer, rational, or real-algebraic sparse polynomial inputs.
        sink: The consumer of projection and depth-first witness events.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: All axes in base-to-top order, fixed axes first; empty chooses their order.
        lazard: Use Lazard projection by default.

    Returns:
        A completion receipt. Events written before an error do not prove a complete decomposition.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return decompose_to[I, P, C, S, F](polynomials, sink, budget, axes=axes, lazard=lazard, fixed=fixed)


def checker[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, F: AnyType](polynomials: I, mut budget: Budget, *, fixed: F, axes: List[Int] = List[Int]()) raises PolyError -> _cf._Checker where _Element[_Nested[I]] == Tuple[P, C]:
    """Prepare independent incremental checking for the requested coordinate slice.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        polynomials: The original input polynomials.
        budget: The shared work allowance for preparation and computation.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        axes: The producer's complete traversal order; empty uses the same default.

    Returns:
        A checker accepting streamed events; finish(receipt, budget) verifies complete traversal.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var prepared = _inputs[I, P, C](polynomials, budget)
    var bindings = _cf._collect(fixed, sparse._width[P](), budget)
    var original_axes = _cf._order[P](bindings, axes, budget)
    var variables = _aa._axes(original_axes, prepared.embedding, budget)
    _cf._augment(prepared.inputs, bindings, budget)
    var core = _cf._fiber(prepared.embedding, bindings, budget).checker(prepared.inputs, budget)
    return _cf._Checker(core^, variables^, False)


def checker[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, F: AnyType](polynomials: I, *, fixed: F, axes: List[Int] = List[Int](), max_steps: UInt64 = UInt64.MAX) raises PolyError -> _cf._Checker where _Element[_Nested[I]] == Tuple[P, C]:
    """Prepare independent incremental checking for the requested coordinate slice.

    Parameters:
        I: The inferred input type.
        P: The inferred power shape.
        C: The inferred input type.
        F: The inferred iterable or tuple of fixed axis/value pairs.

    Args:
        polynomials: The original input polynomials.
        fixed: Distinct axis/value pairs with integer, rational, or exact real-algebraic values.
        axes: The producer's complete traversal order; empty uses the same default.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A checker accepting streamed events; finish(receipt, budget) verifies complete traversal.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return checker[I, P, C, F](polynomials, budget, fixed=fixed, axes=axes)


def fixed_coordinates[P: Copyable & Deinitable](result: Decomposition[P], mut budget: Budget) raises PolyError -> List[Tuple[Int, Real]]:
    """Copy the user-fixed axes and exact values in traversal order.

    Parameters:
        P: The inferred polynomial power shape.

    Args:
        result: The retained CAD result.
        budget: The shared work allowance, including output storage.

    Returns:
        Independent axis/value pairs; empty when no user coordinates were fixed.

    Raises:
        PolyError: LIMIT or CAPACITY when resource limits are reached.
    """
    ref fixed = result._data[].fixed
    _allocation[Tuple[Int, Real]](len(fixed), budget)
    var output = List[Tuple[Int, Real]](capacity=len(fixed))
    for entry in fixed:
        budget.charge(1)
        output.append((entry[0], entry[1].copy()))
    return output^


def fixed_coordinates[P: Copyable & Deinitable](result: Decomposition[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Int, Real]]:
    """Copy the user-fixed axes and exact values in traversal order.

    Parameters:
        P: The inferred polynomial power shape.

    Args:
        result: The retained CAD result.
        max_steps: The complete call's work allowance.

    Returns:
        Independent axis/value pairs; empty when no user coordinates were fixed.

    Raises:
        PolyError: LIMIT or CAPACITY when resource limits are reached.
    """
    var budget = Budget(max_steps)
    return fixed_coordinates(result, budget)


def fixed_coordinates[P: Copyable & Deinitable](result: Solution[P], mut budget: Budget) raises PolyError -> List[Tuple[Int, Real]]:
    """Copy the user-fixed axes and exact values in traversal order.

    Parameters:
        P: The inferred polynomial power shape.

    Args:
        result: The retained CAD result.
        budget: The shared work allowance, including output storage.

    Returns:
        Independent axis/value pairs; empty when no user coordinates were fixed.

    Raises:
        PolyError: LIMIT or CAPACITY when resource limits are reached.
    """
    ref fixed = result._data[].fixed
    _allocation[Tuple[Int, Real]](len(fixed), budget)
    var output = List[Tuple[Int, Real]](capacity=len(fixed))
    for entry in fixed:
        budget.charge(1)
        output.append((entry[0], entry[1].copy()))
    return output^


def fixed_coordinates[P: Copyable & Deinitable](result: Solution[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Int, Real]]:
    """Copy the user-fixed axes and exact values in traversal order.

    Parameters:
        P: The inferred polynomial power shape.

    Args:
        result: The retained CAD result.
        max_steps: The complete call's work allowance.

    Returns:
        Independent axis/value pairs; empty when no user coordinates were fixed.

    Raises:
        PolyError: LIMIT or CAPACITY when resource limits are reached.
    """
    var budget = Budget(max_steps)
    return fixed_coordinates(result, budget)


def fixed_coordinates[P: Copyable & Deinitable](result: Covering[P], mut budget: Budget) raises PolyError -> List[Tuple[Int, Real]]:
    """Copy the user-fixed axes and exact values in traversal order.

    Parameters:
        P: The inferred polynomial power shape.

    Args:
        result: The retained CAD result.
        budget: The shared work allowance, including output storage.

    Returns:
        Independent axis/value pairs; empty when no user coordinates were fixed.

    Raises:
        PolyError: LIMIT or CAPACITY when resource limits are reached.
    """
    ref fixed = result._data[].fixed
    _allocation[Tuple[Int, Real]](len(fixed), budget)
    var output = List[Tuple[Int, Real]](capacity=len(fixed))
    for entry in fixed:
        budget.charge(1)
        output.append((entry[0], entry[1].copy()))
    return output^


def fixed_coordinates[P: Copyable & Deinitable](result: Covering[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Int, Real]]:
    """Copy the user-fixed axes and exact values in traversal order.

    Parameters:
        P: The inferred polynomial power shape.

    Args:
        result: The retained CAD result.
        max_steps: The complete call's work allowance.

    Returns:
        Independent axis/value pairs; empty when no user coordinates were fixed.

    Raises:
        PolyError: LIMIT or CAPACITY when resource limits are reached.
    """
    var budget = Budget(max_steps)
    return fixed_coordinates(result, budget)
