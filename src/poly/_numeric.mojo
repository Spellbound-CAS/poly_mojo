"""Numerical array storage and explicit-context scalar conversion."""
from ._io_meta import _Record, _Dense, _Sparse
from std.iter import Iterable, Iterator, StopIteration, iter, next
from std.memory import ArcPointer
from std.builtin.rebind import downcast
from apn_mojo import Float, Complex, Ball, ComplexBall, ArithmeticContext, ComplexContext, BallContext, FloatFormat, Rational, ExactComplex, Batch
from polynomials import Budget, PolyError, PolyErrorKind, QQ, NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, NumericBox, numeric_bounds, numeric_value, NumericDensePoly, NumericSparsePoly, NumericTerm, MonomialOrder
from ._common import _Kind, _Element, _rational, _size, _allocation, _charge_size, _grow, _compact, _check_batch

from ._algebraic_value import Real as _Real, Complex as _Complex

comptime _IsComplex[T: AnyType] = T == Complex or T == ComplexBall or T == ExactComplex or T == _Complex
comptime _Domain[C: AnyType, T: AnyType] = ComplexDomain if C == ComplexContext else FloatDomain if C == ArithmeticContext else ComplexBallDomain if _IsComplex[T] else BallDomain


def _context[D: NumericDomain](context: D.Context, mut budget: Budget) raises PolyError:
    comptime if D == BallDomain or D == ComplexBallDomain:
        if not rebind[BallContext](context).precision():
            raise PolyError(PolyErrorKind.DOMAIN, "choose an explicit ball precision")
        _charge_size(rebind[BallContext](context).precision().value(), budget)
    elif D == FloatDomain:
        _charge_size(rebind[ArithmeticContext](context).format().precision(), budget)
    else:
        var c = rebind[ComplexContext](context)
        _charge_size(c.real().format().precision(), budget)
        _charge_size(c.imag().format().precision(), budget)


def _bounds[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> NumericBox:
    comptime if T == _Real:
        from . import _algebraic_numeric
        return _algebraic_numeric._box(_algebraic_numeric._parts(rebind[_Real](value), budget), budget)
    elif T == _Complex:
        from . import _algebraic_numeric
        return _algebraic_numeric._box(_algebraic_numeric._parts(rebind[_Complex](value), budget), budget)
    elif T == Float:
        return numeric_bounds[FloatDomain](rebind[Float](value), budget)
    elif T == Complex:
        return numeric_bounds[ComplexDomain](rebind[Complex](value), budget)
    elif T == Ball:
        return numeric_bounds[BallDomain](rebind[Ball](value), budget)
    elif T == ComplexBall:
        return numeric_bounds[ComplexBallDomain](rebind[ComplexBall](value), budget)
    elif T == Float64 or T == Float32 or T == Float16 or T == BFloat16:
        budget.charge(2048)
        var v: Float
        try:
            v = Float(value, context=ArithmeticContext(format=FloatFormat(precision=64)))
        except e:
            raise PolyError(PolyErrorKind.DOMAIN, String(e))
        return numeric_bounds[FloatDomain](v, budget)
    elif T == ExactComplex:
        _charge_size(_size(value), budget)
        var z = rebind[ExactComplex](value)
        return NumericBox.point(QQ(z.real()), QQ(z.imag()), budget)
    else:
        comptime assert _Kind[T] == 0 or _Kind[T] == 1, "numeric coefficients require real/complex numbers, not modular values or Bool"
        _charge_size(_size(value), budget)
        return NumericBox.point(QQ(_rational(value)), QQ(Rational(0)), budget)


def _value[D: NumericDomain, T: Copyable & Deinitable](value: T, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    comptime if T == _Real:
        from . import _algebraic_numeric
        return _algebraic_numeric._value[D](rebind[_Real](value), context, budget)
    elif T == _Complex:
        from . import _algebraic_numeric
        return _algebraic_numeric._value[D](rebind[_Complex](value), context, budget)
    else:
        return numeric_value[D](_bounds(value, budget), context, budget)


struct DenseArray[D: NumericDomain](Copyable, Movable, Iterable, Sized, _Record):
    """An immutable numerical polynomial with shared canonical storage.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.
    """
    comptime _Format = _Dense
    comptime _Item = Self.D.Scalar
    comptime _Power = Int
    var _core: ArcPointer[NumericDensePoly[Self.D]]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _DenseArrayIterator[Self.D, ImmOrigin(iterable_origin)]
    """Borrow the canonical numerical entries.

    Parameters:
        iterable_mut: Whether the borrowing source origin allows mutation.
        iterable_origin: The lifetime of the borrowed source.

    Returns:
        The borrowing iterator type.
    """

    def __init__(out self, var core: NumericDensePoly[Self.D], mut budget: Budget) raises PolyError:
        """Take canonical storage; construct user values through asarray.

        Args:
            core: Owned canonical core storage.
            budget: The shared work allowance, including preparation, computation, and result storage.

        Raises:
            PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
        """
        _compact(core.coefficients, budget)
        _allocation[NumericDensePoly[Self.D]](1, budget)
        self._core = ArcPointer(core^)

    def __len__(self) -> Int:
        """Return the number of canonical stored entries.

        Returns:
            The number of logical entries.
        """
        return len(self._core[].coefficients)

    def __getitem__(self, index: Int) -> Self.D.Scalar:
        """Copy the entry at a valid zero-based index.

        Args:
            index: A valid zero-based index into the logical entries.

        Returns:
            An independent value copy; polynomial and evidence buffers may remain shared.
        """
        return self._core[].coefficients[index].copy()

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over numerical entries in canonical order.

        Returns:
            An iterator borrowing the immutable source storage.
        """
        return _DenseArrayIterator[Self.D, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def to_list(self, mut budget: Budget) raises PolyError -> List[Self.D.Scalar]:
        """Copy numerical entries into an independent mutable list.

        Args:
            budget: The shared work allowance, including preparation, computation, and result storage.

        Returns:
            An independent mutable list of entries.

        Raises:
            PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
        """
        _allocation[Self.D.Scalar](len(self), budget)
        var values = List[Self.D.Scalar](capacity=len(self))
        for i in range(len(self)):
            _ = _bounds(self._core[].coefficients[i], budget)
            values.append(self[i])
        return values^

    def to_list(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Self.D.Scalar]:
        """Materialize independent numerical entries with one work allowance.

        Args:
            max_steps: The complete call's work allowance; defaults to UInt64.MAX.

        Returns:
            An independent mutable list of entries.

        Raises:
            PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
        """
        var budget = Budget(max_steps)
        return self.to_list(budget)


@fieldwise_init
struct _DenseArrayIterator[D: NumericDomain, origin: ImmOrigin](Iterator):
    comptime Element = Self.D.Scalar
    var source: Pointer[DenseArray[Self.D], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.Element:
        if self.position == len(self.source[]):
            raise StopIteration()
        var value = self.source[][self.position]
        self.position += 1
        return value^


struct SparseArray[P: Copyable & Deinitable, D: NumericDomain](Copyable, Movable, Iterable, Sized, _Record):
    """An immutable numerical polynomial with shared canonical storage.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.
        D: The APN numerical domain associated with the stored scalars.
    """
    comptime _Format = _Sparse
    comptime _Item = Self.D.Scalar
    comptime _Power = Self.P
    var _core: ArcPointer[NumericSparsePoly[Self.D]]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _SparseArrayIterator[Self.P, Self.D, ImmOrigin(iterable_origin)]
    """Borrow the canonical numerical entries.

    Parameters:
        iterable_mut: Whether the borrowing source origin allows mutation.
        iterable_origin: The lifetime of the borrowed source.

    Returns:
        The borrowing iterator type.
    """

    def __init__(out self, var core: NumericSparsePoly[Self.D], mut budget: Budget) raises PolyError:
        """Take canonical storage; construct user values through asarray.

        Args:
            core: Owned canonical core storage.
            budget: The shared work allowance, including preparation, computation, and result storage.

        Raises:
            PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
        """
        _compact(core.terms, budget)
        _allocation[NumericSparsePoly[Self.D]](1, budget)
        self._core = ArcPointer(core^)

    def __len__(self) -> Int:
        """Return the number of canonical stored entries.

        Returns:
            The number of logical entries.
        """
        return len(self._core[].terms)

    def __getitem__(self, index: Int) -> Tuple[Self.P, Self.D.Scalar]:
        """Copy the entry at a valid zero-based index.

        Args:
            index: A valid zero-based index into the logical entries.

        Returns:
            An independent value copy; polynomial and evidence buffers may remain shared.
        """
        from .sparse import _powers
        ref term = self._core[].terms[index]
        return (_powers[Self.P](term.monomial), term.coefficient.copy())

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over numerical entries in canonical order.

        Returns:
            An iterator borrowing the immutable source storage.
        """
        return _SparseArrayIterator[Self.P, Self.D, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def to_list(self, mut budget: Budget) raises PolyError -> List[Tuple[Self.P, Self.D.Scalar]]:
        """Copy numerical entries into an independent mutable list.

        Args:
            budget: The shared work allowance, including preparation, computation, and result storage.

        Returns:
            An independent mutable list of entries.

        Raises:
            PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
        """
        _allocation[Tuple[Self.P, Self.D.Scalar]](len(self), budget)
        var values = List[Tuple[Self.P, Self.D.Scalar]](capacity=len(self))
        for i in range(len(self)):
            _ = _bounds(self._core[].terms[i].coefficient, budget)
            values.append(self[i])
        return values^

    def to_list(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[Self.P, Self.D.Scalar]]:
        """Materialize independent numerical entries with one work allowance.

        Args:
            max_steps: The complete call's work allowance; defaults to UInt64.MAX.

        Returns:
            An independent mutable list of entries.

        Raises:
            PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
        """
        var budget = Budget(max_steps)
        return self.to_list(budget)


@fieldwise_init
struct _SparseArrayIterator[P: Copyable & Deinitable, D: NumericDomain, origin: ImmOrigin](Iterator):
    comptime Element = Tuple[Self.P, Self.D.Scalar]
    var source: Pointer[SparseArray[Self.P, Self.D], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.Element:
        if self.position == len(self.source[]):
            raise StopIteration()
        var value = self.source[][self.position]
        self.position += 1
        return value^


def _dense[D: NumericDomain, I: Iterable](values: I, context: D.Context, mut budget: Budget) raises PolyError -> DenseArray[D]:
    _context[D](context, budget)
    _check_batch(values)
    comptime if I == DenseArray[D]:
        budget.charge(1)
        return rebind[DenseArray[D]](values).copy()
    else:
        comptime if I == Batch[Float]:
            if rebind[Batch[Float]](values).ndim() != 1:
                raise PolyError(PolyErrorKind.DOMAIN, "select one coefficient axis before constructing a polynomial")
        comptime if I == Batch[Complex]:
            if rebind[Batch[Complex]](values).ndim() != 1:
                raise PolyError(PolyErrorKind.DOMAIN, "select one coefficient axis before constructing a polynomial")
        comptime if I == Batch[Ball]:
            if rebind[Batch[Ball]](values).ndim() != 1:
                raise PolyError(PolyErrorKind.DOMAIN, "select one coefficient axis before constructing a polynomial")
        comptime if I == Batch[ComplexBall]:
            if rebind[Batch[ComplexBall]](values).ndim() != 1:
                raise PolyError(PolyErrorKind.DOMAIN, "select one coefficient axis before constructing a polynomial")
        var coefficients = List[D.Scalar]()
        var cursor = iter(values)
        while True:
            budget.charge(1)
            var value: downcast[_Element[I], Copyable & Deinitable]
            try:
                value = rebind_var[downcast[_Element[I], Copyable & Deinitable]](next(cursor))
            except StopIteration:
                break
            _grow(coefficients, budget)
            coefficients.append(_value[D](value, context, budget))
        return DenseArray[D](NumericDensePoly[D].from_coefficients(0, coefficients^, budget), budget)


def _sparse[D: NumericDomain, I: Iterable, P: Copyable & Deinitable, T: Copyable & Deinitable](values: I, context: D.Context, mut budget: Budget) raises PolyError -> SparseArray[P, D] where _Element[I] == Tuple[P, T]:
    from .sparse import _monomial, _width
    from ._interop import _variables
    _context[D](context, budget)
    comptime if I == SparseArray[P, D]:
        budget.charge(1)
        return rebind[SparseArray[P, D]](values).copy()
    else:
        var terms = List[NumericTerm[D]]()
        var cursor = iter(values)
        while True:
            budget.charge(1)
            var value: Tuple[P, T]
            try:
                value = rebind_var[Tuple[P, T]](next(cursor))
            except StopIteration:
                break
            _grow(terms, budget)
            terms.append(NumericTerm[D](_monomial(value[0], budget), _value[D](value[1], context, budget)))
        var order = MonomialOrder(MonomialOrder.LEX, _variables(_width[P](), List[Int64](), budget), budget)
        return SparseArray[P, D](NumericSparsePoly[D].from_terms(terms^, order, context, budget), budget)


comptime _StoredDomain[T: AnyType] = FloatDomain if T == Float else ComplexDomain if T == Complex else BallDomain if T == Ball else ComplexBallDomain


def _validate_sparse_axes[D: NumericDomain](p: NumericSparsePoly[D], width: Int, mut budget: Budget) raises PolyError:
    if width < 0 or len(p.order.variables) != width or p.order.kind != MonomialOrder.LEX:
        raise PolyError(PolyErrorKind.DOMAIN, "serialized sparse shape must match ndim and use lexicographic array axes")
    for i in range(width):
        budget.charge(1)
        if p.order.variables[i] != Int64(i):
            raise PolyError(PolyErrorKind.DOMAIN, "serialized numerical axes must be consecutive from zero")


def _copy_dense[D: NumericDomain](p: NumericDensePoly[D], mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    _allocation[D.Scalar](len(p.coefficients), budget)
    for value in p.coefficients:
        _ = numeric_bounds[D](value, budget)
    return NumericDensePoly[D].from_coefficients(0, p.coefficients.copy(), budget)


def _copy_sparse[D: NumericDomain](p: NumericSparsePoly[D], width: Int, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    _validate_sparse_axes(p, width, budget)
    _allocation[Int64](len(p.order.variables), budget)
    _allocation[NumericTerm[D]](len(p.terms), budget)
    for term in p.terms:
        _allocation[Tuple[Int64, UInt32]](len(term.monomial.pairs), budget)
        _ = numeric_bounds[D](term.coefficient, budget)
    return p.copy()


comptime _AnyComplex[previous: Bool, T: Movable] = previous or _IsComplex[T]
comptime _TupleScalar[T: type_of(Tuple)] = Complex if T.Ts.reduce[False, _AnyComplex] else Float


comptime _IsNumeric[T: AnyType] = T == Float or T == Complex or T == Ball or T == ComplexBall or T == Float64 or T == Float32 or T == Float16 or T == BFloat16


# Sequence preparation keeps interpolation constraints and evaluated zeros intact.
def _sequence[D: NumericDomain, I: Iterable](values: I, context: D.Context, mut budget: Budget) raises PolyError -> List[D.Scalar]:
    _context[D](context, budget)
    _check_sequence(values)
    var result = List[D.Scalar]()
    var cursor = iter(values)
    while True:
        budget.charge(1)
        var value: downcast[_Element[I], Copyable & Deinitable]
        try:
            value = rebind_var[downcast[_Element[I], Copyable & Deinitable]](next(cursor))
        except StopIteration:
            break
        _grow(result, budget)
        result.append(_value[D](value, context, budget))
    _compact(result, budget)
    return result^


struct Values[D: NumericDomain](Copyable, Movable, Iterable, Sized):
    """An immutable sequence from numerical samples, retaining every zero.

    Parameters:
        D: The numerical scalar domain.
    """
    var _values: ArcPointer[List[Self.D.Scalar]]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _ValuesIterator[Self.D, ImmOrigin(iterable_origin)]
    """A borrowing iterator over the shared storage.

    Parameters:
        iterable_mut: Whether the borrowing origin is mutable.
        iterable_origin: The lifetime of the borrowed array.
    """

    def __init__(out self, var values: List[Self.D.Scalar], mut budget: Budget) raises PolyError:
        """Take ownership of values and their already validated numerical domain.

        Args:
            values: The numerical values in their original order.
            budget: The caller-owned work allowance, shared through preparation and computation.

        Raises:
            PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
        """
        _compact(values, budget)
        _allocation[List[Self.D.Scalar]](1, budget)
        self._values = ArcPointer(values^)

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over values without copying the backing buffer.

        Returns:
            An iterator borrowing the shared storage.
        """
        return _ValuesIterator[Self.D, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def __len__(self) -> Int:
        """Return the number of stored values.

        Returns:
            The number of stored entries.
        """
        return len(self._values[])

    def __getitem__(self, index: Int) -> Self.D.Scalar:
        """Return the entry at a valid zero-based index.

        Args:
            index: A valid zero-based index into stored entries.

        Returns:
            An independent numerical scalar.
        """
        return self._values[][index].copy()

    def to_list(self, mut budget: Budget) raises PolyError -> List[Self.D.Scalar]:
        """Return an independent mutable list, charging allocation and coefficient work.

        Args:
            budget: The caller-owned work allowance, shared through preparation and computation.

        Returns:
            An independent mutable list in the current entry order.

        Raises:
            PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
        """
        _allocation[Self.D.Scalar](len(self), budget)
        var result = List[Self.D.Scalar](capacity=len(self))
        for index in range(len(self)):
            _ = numeric_bounds[Self.D](self._values[][index], budget)
            result.append(self[index])
        return result^

    def to_list(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Self.D.Scalar]:
        """Return an independent mutable list, charging allocation and coefficient work.

        Args:
            max_steps: The complete call's work allowance; defaults to UInt64.MAX.

        Returns:
            An independent mutable list in the current entry order.

        Raises:
            PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
        """
        var budget = Budget(max_steps)
        return self.to_list(budget)

@fieldwise_init
struct _ValuesIterator[D: NumericDomain, origin: ImmOrigin](Iterator):
    comptime Element = Self.D.Scalar
    var source: Pointer[Values[Self.D], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.D.Scalar:
        if self.position == len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position]
        self.position += 1
        return result^

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        var remaining = len(self.source[]) - self.position
        return remaining, Optional(remaining)



comptime _CombinedScalar[A: AnyType, B: AnyType] = Complex if _IsComplex[A] or _IsComplex[B] else Float


def _check_sequence[I: Iterable](values: I) raises PolyError:
    _check_batch(values)
    comptime if I == Batch[Float]:
        if rebind[Batch[Float]](values).ndim() != 1:
            raise PolyError(PolyErrorKind.DOMAIN, "select one sequence axis")
    comptime if I == Batch[Complex]:
        if rebind[Batch[Complex]](values).ndim() != 1:
            raise PolyError(PolyErrorKind.DOMAIN, "select one sequence axis")
    comptime if I == Batch[Ball]:
        if rebind[Batch[Ball]](values).ndim() != 1:
            raise PolyError(PolyErrorKind.DOMAIN, "select one sequence axis")
    comptime if I == Batch[ComplexBall]:
        if rebind[Batch[ComplexBall]](values).ndim() != 1:
            raise PolyError(PolyErrorKind.DOMAIN, "select one sequence axis")
