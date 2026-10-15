"""Exact promotion, checked allocations, and shared work accounting."""
from std.iter import Iterable, IterableOwned, Iterator, StopIteration, iter, next
from std.builtin.rebind import downcast
from std.memory import ArcPointer
from std.collections import Array
from apn_mojo import Integer, Rational, ExactComplex, Batch
from polynomials import ZZ, QQ, QQi, ModInt, Modulus, NoContext, real_rational, Budget, PolyError, PolyErrorKind
from std.sys import size_of
from ._algebraic_value import Real, Complex, _plain

comptime _Native[T: AnyType] = T == Int or T == Int8 or T == UInt8 or T == Int16 or T == UInt16 or T == Int32 or T == UInt32 or T == Int64 or T == UInt64
comptime _Kind[T: AnyType] = 5 if T == Complex else 4 if T == Real else 3 if T == ModInt else 2 if T == ExactComplex else 1 if T == Rational else 0 if T == Integer or _Native[T] else -1
comptime _CoreKind[C: AnyType] = 5 if C == Complex else 4 if C == Real else 0 if C == ZZ else 1 if C == QQ else 2 if C == QQi else 3 if C == ModInt else -1
comptime _Scalar[k: Int] = Complex if k == 5 else Real if k == 4 else ModInt if k == 3 else ExactComplex if k == 2 else Rational if k == 1 else Integer
comptime _Ring[k: Int] = Complex if k == 5 else Real if k == 4 else ModInt if k == 3 else QQi if k == 2 else QQ if k == 1 else ZZ
comptime _Element[I: Iterable] = I.IteratorType[ImmutAnyOrigin].Element


trait _FixedInput:
    @staticmethod
    def _coefficient_input_length() -> Int: ...


__extension Array(_FixedInput):
    @staticmethod
    def _coefficient_input_length() -> Int:
        return Self.length


comptime _DirectInput[I: Iterable] = I == List[_Element[I]] or conforms_to(I, _FixedInput)


def _input_view[I: Iterable](ref values: I) -> Span[downcast[_Element[I], Copyable & Deinitable], origin_of(values)]:
    comptime T = downcast[_Element[I], Copyable & Deinitable]
    comptime if I == List[T]:
        return Span(rebind[List[T]](values))
    else:
        comptime n = downcast[I, _FixedInput]._coefficient_input_length()
        return Span(rebind[Array[T, n]](values))


comptime _FieldKind[k: Int] = 1 if k == 0 else k
comptime _Join[a: Int, b: Int] = 3 if a == 3 or b == 3 else 5 if (a == 4 and b == 2) or (a == 2 and b == 4) else a if a >= b else b
comptime _JoinRank[previous: Int, T: Movable] = _Join[previous, _Kind[T]]
comptime _TupleKind[T: type_of(Tuple)] = T.Ts.reduce[0, _JoinRank]
comptime _OwnedElement[I: IterableOwned & Deinitable] = downcast[I.IteratorOwnedType.Element, Copyable & Deinitable]

trait _WithModulus:
    def modulus(self) -> UInt64: ...


def _domain[k: Int](modulus: UInt64) -> String:
    comptime if k == 5:
        return "QQbar"
    elif k == 4:
        return "AA"
    elif k == 3:
        return String("GF(", modulus, ")")
    elif k == 2:
        return "QQi"
    elif k == 1:
        return "QQ"
    else:
        return "ZZ"


def _initial_modulus[I: Iterable](values: I) -> UInt64:
    comptime if conforms_to(I, _WithModulus):
        return values.modulus()
    return 0


def _merge_modulus(a: UInt64, b: UInt64) raises PolyError -> UInt64:
    if a and b and a != b:
        raise PolyError(PolyErrorKind.DOMAIN, "coefficients have different moduli")
    return a if a else b


struct _Context[k: Int](Movable):
    """Validate an inferred field once, before converting any coefficient."""
    var _value: Optional[_Ring[Self.k].Context]

    def __init__(out self):
        """Take ownership of values and their already validated field metadata."""
        self._value = None
        comptime if Self.k != 3:
            self._value = rebind_var[_Ring[Self.k].Context](NoContext())

    def observe(mut self, modulus: UInt64, mut budget: Budget) raises PolyError:
        comptime if Self.k == 3:
            if self._value:
                _ = _merge_modulus(rebind[Modulus](self._value.value()).value, modulus)
            else:
                self._value = rebind_var[_Ring[Self.k].Context](Modulus(modulus, budget))

    def get(self) raises PolyError -> _Ring[Self.k].Context:
        if not self._value:
            raise PolyError(PolyErrorKind.DOMAIN, "empty modular input needs modulus= or a result that retains its modulus")
        return self._value.value().copy()


def _grow[T: Copyable & Deinitable](mut values: List[T], mut budget: Budget) raises PolyError:
    """Charge checked geometric growth before List can reallocate."""
    if len(values) == values.capacity():
        var capacity = max(1, _product(values.capacity(), 2))
        _allocation[T](capacity, budget)
        _charge_size(len(values), budget)
        values.reserve(capacity)


def _compact[T: Copyable & Deinitable](mut values: List[T], mut budget: Budget) raises PolyError:
    """Release empty or mostly unused buffers while moving, never copying, entries."""
    var count = len(values)
    if count == 0:
        values = List[T]()
    elif values.capacity() - count > count:
        _allocation[T](count, budget)
        _charge_size(_product(count, 2), budget)
        var compacted = List[T](capacity=count)
        # Pop moves owned entries out of the old buffer. Reverse once to restore
        # coefficient/term order without copying arbitrary-precision payloads.
        while len(values):
            compacted.append(values.pop())
        compacted.reverse()
        values = compacted^


def _check_batch[I: AnyType](values: I, *, message: String = "dense input must have one coefficient axis; select a row or flatten explicitly") raises PolyError:
    comptime if I == Batch[Integer]:
        if rebind[Batch[Integer]](values).ndim() != 1:
            raise PolyError(PolyErrorKind.DOMAIN, message)
    comptime if I == Batch[Rational]:
        if rebind[Batch[Rational]](values).ndim() != 1:
            raise PolyError(PolyErrorKind.DOMAIN, message)


def _size[T: Movable](value: T) raises PolyError -> Int:
    comptime if T == Real or T == Complex:
        return rebind[Real if T == Real else Complex](value).size()
    elif T == Integer:
        return _sum(rebind[Integer](value).magnitude_bit_length() // 8, 1)
    elif T == Rational:
        ref q = rebind[Rational](value)
        return _sum(_size(q.numerator()), _size(q.denominator()))
    elif T == ExactComplex:
        ref z = rebind[ExactComplex](value)
        return _sum(_size(z.real()), _size(z.imag()))
    else:
        return 8


def _integer[T: Movable](value: T) -> Integer:
    comptime assert _Kind[T] == 0, "expected an exact integer; floats and Bool require explicit conversion"
    comptime if T == Integer:
        return rebind[Integer](value)
    elif T == Int:
        return Integer(rebind[Int](value))
    elif T == Int8:
        return Integer(rebind[Int8](value))
    elif T == UInt8:
        return Integer(rebind[UInt8](value))
    elif T == Int16:
        return Integer(rebind[Int16](value))
    elif T == UInt16:
        return Integer(rebind[UInt16](value))
    elif T == Int32:
        return Integer(rebind[Int32](value))
    elif T == UInt32:
        return Integer(rebind[UInt32](value))
    elif T == Int64:
        return Integer(rebind[Int64](value))
    else:
        return Integer(rebind[UInt64](value))


def _rational[T: Movable](value: T) -> Rational:
    comptime if T == Rational:
        return rebind[Rational](value)
    else:
        return Rational(_integer(value))


def _convert[k: Int, T: Movable](value: T, ctx: _Ring[k].Context, mut budget: Budget) raises PolyError -> _Scalar[k]:
    comptime assert _Kind[T] >= 0, "exact coefficients must be integers, Rational, ExactComplex or ModInt; convert floats explicitly"
    var size = _size(value)
    # APN calls are indivisible; precharge promotion and modular reduction.
    comptime if k != _Kind[T]:
        size = _sum(size, 8)
        size = _product(size, size)
    _charge_size(size, budget)
    comptime if k >= 4:
        comptime assert T != ModInt, "modular coefficients have no implicit algebraic conversion"
        comptime if k == 4:
            comptime assert _Kind[T] == 0 or _Kind[T] == 1 or T == Real, "real algebraic coefficients require real exact values"
            comptime if T == Real:
                return rebind_var[_Scalar[k]](rebind[Real](value).copy())
            else:
                return rebind_var[_Scalar[k]](_plain(real_rational(QQ(_rational(value)), budget), budget))
        else:
            comptime if T == Complex:
                return rebind_var[_Scalar[k]](rebind[Complex](value).copy())
            elif T == ExactComplex:
                ref z = rebind[ExactComplex](value)
                return rebind_var[_Scalar[k]](Complex(_convert[4](z.real(), NoContext(), budget), _convert[4](z.imag(), NoContext(), budget)))
            else:
                return rebind_var[_Scalar[k]](Complex(_convert[4](value, NoContext(), budget), Real.zero(NoContext())))
    elif k == 3:
        comptime assert _Kind[T] == 0 or T == ModInt, "modular coefficients cannot mix implicitly with rational or complex coefficients"
        ref field = rebind[Modulus](ctx)
        comptime if T == ModInt:
            ref source = rebind[ModInt](value)
            if source.modulus != field.value:
                raise PolyError(PolyErrorKind.DOMAIN, "coefficients have different moduli")
            if source.value >= field.value:
                raise PolyError(PolyErrorKind.DOMAIN, "modular coefficient is not a canonical residue")
            return rebind_var[_Scalar[k]](source.copy())
        else:
            try:
                var reduced = _integer(value) % Integer(field.value)
                if reduced < 0:
                    reduced += Integer(field.value)
                return rebind_var[_Scalar[k]](ModInt(reduced.to_native_exact[DType.uint64](), field))
            except error:
                raise PolyError(PolyErrorKind.CAPACITY, String(error))
    elif k == 2:
        comptime assert T != ModInt, "modular coefficients have no implicit complex conversion"
        comptime if T == ExactComplex:
            return rebind_var[_Scalar[k]](rebind[ExactComplex](value).copy())
        else:
            return rebind_var[_Scalar[k]](ExactComplex(_rational(value)))
    elif k == 1:
        return rebind_var[_Scalar[k]](_rational(value))
    else:
        return rebind_var[_Scalar[k]](_integer(value))


def _wrap[k: Int](value: _Scalar[k]) -> _Ring[k]:
    comptime if k >= 4:
        return rebind_var[_Ring[k]](value.copy())
    elif k == 3:
        return rebind_var[_Ring[k]](rebind[ModInt](value).copy())
    elif k == 2:
        return rebind_var[_Ring[k]](QQi(rebind[ExactComplex](value)))
    elif k == 1:
        return rebind_var[_Ring[k]](QQ(rebind[Rational](value)))
    else:
        return rebind_var[_Ring[k]](ZZ(rebind[Integer](value)))


@always_inline
def _move_coefficient[k: Int](var value: _Scalar[k]) -> _Ring[k]:
    # The caller owns a canonical scalar and has already charged its transfer.
    comptime if k >= 3:
        return rebind_var[_Ring[k]](value^)
    elif k == 2:
        return rebind_var[_Ring[k]](QQi(rebind_var[ExactComplex](value^)))
    elif k == 1:
        return rebind_var[_Ring[k]](QQ(rebind_var[Rational](value^)))
    else:
        return rebind_var[_Ring[k]](ZZ(rebind_var[Integer](value^)))


@always_inline
def _take_coefficient[k: Int, T: Copyable & Deinitable](var value: T, ctx: _Ring[k].Context, mut budget: Budget) raises PolyError -> _Ring[k]:
    """Move an already canonical scalar; only promotion inspects its payload."""
    comptime if T == _Scalar[k] and k != 3:
        # APN scalars own shared payloads. Moving the handle is constant work,
        # independent of magnitude; there is no big-integer copy to precharge.
        budget.charge(1)
        return _move_coefficient[k](rebind_var[_Scalar[k]](value^))
    else:
        return _wrap[k](_convert[k](value, ctx, budget))


def _unwrap[k: Int](value: _Ring[k]) -> _Scalar[k]:
    comptime if k >= 4:
        return rebind_var[_Scalar[k]](value.copy())
    elif k == 3:
        return rebind_var[_Scalar[k]](rebind[ModInt](value).copy())
    elif k == 2:
        return rebind_var[_Scalar[k]](rebind[QQi](value).value.copy())
    elif k == 1:
        return rebind_var[_Scalar[k]](rebind[QQ](value).value.copy())
    else:
        return rebind_var[_Scalar[k]](rebind[ZZ](value).value.copy())


struct Values[T: Copyable & Deinitable, kind: Int](Copyable, Movable, Iterable, Sized, _WithModulus):
    """An immutable sequence from coefficients or evaluate_many, retaining every zero.

    Parameters:
        T: The sequence's exact scalar type.
        kind: The inferred exact coefficient kind retained by this sequence.
    """
    var _values: ArcPointer[List[Self.T]]
    var _modulus: UInt64
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _ValuesIterator[Self.T, Self.kind, ImmOrigin(iterable_origin)]
    """A borrowing iterator over the shared storage.

    Parameters:
        iterable_mut: Whether the borrowing origin is mutable.
        iterable_origin: The lifetime of the borrowed array.
    """

    def __init__(out self, var values: List[Self.T], modulus: UInt64, mut budget: Budget) raises PolyError:
        """Take ownership of values and their already validated field metadata.

        Args:
            values: The coefficient sequence, sparse terms, or prepared array.
            modulus: The prime modulus; supported primes are below 2**62.
            budget: The caller-owned work allowance, shared through preparation and computation.

        Raises:
            PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
        """
        _allocation[List[Self.T]](1, budget)
        self._values = ArcPointer(values^)
        self._modulus = modulus

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over values without copying the backing buffer.

        Returns:
            An iterator borrowing the shared storage.
        """
        return _ValuesIterator[Self.T, Self.kind, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def __len__(self) -> Int:
        """Return the number of stored values.

        Returns:
            The number of stored entries.
        """
        return len(self._values[])

    def __getitem__(self, index: Int) -> Self.T:
        """Return the entry at a valid zero-based index.

        Args:
            index: A valid zero-based index into stored entries.

        Returns:
            The scalar or sparse term at the requested index.
        """
        return self._values[][index].copy()

    def to_list(self, mut budget: Budget) raises PolyError -> List[Self.T]:
        """Return an independent mutable list, charging allocation and coefficient work.

        Args:
            budget: The caller-owned work allowance, shared through preparation and computation.

        Returns:
            An independent mutable list in the current entry order.

        Raises:
            PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
        """
        _allocation[Self.T](len(self), budget)
        var result = List[Self.T](capacity=len(self))
        for index in range(len(self)):
            _charge_size(_size(self._values[][index]), budget)
            result.append(self[index])
        return result^

    def to_list(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Self.T]:
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

    def modulus(self) -> UInt64:
        """Return the retained prime modulus, or zero outside a prime field.

        Returns:
            The prime modulus for modular data, or zero otherwise.
        """
        return self._modulus

    def domain(self) -> String:
        """Return ZZ, QQ, QQi, GF(p), AA, or QQbar, including for an empty result.

        Returns:
            ZZ, QQ, QQi, GF(p), AA, or QQbar, including for an empty result.
        """
        return _domain[Self.kind](self._modulus)


@fieldwise_init
struct _ValuesIterator[T: Copyable & Deinitable, kind: Int, origin: ImmOrigin](Iterator):
    comptime Element = Self.T
    var source: Pointer[Values[Self.T, Self.kind], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.T:
        if self.position == len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position]
        self.position += 1
        return result^

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        var remaining = len(self.source[]) - self.position
        return remaining, Optional(remaining)


def coefficients[*Ts: Copyable & Deinitable](values: Tuple[*Ts], mut budget: Budget) raises PolyError -> Values[_Scalar[_TupleKind[Tuple[*Ts]]], _TupleKind[Tuple[*Ts]]]:
    """Promote a heterogeneous tuple using the caller's work allowance.

    Parameters:
        Ts: The inferred element types of the input tuple.

    Args:
        values: The coefficient sequence, sparse terms, or prepared array.
        budget: The caller-owned work allowance, shared through preparation and computation.

    Returns:
        An immutable sequence in the common exact domain, with zeros retained.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    comptime T = Tuple[*Ts]
    comptime k = _TupleKind[T]
    var context = _Context[k]()
    var modulus = UInt64(0)
    comptime for i in range(T.Ts.length):
        budget.charge(1)
        comptime if T.Ts[i] == ModInt:
            modulus = _merge_modulus(modulus, rebind[ModInt](values[i]).modulus)
            context.observe(modulus, budget)
    var ctx = context.get()
    _allocation[_Scalar[k]](T.Ts.length, budget)
    var result = List[_Scalar[k]](capacity=T.Ts.length)
    comptime for i in range(T.Ts.length):
        result.append(_convert[k](values[i], ctx, budget))
    return Values[_Scalar[k], k](result^, modulus, budget)


def coefficients[*Ts: Copyable & Deinitable](values: Tuple[*Ts], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Values[_Scalar[_TupleKind[Tuple[*Ts]]], _TupleKind[Tuple[*Ts]]]:
    """Promote mixed tuple entries to one exact domain with an optional work limit.

    Parameters:
        Ts: The inferred element types of the input tuple.

    Args:
        values: The coefficient sequence, sparse terms, or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An immutable sequence in the common exact domain, with zeros retained.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    var budget = Budget(max_steps)
    return coefficients(values, budget)


def modint[T: Movable](value: T, modulus: UInt64, mut budget: Budget) raises PolyError -> ModInt:
    """Reduce an integer into a validated prime field.

    Parameters:
        T: The sequence's exact scalar type.

    Args:
        value: An exact integer to reduce modulo the prime.
        modulus: The prime modulus; supported primes are below 2**62.
        budget: The caller-owned work allowance, shared through preparation and computation.

    Returns:
        A canonical ModInt residue in the validated prime field.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    var ctx = Modulus(modulus, budget)
    return _convert[3](value, ctx, budget)


def modint[T: Movable](value: T, modulus: UInt64, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> ModInt:
    """Reduce an exact integer into a prime field with an optional work limit.

    Parameters:
        T: The sequence's exact scalar type.

    Args:
        value: An exact integer to reduce modulo the prime.
        modulus: The prime modulus; supported primes are below 2**62.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A canonical ModInt residue in the validated prime field.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    var budget = Budget(max_steps)
    return modint(value, modulus, budget)


def _sum(a: Int, b: Int) raises PolyError -> Int:
    if a < 0 or b < 0 or a > Int.MAX - b:
        raise PolyError(PolyErrorKind.CAPACITY, "size addition overflow")
    return a + b


def _product(a: Int, b: Int) raises PolyError -> Int:
    if a < 0 or b < 0 or (b != 0 and a > Int.MAX // b):
        raise PolyError(PolyErrorKind.CAPACITY, "size multiplication overflow")
    return a * b


def _charge_size(count: Int, mut budget: Budget) raises PolyError:
    if count < 0:
        raise PolyError(PolyErrorKind.CAPACITY, "negative size")
    budget.charge(UInt64(count))


def _allocation[T: AnyType](count: Int, mut budget: Budget) raises PolyError:
    _ = _product(count, size_of[T]())
    _charge_size(count, budget)


def _sequence[I: Iterable](values: I, mut budget: Budget) raises PolyError -> Values[_Scalar[_Kind[_Element[I]]], _Kind[_Element[I]]]:
    comptime k = _Kind[_Element[I]]
    comptime assert k >= 0, "sequence entries must be exact numbers; convert floats explicitly"
    _check_batch(values, message="scalar data must form a one-dimensional sequence")
    var context = _Context[k]()
    var modulus = _initial_modulus(values)
    if modulus:
        context.observe(modulus, budget)
    var cursor = iter(values)
    comptime T = downcast[type_of(cursor).Element, Copyable & Deinitable]
    var result = List[_Scalar[k]]()
    while True:
        budget.charge(1)
        var raw: type_of(cursor).Element
        try:
            raw = next(cursor)
        except StopIteration:
            break
        var value = rebind_var[T](raw^)
        comptime if T == ModInt:
            context.observe(rebind[ModInt](value).modulus, budget)
        var ctx = context.get()
        _grow(result, budget)
        result.append(_convert[k](value, ctx, budget))
    comptime if k == 3:
        if len(result) or modulus:
            modulus = rebind[Modulus](context.get()).value
    return Values[_Scalar[k], k](result^, modulus, budget)


def _sequence_coefficients[k: Int, T: Copyable & Deinitable, source: Int](values: Values[T, source], ctx: _Ring[k].Context, mut budget: Budget) raises PolyError -> List[_Ring[k]]:
    _allocation[_Ring[k]](len(values), budget)
    var result = List[_Ring[k]](capacity=len(values))
    for index in range(len(values)):
        budget.charge(1)
        result.append(_wrap[k](_convert[k](values._values[][index], ctx, budget)))
    return result^


def _values_from_core[k: Int](var values: List[_Ring[k]], ctx: _Ring[k].Context, mut budget: Budget) raises PolyError -> Values[_Scalar[k], k]:
    _allocation[_Scalar[k]](len(values), budget)
    var result = List[_Scalar[k]](capacity=len(values))
    for c in values:
        _charge_size(c.size(), budget)
        result.append(_unwrap[k](c))
    var modulus = UInt64(0)
    comptime if k == 3:
        modulus = rebind[Modulus](ctx).value
    return Values[_Scalar[k], k](result^, modulus, budget)


comptime _Nested[I: Iterable] = downcast[_Element[I], Iterable & Copyable & Deinitable]


def _data_context[k: Int](a: UInt64, b: UInt64, mut budget: Budget) raises PolyError -> _Ring[k].Context:
    var context = _Context[k]()
    var modulus = _merge_modulus(a, b)
    if modulus:
        context.observe(modulus, budget)
    return context.get()


def _rows[I: Iterable](values: I, mut budget: Budget) raises PolyError -> List[Values[_Scalar[_Kind[_Element[_Nested[I]]]], _Kind[_Element[_Nested[I]]]]]:
    comptime k = _Kind[_Element[_Nested[I]]]
    var result = List[Values[_Scalar[k], k]]()
    var cursor = iter(values)
    comptime T = downcast[type_of(cursor).Element, Iterable & Copyable & Deinitable]
    while True:
        budget.charge(1)
        var raw: type_of(cursor).Element
        try:
            raw = next(cursor)
        except StopIteration:
            break
        var row = rebind_var[_Nested[I]](raw^)
        _grow(result, budget)
        result.append(_sequence(row, budget))
    return result^
