"""Functional maps, broadcasting and folds over shaped collections of values."""
from ._io_meta import _Record, _Batch
from std.iter import Iterable, Iterator, StopIteration, iter, next
from std.memory import ArcPointer
from std.builtin.rebind import downcast
from polynomials import Budget, PolyError, PolyErrorKind
from ._common import _Element, _allocation, _grow, _product, _charge_size, _compact, _sum


@fieldwise_init
struct _Layout(Copyable, Movable):
    var shape: List[Int]
    var strides: List[Int]
    var count: Int


def _layout(shape: List[Int], mut budget: Budget) raises PolyError -> _Layout:
    _allocation[Int](len(shape), budget)
    _allocation[Int](len(shape), budget)
    var strides = List[Int](length=len(shape), fill=0)
    var empty = False
    for dimension in shape:
        budget.charge(1)
        if dimension < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "batch dimensions must be nonnegative")
        empty = empty or dimension == 0
    var count = 0 if empty else 1
    if not empty:
        for i in range(len(shape) - 1, -1, -1):
            budget.charge(1)
            strides[i] = count
            count = _product(count, shape[i])
    return _Layout(shape.copy(), strides^, count)


@fieldwise_init
struct Array[T: Copyable & Deinitable](Copyable, Movable, Iterable, Sized, _Record):
    """An immutable row-major batch; each item can be an entire polynomial.

    Parameters:
        T: The scalar or batch item type.
    """
    comptime _Format = _Batch
    comptime _Item = Self.T
    comptime _Power = Int
    var _items: ArcPointer[List[Self.T]]
    var _layout: ArcPointer[_Layout]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _Iterator[Self.T, ImmOrigin(iterable_origin)]
    """A row-major iterator over the logical batch, including broadcast views.

    Parameters:
        iterable_mut: Whether the borrowing source origin allows mutation.
        iterable_origin: The lifetime of the borrowed source.

    Returns:
        The borrowing iterator type.
    """

    def __len__(self) -> Int:
        """Return the logical number of batch items.

        Returns:
            The number of logical entries.
        """
        return self._layout[].count

    def ndim(self) -> Int:
        """Return the number of batch dimensions.

        Returns:
            The number of batch axes.
        """
        return len(self._layout[].shape)

    def shape(self) -> ref[ImmOrigin(origin_of(self._layout[].shape))] List[Int]:
        """Borrow the immutable batch shape.

        Returns:
            An immutable reference to the logical dimensions.
        """
        return self._layout[].shape

    def _index(self, var index: Int) -> Int:
        var result = 0
        for i in range(self.ndim() - 1, -1, -1):
            result += (index % self._layout[].shape[i]) * self._layout[].strides[i]
            index //= self._layout[].shape[i]
        return result

    def __getitem__(self, index: Int) -> Self.T:
        """Copy the item at a valid flat index in row-major order.

        Args:
            index: A valid zero-based index into the logical entries.

        Returns:
            An independent value copy; polynomial and evidence buffers may remain shared.
        """
        return self._items[][self._index(index)].copy()

    def item(self, mut budget: Budget) raises PolyError -> Self.T:
        """Copy the only item of a batch containing exactly one value.

        Args:
            budget: The shared allowance for locating and copying the item.

        Returns:
            The value, independent of the batch rank.

        Raises:
            PolyError: DOMAIN unless the batch has one item, or LIMIT on exhausted work.
        """
        _charge_size(_sum(self.ndim(), 1), budget)
        if len(self) != 1:
            raise PolyError(PolyErrorKind.DOMAIN, "item requires exactly one batch value")
        return self[0]

    def item(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Self.T:
        """Copy the only item, including the result of an all-axis reduction.

        Args:
            max_steps: The copying allowance; defaults to UInt64.MAX.

        Returns:
            The value, independent of the batch rank.

        Raises:
            PolyError: DOMAIN unless the batch has one item, or LIMIT on exhausted work.
        """
        var budget = Budget(max_steps)
        return self.item(budget)

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over logical batch items in row-major order.

        Returns:
            An iterator borrowing the immutable source storage.
        """
        return _Iterator[Self.T, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)


@fieldwise_init
struct _Iterator[T: Copyable & Deinitable, origin: ImmOrigin](Iterator):
    comptime Element = Self.T
    var source: Pointer[Array[Self.T], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.T:
        if self.position == len(self.source[]):
            raise StopIteration()
        var value = self.source[][self.position]
        self.position += 1
        return value^


def _take[T: Copyable & Deinitable](var items: List[T], var layout: _Layout, mut budget: Budget) raises PolyError -> Array[T]:
    _compact(items, budget)
    _allocation[List[T]](1, budget)
    _allocation[_Layout](1, budget)
    return Array[T](ArcPointer(items^), ArcPointer(layout^))


def asarray[I: Iterable](values: I, mut budget: Budget, *, shape: Optional[List[Int]] = None) raises PolyError -> Array[downcast[_Element[I], Copyable & Deinitable]]:
    """Collect one pass of batch items; shape defaults to a one-dimensional batch.

    Parameters:
        I: The inferred input iterable type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        budget: The shared work allowance, including preparation, computation, and result storage.
        shape: The row-major batch dimensions; the product must match the logical item count.

    Returns:
        An immutable prepared array in the selected domain or batch shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    comptime T = downcast[_Element[I], Copyable & Deinitable]
    comptime if I == Array[T]:
        if not shape:
            budget.charge(1)
            return rebind[Array[T]](values).copy()
    var items = List[T]()
    var cursor = iter(values)
    while True:
        budget.charge(1)
        var value: T
        try:
            value = rebind_var[T](next(cursor))
        except StopIteration:
            break
        _grow(items, budget)
        items.append(rebind_var[T](value^))
    var layout: _Layout
    if shape:
        layout = _layout(shape.value(), budget)
    else:
        _allocation[Int](1, budget)
        layout = _layout([len(items)], budget)
    if layout.count != len(items):
        raise PolyError(PolyErrorKind.DOMAIN, "batch shape does not match the item count")
    return _take(items^, layout^, budget)


def scalar[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> Array[T]:
    """Make a zero-dimensional batch that broadcasts one value.

    Parameters:
        T: The scalar or batch item type.

    Args:
        value: The value to convert or inspect.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        A zero-dimensional batch sharing one logical item.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    _allocation[T](1, budget)
    return _take(List[T](length=1, fill=value), _layout(List[Int](), budget), budget)


def broadcast_to[T: Copyable & Deinitable](values: Array[T], shape: List[Int], mut budget: Budget) raises PolyError -> Array[T]:
    """Create a shared view using NumPy's trailing-axis broadcasting rules.

    Parameters:
        T: The scalar or batch item type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        shape: The row-major batch dimensions; the product must match the logical item count.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        A shared broadcast view with the requested shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    if len(shape) < values.ndim():
        raise PolyError(PolyErrorKind.DOMAIN, "broadcasting cannot remove batch axes")
    var layout = _layout(shape, budget)
    for i in range(len(shape)):
        budget.charge(1)
        var source = i - (len(shape) - values.ndim())
        if source < 0:
            layout.strides[i] = 0
        elif values.shape()[source] == shape[i]:
            layout.strides[i] = values._layout[].strides[source]
        elif values.shape()[source] == 1:
            layout.strides[i] = 0
        else:
            raise PolyError(PolyErrorKind.DOMAIN, "incompatible batch shapes")
    _allocation[_Layout](1, budget)
    return Array[T](values._items.copy(), ArcPointer(layout^))


def _broadcast_shape(a: List[Int], b: List[Int], mut budget: Budget) raises PolyError -> List[Int]:
    var rank = max(len(a), len(b))
    _allocation[Int](rank, budget)
    var shape = List[Int](length=rank, fill=1)
    for i in range(rank):
        budget.charge(1)
        var ai = i - (rank - len(a))
        var bi = i - (rank - len(b))
        var x = a[ai] if ai >= 0 else 1
        var y = b[bi] if bi >= 0 else 1
        if x != y and x != 1 and y != 1:
            raise PolyError(PolyErrorKind.DOMAIN, "incompatible batch shapes")
        shape[i] = y if x == 1 else x
    return shape^


def vmap[T: Copyable & Deinitable, R: Copyable & Deinitable](function: def(x: T, mut budget: Budget) raises PolyError thin -> R, values: Array[T], mut budget: Budget) raises PolyError -> Array[R]:
    """Apply a budget-aware function to every batch item, preserving shape.

    Parameters:
        T: The scalar or batch item type.
        R: The callback result type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Callback results in row-major order and the resulting broadcast shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var layout = _layout(values.shape(), budget)
    _allocation[R](len(values), budget)
    var output = List[R](capacity=len(values))
    for i in range(len(values)):
        _charge_size(values.ndim() + 1, budget)
        output.append(function(values._items[][values._index(i)], budget))
    return _take(output^, layout^, budget)


def vmap[A: Copyable & Deinitable, B: Copyable & Deinitable, R: Copyable & Deinitable](function: def(a: A, b: B, mut budget: Budget) raises PolyError thin -> R, a: Array[A], b: Array[B], mut budget: Budget) raises PolyError -> Array[R]:
    """Apply a binary function after broadcasting both batch shapes.

    Parameters:
        A: The first input type.
        B: The second input type.
        R: The callback result type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        a: The first input value.
        b: The second input value.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Callback results in row-major order and the resulting broadcast shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var shape = _broadcast_shape(a.shape(), b.shape(), budget)
    var left = broadcast_to(a, shape, budget)
    var right = broadcast_to(b, shape, budget)
    var layout = _layout(shape, budget)
    _allocation[R](len(left), budget)
    var output = List[R](capacity=len(left))
    for i in range(len(left)):
        _charge_size(2 * len(shape) + 1, budget)
        output.append(function(left._items[][left._index(i)], right._items[][right._index(i)], budget))
    return _take(output^, layout^, budget)


def outer[A: Copyable & Deinitable, B: Copyable & Deinitable, R: Copyable & Deinitable](function: def(a: A, b: B, mut budget: Budget) raises PolyError thin -> R, a: Array[A], b: Array[B], mut budget: Budget) raises PolyError -> Array[R]:
    """Apply a binary function to every pair; the output shape concatenates both shapes.

    Parameters:
        A: The first input type.
        B: The second input type.
        R: The callback result type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        a: The first input value.
        b: The second input value.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Every pairwise result with the concatenated input shapes.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    _allocation[Int](a.ndim() + b.ndim(), budget)
    var shape = List[Int](capacity=a.ndim() + b.ndim())
    for dimension in a.shape():
        budget.charge(1)
        shape.append(dimension)
    for dimension in b.shape():
        budget.charge(1)
        shape.append(dimension)
    var layout = _layout(shape, budget)
    _allocation[R](layout.count, budget)
    var output = List[R](capacity=layout.count)
    for i in range(len(a)):
        for j in range(len(b)):
            _charge_size(a.ndim() + b.ndim() + 1, budget)
            output.append(function(a._items[][a._index(i)], b._items[][b._index(j)], budget))
    return _take(output^, layout^, budget)


def _mark_axis(var axis: Int, mut selected: List[Bool], mut budget: Budget) raises PolyError:
    budget.charge(1)
    var rank = len(selected)
    if axis < -rank or axis >= rank:
        raise PolyError(PolyErrorKind.DOMAIN, "reduction axis is outside the batch rank")
    if axis < 0:
        axis += rank
    if selected[axis]:
        raise PolyError(PolyErrorKind.DOMAIN, "reduction axes must be distinct")
    selected[axis] = True


trait _TupleAxes:
    def _mark_reduce_axes(self, mut selected: List[Bool], mut budget: Budget) raises PolyError: ...


__extension Tuple(_TupleAxes):
    def _mark_reduce_axes(self, mut selected: List[Bool], mut budget: Budget) raises PolyError:
        comptime for i in range(Self.Ts.length):
            comptime assert Self.Ts[i] == Int, "reduction axes must be Int"
            _mark_axis(rebind[Int](self[i]), selected, budget)


def _axes[A: AnyType](axis: A, rank: Int, mut budget: Budget) raises PolyError -> List[Bool]:
    _allocation[Bool](rank, budget)
    var selected = List[Bool](length=rank, fill=False)
    comptime if A == NoneType or A == type_of(None):
        for i in range(rank):
            budget.charge(1)
            selected[i] = True
    elif A == Int:
        _mark_axis(rebind[Int](axis), selected, budget)
    elif A == Optional[Int]:
        ref choice = rebind[Optional[Int]](axis)
        if choice:
            _mark_axis(choice.value(), selected, budget)
        else:
            for i in range(rank):
                budget.charge(1)
                selected[i] = True
    elif conforms_to(A, _TupleAxes):
        rebind[downcast[A, _TupleAxes]](axis)._mark_reduce_axes(selected, budget)
    elif conforms_to(A, Iterable):
        comptime I = downcast[A, Iterable]
        comptime assert _Element[I] == Int, "reduction axes must be Int"
        var cursor = iter(rebind[I](axis))
        while True:
            budget.charge(1)
            var value: Int
            try:
                value = rebind_var[Int](next(cursor))
            except StopIteration:
                break
            _mark_axis(value, selected, budget)
    else:
        comptime assert False, "axis must be None, Int, Optional[Int], or Int axes"
    return selected^


def _sub_layout(layout: _Layout, selected: List[Bool], reduced: Bool, mut budget: Budget) raises PolyError -> _Layout:
    var rank = len(layout.shape)
    _allocation[Int](rank, budget)
    _allocation[Int](rank, budget)
    var shape = List[Int](capacity=rank)
    var strides = List[Int](capacity=rank)
    var empty = False
    for i in range(rank):
        budget.charge(1)
        if selected[i] == reduced:
            shape.append(layout.shape[i])
            strides.append(layout.strides[i])
            empty = empty or layout.shape[i] == 0
    var count = 0 if empty else 1
    if not empty:
        for dimension in shape:
            budget.charge(1)
            count = _product(count, dimension)
    return _Layout(shape^, strides^, count)


def _offset(layout: _Layout, var index: Int) -> Int:
    # Called only for nonempty layouts. These strides belong to disjoint axes
    # of a validated source, so their summed offsets stay inside its storage.
    var position = 0
    for i in range(len(layout.shape) - 1, -1, -1):
        position += (index % layout.shape[i]) * layout.strides[i]
        index //= layout.shape[i]
    return position


def _reduce_axes[T: Copyable & Deinitable, A: AnyType](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, axis: A, keepdims: Bool, mut budget: Budget) raises PolyError -> Array[T]:
    budget.charge(1)
    var selected = _axes(axis, values.ndim(), budget)
    var kept = _sub_layout(values._layout[], selected, False, budget)
    _allocation[Int](values.ndim(), budget)
    var shape = List[Int](capacity=values.ndim())
    for i in range(values.ndim()):
        budget.charge(1)
        if not selected[i]:
            shape.append(values.shape()[i])
        elif keepdims:
            shape.append(1)
    var layout = _layout(shape, budget)
    _allocation[T](layout.count, budget)
    var output = List[T](capacity=layout.count)
    if kept.count == 0:
        return _take(output^, layout^, budget)
    var reduced = _sub_layout(values._layout[], selected, True, budget)
    for i in range(kept.count):
        _charge_size(_sum(len(kept.shape), 1), budget)
        var base = _offset(kept, i)
        var current = initial.copy()
        for j in range(reduced.count):
            _charge_size(_sum(len(reduced.shape), 1), budget)
            current = function(current, values._items[][base + _offset(reduced, j)], budget)
        output.append(current^)
    return _take(output^, layout^, budget)


def reduce[T: Copyable & Deinitable](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, mut budget: Budget) raises PolyError -> T:
    """Fold row-major batch items from an explicit initial value, including empty batches.

    Parameters:
        T: The scalar or batch item type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        initial: The initial accumulator, also returned for an empty reduction.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        The final left-fold accumulator.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    budget.charge(1)
    var result = initial.copy()
    for i in range(len(values)):
        _charge_size(values.ndim() + 1, budget)
        result = function(result, values._items[][values._index(i)], budget)
    return result^


def accumulate[T: Copyable & Deinitable](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, mut budget: Budget) raises PolyError -> Array[T]:
    """Return each running row-major fold result in the original batch shape.

    Parameters:
        T: The scalar or batch item type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        initial: The initial accumulator, also returned for an empty reduction.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Each running fold result in the original batch shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var layout = _layout(values.shape(), budget)
    _allocation[T](len(values), budget)
    var output = List[T](capacity=len(values))
    var current = initial.copy()
    for i in range(len(values)):
        _charge_size(values.ndim() + 1, budget)
        current = function(current, values._items[][values._index(i)], budget)
        output.append(current.copy())
    return _take(output^, layout^, budget)


def asarray[I: Iterable](values: I, *, shape: Optional[List[Int]] = None, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[downcast[_Element[I], Copyable & Deinitable]]:
    """Collect one pass of batch items; shape defaults to a one-dimensional batch.

    Parameters:
        I: The inferred input iterable type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        shape: The row-major batch dimensions; the product must match the logical item count.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An immutable prepared array in the selected domain or batch shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return asarray(values, budget, shape=shape)


def scalar[T: Copyable & Deinitable](value: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[T]:
    """Make a zero-dimensional batch that broadcasts one value.

    Parameters:
        T: The scalar or batch item type.

    Args:
        value: The value to convert or inspect.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A zero-dimensional batch sharing one logical item.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return scalar(value, budget)


def broadcast_to[T: Copyable & Deinitable](values: Array[T], shape: List[Int], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[T]:
    """Create a shared view using NumPy's trailing-axis broadcasting rules.

    Parameters:
        T: The scalar or batch item type.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        shape: The row-major batch dimensions; the product must match the logical item count.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A shared broadcast view with the requested shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return broadcast_to(values, shape, budget)


def vmap[T: Copyable & Deinitable, R: Copyable & Deinitable](function: def(x: T, mut budget: Budget) raises PolyError thin -> R, values: Array[T], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[R]:
    """Apply a budget-aware function to every batch item, preserving shape.

    Parameters:
        T: The scalar or batch item type.
        R: The callback result type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Callback results in row-major order and the resulting broadcast shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return vmap[T, R](function, values, budget)


def vmap[A: Copyable & Deinitable, B: Copyable & Deinitable, R: Copyable & Deinitable](function: def(a: A, b: B, mut budget: Budget) raises PolyError thin -> R, a: Array[A], b: Array[B], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[R]:
    """Apply a binary function after broadcasting both batch shapes.

    Parameters:
        A: The first input type.
        B: The second input type.
        R: The callback result type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        a: The first input value.
        b: The second input value.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Callback results in row-major order and the resulting broadcast shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return vmap[A, B, R](function, a, b, budget)


def outer[A: Copyable & Deinitable, B: Copyable & Deinitable, R: Copyable & Deinitable](function: def(a: A, b: B, mut budget: Budget) raises PolyError thin -> R, a: Array[A], b: Array[B], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[R]:
    """Apply a binary function to every pair; the output shape concatenates both shapes.

    Parameters:
        A: The first input type.
        B: The second input type.
        R: The callback result type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        a: The first input value.
        b: The second input value.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Every pairwise result with the concatenated input shapes.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return outer[A, B, R](function, a, b, budget)


def reduce[T: Copyable & Deinitable](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> T:
    """Fold row-major batch items from an explicit initial value, including empty batches.

    Parameters:
        T: The scalar or batch item type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        initial: The initial accumulator, also returned for an empty reduction.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The final left-fold accumulator.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return reduce(function, values, initial, budget)


def accumulate[T: Copyable & Deinitable](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[T]:
    """Return each running row-major fold result in the original batch shape.

    Parameters:
        T: The scalar or batch item type.

    Args:
        function: A thin Mojo callback with a final shared Budget argument and typed PolyError failures.
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        initial: The initial accumulator, also returned for an empty reduction.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Each running fold result in the original batch shape.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return accumulate(function, values, initial, budget)


def reduce[T: Copyable & Deinitable, A: AnyType](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, mut budget: Budget, *, axis: A, keepdims: Bool = False) raises PolyError -> Array[T]:
    """Reduce selected batch axes with one ordered fold per output item.

    Parameters:
        T: The batch item and accumulator type.
        A: The inferred axis argument type.

    Args:
        function: A thin callback taking two items and the shared Budget, raising PolyError.
        values: The input batch; polynomial variables are not batch axes.
        initial: A fresh starting value for each fold, returned for an empty group.
        budget: The shared allowance for axis validation, storage, and callbacks.
        axis: None for all axes, an Int, or a tuple/iterable of Int axes. Negative axes count from the end; Optional[Int] is also accepted. Duplicate or out-of-range axes raise DOMAIN.
        keepdims: Keep reduced axes as dimensions of size one.

    Returns:
        A batch, including a zero-dimensional batch when all axes are removed.
        Outputs and each fold visit their axes in row-major order, independent
        of axis-list order. An empty axis list applies function(initial, item)
        to every input. Use item() to extract a result containing one value.

    Raises:
        PolyError: DOMAIN for invalid axes or callback failures, LIMIT on
            exhausted work, or CAPACITY for unrepresentable output storage.
    """
    return _reduce_axes(function, values, initial, axis, keepdims, budget)


def reduce[T: Copyable & Deinitable](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, mut budget: Budget, *, keepdims: Bool) raises PolyError -> Array[T]:
    """Reduce all batch axes and return a batch with the requested dimensions.

    Parameters:
        T: The batch item and accumulator type.

    Args:
        function: A thin callback taking two items and the shared Budget, raising PolyError.
        values: The input batch; polynomial variables are not batch axes.
        initial: A fresh starting value for each fold, returned for an empty group.
        budget: The shared allowance for axis validation, storage, and callbacks.
        keepdims: Keep reduced axes as dimensions of size one.

    Returns:
        A batch, including a zero-dimensional batch when all axes are removed.
        Outputs and each fold visit their axes in row-major order, independent
        of axis-list order. An empty axis list applies function(initial, item)
        to every input. Use item() to extract a result containing one value.

    Raises:
        PolyError: DOMAIN for invalid axes or callback failures, LIMIT on
            exhausted work, or CAPACITY for unrepresentable output storage.
    """
    return _reduce_axes(function, values, initial, None, keepdims, budget)


def reduce[T: Copyable & Deinitable, A: AnyType](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, *, axis: A, keepdims: Bool = False, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[T]:
    """Reduce selected batch axes with one ordered fold per output item.

    Parameters:
        T: The batch item and accumulator type.
        A: The inferred axis argument type.

    Args:
        function: A thin callback taking two items and the shared Budget, raising PolyError.
        values: The input batch; polynomial variables are not batch axes.
        initial: A fresh starting value for each fold, returned for an empty group.
        axis: None for all axes, an Int, or a tuple/iterable of Int axes. Negative axes count from the end; Optional[Int] is also accepted. Duplicate or out-of-range axes raise DOMAIN.
        keepdims: Keep reduced axes as dimensions of size one.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A batch, including a zero-dimensional batch when all axes are removed.
        Outputs and each fold visit their axes in row-major order, independent
        of axis-list order. An empty axis list applies function(initial, item)
        to every input. Use item() to extract a result containing one value.

    Raises:
        PolyError: DOMAIN for invalid axes or callback failures, LIMIT on
            exhausted work, or CAPACITY for unrepresentable output storage.
    """
    var budget = Budget(max_steps)
    return _reduce_axes(function, values, initial, axis, keepdims, budget)


def reduce[T: Copyable & Deinitable](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, *, keepdims: Bool, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[T]:
    """Reduce all batch axes and return a batch with the requested dimensions.

    Parameters:
        T: The batch item and accumulator type.

    Args:
        function: A thin callback taking two items and the shared Budget, raising PolyError.
        values: The input batch; polynomial variables are not batch axes.
        initial: A fresh starting value for each fold, returned for an empty group.
        keepdims: Keep reduced axes as dimensions of size one.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A batch, including a zero-dimensional batch when all axes are removed.
        Outputs and each fold visit their axes in row-major order, independent
        of axis-list order. An empty axis list applies function(initial, item)
        to every input. Use item() to extract a result containing one value.

    Raises:
        PolyError: DOMAIN for invalid axes or callback failures, LIMIT on
            exhausted work, or CAPACITY for unrepresentable output storage.
    """
    var budget = Budget(max_steps)
    return _reduce_axes(function, values, initial, None, keepdims, budget)


def _accumulate_axis[T: Copyable & Deinitable, A: AnyType](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, axis: A, mut budget: Budget) raises PolyError -> Array[T]:
    budget.charge(1)
    var choice = Optional[Int]()
    comptime if A == Int:
        choice = rebind[Int](axis)
    elif A == Optional[Int]:
        choice = rebind[Optional[Int]](axis).copy()
    else:
        comptime assert A == NoneType or A == type_of(None), "accumulation axis must be None, Int, or Optional[Int]"
    if not choice:
        return accumulate(function, values, initial, budget)
    var index = choice.value()
    var rank = values.ndim()
    if index < -rank or index >= rank:
        raise PolyError(PolyErrorKind.DOMAIN, "accumulation axis is outside the batch rank")
    if index < 0:
        index += rank
    var layout = _layout(values.shape(), budget)
    _allocation[T](len(values), budget)
    var output = List[T](capacity=len(values))
    for i in range(len(values)):
        _charge_size(_sum(rank, 1), budget)
        var stride = layout.strides[index]
        var current: T
        if (i // stride) % layout.shape[index] == 0:
            current = initial.copy()
        else:
            current = output[i - stride].copy()
        current = function(current, values._items[][values._index(i)], budget)
        output.append(current^)
    return _take(output^, layout^, budget)


def accumulate[T: Copyable & Deinitable, A: AnyType](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, mut budget: Budget, *, axis: A) raises PolyError -> Array[T]:
    """Return running folds along one batch axis, preserving the input shape.

    Parameters:
        T: The batch item and accumulator type.
        A: The inferred axis argument type.

    Args:
        function: A thin callback taking two items and the shared Budget.
        values: The input batch, including broadcast views.
        initial: A fresh starting value for each line along the selected axis.
        budget: The shared allowance for validation, storage, and callbacks.
        axis: An Int or Optional[Int]; negative axes count from the end. None scans the flattened batch and preserves its shape.

    Returns:
        Every running fold, in the original shape. Callbacks visit items in
        row-major order; each line starts at initial. Empty batches remain empty.

    Raises:
        PolyError: DOMAIN for invalid axes or callback failures, LIMIT on
            exhausted work, or CAPACITY for unrepresentable output storage.
    """
    return _accumulate_axis(function, values, initial, axis, budget)


def accumulate[T: Copyable & Deinitable, A: AnyType](function: def(a: T, b: T, mut budget: Budget) raises PolyError thin -> T, values: Array[T], initial: T, *, axis: A, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[T]:
    """Return running folds along one batch axis, preserving the input shape.

    Parameters:
        T: The batch item and accumulator type.
        A: The inferred axis argument type.

    Args:
        function: A thin callback taking two items and the shared Budget.
        values: The input batch, including broadcast views.
        initial: A fresh starting value for each line along the selected axis.
        axis: An Int or Optional[Int]; negative axes count from the end. None scans the flattened batch and preserves its shape.
        max_steps: The complete call's work allowance.

    Returns:
        Every running fold, in the original shape. Callbacks visit items in
        row-major order; each line starts at initial. Empty batches remain empty.

    Raises:
        PolyError: DOMAIN for invalid axes or callback failures, LIMIT on
            exhausted work, or CAPACITY for unrepresentable output storage.
    """
    var budget = Budget(max_steps)
    return _accumulate_axis(function, values, initial, axis, budget)
