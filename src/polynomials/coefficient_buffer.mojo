"""Owned coefficient storage and allocation-free views of dense inputs."""
from std.iter import Iterable, Iterator, StopIteration
from .errors import Budget, PolyError, _allocation, _charge_size, _product


comptime _FIRST = 4
comptime _SECOND = 1024
comptime _PREFIX = _FIRST + _SECOND
comptime _BLOCK = 32768


struct CoefficientBuffer[T: Copyable & Deinitable](Copyable, Movable, Sized, Iterable, Equatable where conforms_to(T, Equatable), Writable where conforms_to(T, Writable)):
    """Ascending coefficients in contiguous or chunked owned storage.

    Copies own independent storage. Polynomial arrays share their enclosing
    immutable polynomial handle instead. Indexing never flattens chunks.

    Parameters:
        T: The coefficient type.
    """
    comptime _block_size = _BLOCK
    comptime _initial_size = _FIRST
    comptime _second_size = _SECOND
    var _flat: List[Self.T]
    var _chunks: List[List[Self.T]]
    var _count: Int
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _BufferIterator[Self.T, ImmOrigin(iterable_origin)]
    """The iterator that borrows the owner and yields independent values.

    Parameters:
        iterable_mut: Whether the source origin is mutable.
        iterable_origin: The source lifetime.
    """

    def __init__(out self):
        """Construct an empty contiguous buffer."""
        self._flat = List[Self.T]()
        self._chunks = List[List[Self.T]]()
        self._count = 0

    @implicit
    def __init__(out self, var values: List[Self.T]):
        """Adopt a list without copying its coefficient storage.

        Args:
            values: The owned ascending coefficients.
        """
        self._count = len(values)
        self._flat = values^
        self._chunks = List[List[Self.T]]()

    def __init__(out self, *, copy: Self):
        """Copy coefficients into independently mutable storage.

        Args:
            copy: The source buffer.
        """
        self._flat = List[Self.T](capacity=len(copy))
        if len(copy._chunks):
            for chunk in copy._chunks:
                self._flat.extend(Span(chunk))
        else:
            self._flat.extend(Span(copy._flat))
        self._chunks = List[List[Self.T]]()
        self._count = len(self._flat)

    @staticmethod
    def _from_chunks(var chunks: List[List[Self.T]], count: Int) -> Self:
        # Only the checked collector calls this: all but the last chunk have
        # their scheduled capacity (_FIRST, _SECOND, then _BLOCK). The last is
        # nonempty, and count is their checked sum.
        var result = Self()
        if len(chunks) == 1:
            result._flat = chunks.pop()
        else:
            result._chunks = chunks^
        result._count = count
        return result^

    def __len__(self) -> Int:
        """Return the initialized coefficient count.

        Returns:
            The number of coefficients.
        """
        return self._count

    @__unsafe_nested_origins_read_only
    def __getitem__(ref self, index: Int) -> ref[origin_of(self._flat)._get_owned_interior["element"], origin_of(self._chunks)._get_owned_interior["element"]._get_owned_interior["element"]] Self.T:
        """Borrow a coefficient, without copying or flattening.

        Args:
            index: A nonnegative index smaller than the buffer length.

        Returns:
            The coefficient reference, with the source mutability and lifetime.
        """
        debug_assert(index >= 0 and index < self._count, "coefficient index out of bounds")
        if len(self._chunks):
            if index < _FIRST:
                return self._chunks[0][index]
            if index < _PREFIX:
                return self._chunks[1][index - _FIRST]
            return self._chunks[2 + (index - _PREFIX) // _BLOCK][(index - _PREFIX) % _BLOCK]
        return self._flat[index]

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Borrow the buffer and yield independent coefficient values.

        Returns:
            An iterator over ascending values; coefficient copies use their usual value semantics.
        """
        return _BufferIterator[Self.T, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def __eq__(self, other: Self) -> Bool where conforms_to(Self.T, Equatable):
        """Compare coefficient values independently of storage layout.

        Args:
            other: The other buffer.

        Returns:
            Whether lengths and every coefficient agree.
        """
        if len(self) != len(other):
            return False
        for i in range(len(self)):
            if self[i] != other[i]:
                return False
        return True

    def write_to[W: Writer](self, mut writer: W) where conforms_to(Self.T, Writable):
        """Write the coefficients in list notation.

        Parameters:
            W: The output writer type.

        Args:
            writer: The destination receiving the representation.
        """
        writer.write("[")
        for i in range(len(self)):
            if i:
                writer.write(", ")
            writer.write(self[i])
        writer.write("]")

    def capacity(self) -> Int:
        """Return the total allocated coefficient capacity.

        Returns:
            Initialized and spare slots, excluding chunk metadata.
        """
        if len(self._chunks):
            var complete = _FIRST if len(self._chunks) == 2 else _PREFIX + (len(self._chunks) - 3) * _BLOCK
            return complete + self._chunks[len(self._chunks) - 1].capacity()
        return self._flat.capacity()

    def is_contiguous(self) -> Bool:
        """Test whether all coefficients occupy one list.

        Returns:
            True for contiguous storage, including an empty buffer.
        """
        return len(self._chunks) == 0

    def clear(mut self):
        """Release all coefficients and chunk metadata, leaving an empty buffer."""
        self = Self()

    def pop(mut self) -> Self.T:
        """Remove and return the final coefficient.

        Returns:
            The owned coefficient. The buffer must be nonempty.
        """
        debug_assert(self._count > 0, "cannot pop an empty coefficient buffer")
        self._count -= 1
        if len(self._chunks):
            var value = self._chunks[len(self._chunks) - 1].pop()
            if len(self._chunks[len(self._chunks) - 1]) == 0:
                _ = self._chunks.pop()
            if len(self._chunks) == 1:
                self._flat = self._chunks.pop()
                self._chunks = List[List[Self.T]]()
            return value^
        return self._flat.pop()

    def append(mut self, var value: Self.T):
        """Append a coefficient, preserving the current storage form.

        Args:
            value: The owned coefficient. Callers must charge growth first.
        """
        if len(self._chunks):
            var limit = _SECOND if len(self._chunks) == 2 else _BLOCK
            if len(self._chunks[len(self._chunks) - 1]) == limit:
                self._chunks.append(List[Self.T](capacity=_BLOCK))
            else:
                var capacity = self._chunks[len(self._chunks) - 1].capacity()
                if len(self._chunks[len(self._chunks) - 1]) == capacity:
                    # A compacted tail may have a non-power-of-two capacity.
                    # Cap regrowth so completed chunks keep their exact extent.
                    self._chunks[len(self._chunks) - 1].reserve(min(limit, max(1, capacity * 2)))
            self._chunks[len(self._chunks) - 1].append(value^)
        else:
            self._flat.append(value^)
        self._count += 1

    def to_list(self, mut budget: Budget) raises PolyError -> List[Self.T]:
        """Copy into a mutable contiguous list with explicit allocation work.

        Args:
            budget: The shared work allowance for allocation and handle copies.

        Returns:
            An independent list in coefficient order.

        Raises:
            PolyError: LIMIT or CAPACITY before unbudgeted work or allocation.
        """
        _allocation[Self.T](len(self), budget)
        _charge_size(len(self), budget)
        var result = List[Self.T](capacity=len(self))
        if len(self._chunks):
            for chunk in self._chunks:
                result.extend(Span(chunk))
        else:
            result.extend(Span(self._flat))
        return result^

    def _compact(mut self, mut budget: Budget) raises PolyError:
        # Pop has already released empty chunks. Only the last allocation can
        # retain spare coefficient capacity after canonical zero trimming.
        if self._count == 0:
            self = Self()
            return
        if len(self._chunks):
            _compact_list(self._chunks[len(self._chunks) - 1], budget)
            _compact_list(self._chunks, budget)
        else:
            _compact_list(self._flat, budget)


def _compact_list[T: Copyable & Deinitable](mut values: List[T], mut budget: Budget) raises PolyError:
    var count = len(values)
    if count == 0:
        values = List[T]()
    elif values.capacity() - count > count:
        _allocation[T](count, budget)
        _charge_size(_product(count, 2), budget)
        var compacted = List[T](capacity=count)
        compacted.extend(values^)
        values = compacted^


@fieldwise_init
struct _BufferIterator[T: Copyable & Deinitable, origin: ImmOrigin](Iterator):
    # Keep a pointer to the live owner, never a cached element address or extent.
    # This supports std.iter.iter on Mojo 1.1 without erasing interior lifetimes.
    # An emitted coefficient owns its value even if the source later changes.
    comptime Element = Self.T
    var source: Pointer[CoefficientBuffer[Self.T], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.T:
        if self.position >= len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position].copy()
        self.position += 1
        return result^

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        var remaining = max(0, len(self.source[]) - self.position)
        return remaining, Optional(remaining)


@fieldwise_init
struct _CoefficientView[T: Copyable & Deinitable, origin: ImmOrigin](Copyable, Movable, Sized, Iterable):
    # Interior origins invalidate the view on owner mutation, including changes
    # that free/reallocate nested lists. Whole-object origins are insufficient:
    # the owner can remain alive after its coefficient storage has changed.
    var _flat: Pointer[List[Self.T], Self.origin]
    var _chunks: Optional[Pointer[List[List[Self.T]], Self.origin]]
    var _count: Int
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _ViewIterator[Self.T, Self.origin]

    def __len__(self) -> Int:
        return self._count

    def __getitem__(self, index: Int) -> ref[Self.origin] Self.T:
        debug_assert(index >= 0 and index < self._count, "coefficient index out of bounds")
        if self._chunks:
            var chunk = 0 if index < _FIRST else 1 if index < _PREFIX else 2 + (index - _PREFIX) // _BLOCK
            var offset = index if index < _FIRST else index - _FIRST if index < _PREFIX else (index - _PREFIX) % _BLOCK
            return Pointer(to=self._chunks.value()[][chunk][offset]).unsafe_origin_cast[Self.origin]()[]
        return Pointer(to=self._flat[][index]).unsafe_origin_cast[Self.origin]()[]

    def _run(self, index: Int) -> Span[Self.T, Self.origin]:
        # Return the remainder of the containing allocation. The view's origin
        # prevents structural mutation until the borrowed run is no longer used.
        debug_assert(index >= 0 and index < self._count, "coefficient index out of bounds")
        if self._chunks:
            var chunk = 0 if index < _FIRST else 1 if index < _PREFIX else 2 + (index - _PREFIX) // _BLOCK
            var offset = index if index < _FIRST else index - _FIRST if index < _PREFIX else (index - _PREFIX) % _BLOCK
            ref values = self._chunks.value()[][chunk]
            return Span(unsafe_ptr=Pointer(to=values[offset]).unsafe_origin_cast[Self.origin](), length=len(values) - offset)
        return Span(unsafe_ptr=Pointer(to=self._flat[][index]).unsafe_origin_cast[Self.origin](), length=self._count - index)

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        return _ViewIterator[Self.T, Self.origin](self.copy(), 0)


@fieldwise_init
struct _ViewIterator[T: Copyable & Deinitable, origin: ImmOrigin](Iterator):
    comptime Element = Self.T
    var view: _CoefficientView[Self.T, Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> ref[Self.origin] Self.T:
        if self.position == len(self.view):
            raise StopIteration()
        self.position += 1
        return self.view[self.position - 1]

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        var remaining = len(self.view) - self.position
        return remaining, Optional(remaining)


@__unsafe_nested_origins_read_only
def _view[T: Copyable & Deinitable](values: List[T]) -> _CoefficientView[T, origin_of(values)._get_owned_interior["coefficient_view"]]:
    comptime view_origin = origin_of(values)._get_owned_interior["coefficient_view"]
    return _CoefficientView[T, view_origin](Pointer(to=Pointer(to=values)._get_ref_with_unsafe_interior_origin["coefficient_view", origin_of(values)]()), None, len(values))


@__unsafe_nested_origins_read_only
def _view[T: Copyable & Deinitable](values: CoefficientBuffer[T]) -> _CoefficientView[T, origin_of(values)._get_owned_interior["coefficient_view"]]:
    comptime view_origin = origin_of(values)._get_owned_interior["coefficient_view"]
    ref owner = Pointer(to=values)._get_ref_with_unsafe_interior_origin["coefficient_view", origin_of(values)]()
    var chunks = Optional[Pointer[List[List[T]], view_origin]]()
    if len(owner._chunks):
        chunks = Pointer(to=owner._chunks).unsafe_origin_cast[view_origin]()
    return _CoefficientView[T, view_origin](Pointer(to=owner._flat).unsafe_origin_cast[view_origin](), chunks^, len(owner))


@__unsafe_nested_origins_read_only
def _view[T: Copyable & Deinitable, origin: ImmOrigin](values: _CoefficientView[T, origin]) -> _CoefficientView[T, origin]:
    return values.copy()
