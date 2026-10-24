"""Exact SIMD arithmetic on independent native-word coefficients."""

comptime _WORD_LANES = 8


@always_inline
def _mul_wide[W: Int](a: SIMD[DType.uint64, W], b: SIMD[DType.uint64, W]) -> Tuple[SIMD[DType.uint64, W], SIMD[DType.uint64, W]]:
    """The low and high words of each exact unsigned 64-by-64-bit product.

    Every partial product has 32-bit operands. The middle sum is below
    3 * 2^32, and the final high word fits by the 128-bit product bound.
    """
    var mask = SIMD[DType.uint64, W](0xFFFFFFFF)
    var al = a & mask
    var bl = b & mask
    var ah = a >> 32
    var bh = b >> 32
    var ll = al * bl
    var lh = al * bh
    var hl = ah * bl
    var middle = (ll >> 32) + (lh & mask) + (hl & mask)
    return ((ll & mask) | (middle << 32), ah * bh + (lh >> 32) + (hl >> 32) + (middle >> 32))


@always_inline
def _mul_high[W: Int](a: SIMD[DType.uint64, W], b: SIMD[DType.uint64, W]) -> SIMD[DType.uint64, W]:
    """The high word of each exact unsigned 64-by-64-bit product."""
    return _mul_wide[W](a, b)[1]


@always_inline
def _sum_products(a: Pointer[mut=False, UInt64, _], b: Pointer[mut=False, UInt64, _], count: Int) -> Tuple[UInt128, UInt64]:
    """The exact 192-bit sum of products of count pairs of words below 2^62.

    Sixteen products per lane sum below 2^128. Between blocks, fold lanes
    into a 128-bit low part and a carry word. Since count fits Int, the
    complete sum is below 2^187, so the carry word cannot overflow.
    """
    comptime W = _WORD_LANES
    comptime Lane = SIMD[DType.uint64, W]
    var low = UInt128(0)
    var carry = UInt64(0)
    var i = 0
    while count - i >= W:
        var stop = i + min(count - i, 16 * W)
        var lo = Lane(0)
        var hi = Lane(0)
        while stop - i >= W:
            var product = _mul_wide[W](a.unsafe_offset(i).unsafe_load[width=W](), b.unsafe_offset(i).unsafe_load[width=W]())
            var next = lo + product[0]
            hi += product[1] + next.lt(lo).cast[DType.uint64]()
            lo = next
            i += W
        # Sum base-2^32 limbs horizontally: eight limbs and their incoming
        # carry fit comfortably in a word. This avoids extracting every
        # 128-bit lane separately while preserving all carries between them.
        var mask = Lane(0xFFFFFFFF)
        var limb0 = (lo & mask).reduce_add()
        var limb1 = (lo >> 32).reduce_add() + (limb0 >> 32)
        var limb2 = (hi & mask).reduce_add() + (limb1 >> 32)
        var limb3 = (hi >> 32).reduce_add() + (limb2 >> 32)
        var bottom = (limb0 & 0xFFFFFFFF) | (limb1 << 32)
        var top = (limb2 & 0xFFFFFFFF) | (limb3 << 32)
        var next = low + ((UInt128(top) << 64) | UInt128(bottom))
        carry += (limb3 >> 32) + (UInt64(1) if next < low else UInt64(0))
        low = next
    while i < count:
        var next = low + UInt128(a.unsafe_offset(i).unsafe_load()) * UInt128(b.unsafe_offset(i).unsafe_load())
        carry += UInt64(1) if next < low else UInt64(0)
        low = next
        i += 1
    return (low, carry)


@always_inline
def _column_axpy[vectorized: Bool = True](data: Pointer[mut=True, UInt64, _], column: Pointer[mut=False, UInt32, _], count: Int, c: UInt64, square: UInt64):
    """Add c times a column, keeping accumulators below square = p*p.

    Requires p*p < 2^63, c and column entries below p, and disjoint buffers.
    Each intermediate is below 2*p*p < 2^64. This includes the primes
    just above 2^31 used by rational FGLM. The caller charges the work.
    """
    var i = 0
    comptime if vectorized:
        comptime W = _WORD_LANES
        var factor = SIMD[DType.uint64, W](c) & SIMD[DType.uint64, W](0xFFFFFFFF)
        var bound = SIMD[DType.uint64, W](square)
        while i + W <= count:
            var values = column.unsafe_offset(i).unsafe_load[width=W]().cast[DType.uint64]()
            var x = data.unsafe_offset(i).unsafe_load[width=W]() + factor * values
            data.unsafe_offset(i).unsafe_store[width=W](x.ge(bound).select(x - bound, x))
            i += W
    while i < count:
        var x = data.unsafe_offset(i).unsafe_load() + c * UInt64(column.unsafe_offset(i).unsafe_load())
        data.unsafe_offset(i).unsafe_store(x - square if x >= square else x)
        i += 1
