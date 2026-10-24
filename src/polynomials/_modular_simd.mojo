"""Exact native SIMD arithmetic for modular Bezout and cofactor interpolation.

Ordinary residues satisfy 2 <= p < 2^32. Products fit UInt64; exact
reciprocal reduction retains their high words with native integer lanes.
"""
from ._word_simd import _mul_high, _WORD_LANES


@always_inline
def _mod_add[W: Int](a: SIMD[DType.uint64, W], b: SIMD[DType.uint64, W], p: UInt64) -> SIMD[DType.uint64, W]:
    var t = a + b
    return t.ge(p).select(t - p, t)


@always_inline
def _mod_sub[W: Int](a: SIMD[DType.uint64, W], b: SIMD[DType.uint64, W], p: UInt64) -> SIMD[DType.uint64, W]:
    return a.ge(b).select(a - b, a + p - b)


@fieldwise_init
struct _Prime32(ImplicitlyCopyable, Movable):
    """Ordinary residues modulo 2 <= p < 2^32, reduced by a fixed reciprocal."""
    var p: UInt64
    var reciprocal: UInt64

    @staticmethod
    def of(p: UInt64) -> Self:
        return Self(p, UInt64((UInt128(1) << 64) // UInt128(p)))

    @always_inline
    def mul[W: Int](self, a: SIMD[DType.uint64, W], b: SIMD[DType.uint64, W]) -> SIMD[DType.uint64, W]:
        """Canonical product, with every lane of a,b below p."""
        # Masks expose the 32-by-32 product to instruction selection.
        var product = (a & 0xFFFFFFFF) * (b & 0xFFFFFFFF)
        var q: SIMD[DType.uint64, W]
        comptime if W == 1:
            # Explicit unsigned widening also covers products with bit 63 set.
            q = UInt64((product.cast[DType.uint128]() * self.reciprocal.cast[DType.uint128]()) >> 64)
        else:
            q = _mul_high[W](product, SIMD[DType.uint64, W](self.reciprocal))
        # q is floor(product/p) or one less, since product < 2^64.
        # Therefore product-q*p is nonnegative and below 2p < 2^33.
        var r = product - q * self.p
        return r.ge(self.p).select(r - self.p, r)

    @always_inline
    def power[W: Int](self, a: SIMD[DType.uint64, W], exponent: UInt64) -> SIMD[DType.uint64, W]:
        var result = SIMD[DType.uint64, W](1)
        var base = a
        var e = exponent
        while e:
            if e & 1:
                result = self.mul[W](result, base)
            e >>= 1
            if e:
                base = self.mul[W](base, base)
        return result


def _inverse32[vectorized: Bool = True](mut values: List[UInt64], field: _Prime32):
    """Invert independent nonzero residues modulo a prime, in place."""
    var i = 0
    var data = values.unsafe_ptr()
    comptime if vectorized:
        comptime W = _WORD_LANES
        while len(values) - i >= W:
            var row = data.unsafe_offset(i)
            row.unsafe_store[width=W](field.power[W](row.unsafe_load[width=W](), field.p - 2))
            i += W
    while i < len(values):
        values[i] = field.power[1](values[i], field.p - 2)
        i += 1


def _subtract32[vectorized: Bool = True](destination: Pointer[mut=True, UInt64, MutUntrackedOrigin], source: Pointer[mut=False, UInt64, _], count: Int, factor: UInt64, field: _Prime32):
    """out -= factor*source modulo p; exact aliasing or disjoint buffers."""
    var i = 0
    comptime if vectorized:
        comptime W = _WORD_LANES
        while count - i >= W:
            var t = field.mul[W](source.unsafe_offset(i).unsafe_load[width=W](), SIMD[DType.uint64, W](factor))
            var target = destination.unsafe_offset(i)
            target.unsafe_store[width=W](_mod_sub[W](target.unsafe_load[width=W](), t, field.p))
            i += W
    while i < count:
        var target = destination.unsafe_offset(i)
        target.unsafe_store(_mod_sub[1](target.unsafe_load(), field.mul[1](source.unsafe_offset(i).unsafe_load(), factor), field.p))
        i += 1
