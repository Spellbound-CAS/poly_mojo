"""Packed-exponent products and exact quotients of sparse ZZ and GF(p) polynomials in LEX order.

A monomial over variables 0..n-1 becomes one machine word by Kronecker
substitution with mixed radices: radix_i exceeds the largest exponent of
variable i that the result can hold, and variable 0 is the most significant
digit. Within those bounds the map is additive (monomial products are key
sums, with no carries) and increasing (LEX comparison is integer comparison).

Products split the key into a prefix of the first L variables and a low part
of the rest, with L the smallest count whose low range fits a dense window.
Each output prefix collects the pairs of operand groups whose prefixes sum to
it, accumulates their products in the dense window, and emits the window's
touched slots in increasing order, so the output arrives sorted without a
heap or a sort of the terms. ZZ coefficients of at most 63 bits accumulate in
192-bit two's complement words; GF(p) residues accumulate in 128-bit words.

Exact quotients use the same chunks from the top prefix down: with h the
divisor's top prefix, the dividend chunk at prefix t minus the known products
of earlier quotient chunks equals the quotient chunk at t - h times the
divisor's top chunk. When that top chunk is a single term, each quotient
chunk follows term by term.

Every choice depends on sizes and exponents only, never on time.
"""
from apn_mojo import Integer
from .coefficients import CoefficientRing, ZZ, ModInt, _mod_mul
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _product, _sum
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term
from std.collections import Dict
from std.math import log2
from ._packed_poly import _PackedPoly

comptime _WINDOW = 1 << 18
"""The largest dense accumulation window, in slots."""
comptime _PRODUCT_WINDOW = 1 << 17
"""The largest product window: 2 MiB of 128-bit slots, near the L2 size, so a
window's random adds stay cached; a box of a quarter-million slots or more
splits on its first variable. Division keeps _WINDOW."""
comptime _POWER_GATE = 1
"""The power recurrence runs when its heap steps stay below this many times the products of the last square."""
comptime _POWER_SPREAD = UInt64(2)
"""The largest base exponent for which the power recurrence runs."""
comptime _KEY_LIMIT = UInt64(1) << 62
"""Keys stay below this, so prefix and window arithmetic never wraps."""


struct _Packing(Copyable, Movable):
    """Mixed radices over variables 0..n-1 and the prefix/window split."""
    var n: Int
    var radix: List[UInt64]
    var weight: List[UInt64]
    """weight[i] = radix[i+1] * ... * radix[n-1]: the place value of variable i."""
    var prefix_count: Int
    """L, the number of leading variables in the prefix."""
    var window: UInt64
    """The low range: radix[L] * ... * radix[n-1]."""

    def __init__(out self, var radix: List[UInt64], limit: Int = _WINDOW):
        self.n = len(radix)
        self.weight = List[UInt64](length=self.n, fill=1)
        var i = self.n - 2
        while i >= 0:
            self.weight[i] = self.weight[i + 1] * radix[i + 1]
            i -= 1
        self.prefix_count = 0
        self.window = (self.weight[0] * radix[0]) if self.n else 1
        while self.prefix_count < self.n and self.window > UInt64(limit):
            self.window = self.weight[self.prefix_count]
            self.prefix_count += 1
        self.radix = radix^

    def split(mut self, prefix_count: Int):
        """Use the first prefix_count variables as the prefix (at least the current count)."""
        while self.prefix_count < prefix_count and self.prefix_count < self.n:
            self.window = self.weight[self.prefix_count]
            self.prefix_count += 1

    def key(self, m: Monomial) -> UInt64:
        var k = UInt64(0)
        for pair in m.pairs:
            k += UInt64(pair[1]) * self.weight[Int(pair[0])]
        return k

    def monomial(self, key: UInt64, mut budget: Budget) raises PolyError -> Monomial:
        var count = 0
        var rest = key
        for i in range(self.n):
            var e = rest // self.weight[i]
            rest -= e * self.weight[i]
            count += Int(e != 0)
        _allocation[Tuple[Int64, UInt32]](count, budget)
        var pairs = List[Tuple[Int64, UInt32]](capacity=count)
        rest = key
        for i in range(self.n):
            var e = rest // self.weight[i]
            rest -= e * self.weight[i]
            if e:
                pairs.append((Int64(i), UInt32(e)))
        var result = Monomial()
        result.pairs = pairs^
        return result^

    def exponent(self, key: UInt64, i: Int) -> UInt64:
        return (key // self.weight[i]) % self.radix[i]


def _lex_identity[C: CoefficientRing](p: SparsePoly[C]) -> Bool:
    """Whether p is in LEX over variables 0..n-1 in priority order, using only those."""
    if p.order.kind != MonomialOrder.LEX or not p.order._identity():
        return False
    for term in p.terms:
        if not p.order._within(term.monomial):
            return False
    return True


def _maxima[C: CoefficientRing](p: SparsePoly[C], n: Int) -> List[UInt64]:
    var result = List[UInt64](length=n, fill=0)
    for term in p.terms:
        for pair in term.monomial.pairs:
            result[Int(pair[0])] = max(result[Int(pair[0])], UInt64(pair[1]))
    return result^


def _radices(a: List[UInt64], b: List[UInt64]) -> Optional[List[UInt64]]:
    """radix_i = a_i + b_i + 1, or None when the key range reaches _KEY_LIMIT."""
    var radix = List[UInt64](capacity=len(a))
    var total = UInt128(1)
    for i in range(len(a)):
        var r = a[i] + b[i] + 1
        total *= UInt128(r)
        if total >= UInt128(_KEY_LIMIT):
            return None
        radix.append(r)
    return radix^


# ----------------------------------------------------------- coefficients

@always_inline
def _small[C: CoefficientRing](c: C) -> Bool:
    comptime if C == ZZ:
        return rebind[ZZ](c).value.magnitude_bit_length() <= 63
    else:
        return True


def _signed[C: CoefficientRing](c: C) -> Int64:
    """A ZZ coefficient of at most 63 bits, or a residue (as its bits)."""
    comptime if C == ZZ:
        try:
            return rebind[ZZ](c).value.to_native_exact[DType.int64]()
        except:
            return 0
    elif C == ModInt:
        return Int64(rebind[ModInt](c).value)
    else:
        return 0


struct _Accumulator[narrow: Bool = False](Movable):
    """Dense window slots.

    ZZ slots hold two's complement sums: 128 bits in `lows`, extended to 192
    bits by `highs` when the coefficient bound needs it (`wide`), or 64 bits
    in `words` alone when the bound stays below 2^63 (`narrow`, ZZ only).
    GF(p) slots hold unreduced 128-bit sums. Tracking records the touched
    slots of sparse prefixes; dense prefixes scan the window instead.
    """
    var lows: List[UInt128]
    var highs: List[UInt64]
    var words: List[Int64]
    var touched_flag: List[UInt8]
    var touched: List[Int]

    def __init__(out self, slots: Int, wide: Bool):
        comptime if Self.narrow:
            self.lows = List[UInt128]()
            self.words = List[Int64](length=slots, fill=0)
        else:
            self.lows = List[UInt128](length=slots, fill=0)
            self.words = List[Int64]()
        self.highs = List[UInt64](length=slots if wide else 0, fill=0)
        self.touched_flag = List[UInt8](length=slots, fill=0)
        self.touched = List[Int]()

    @always_inline
    def touch(mut self, slot: Int):
        if not self.touched_flag.unsafe_get(slot):
            self.touched_flag.unsafe_get(slot) = 1
            self.touched.append(slot)

    @always_inline
    def add_signed[wide: Bool, track: Bool](mut self, slot: Int, a: Int64, b: Int64):
        """Add a * b for |a|, |b| < 2^63."""
        comptime if track:
            self.touch(slot)
        comptime if Self.narrow:
            # The coefficient bound keeps every partial sum below 2^63.
            self.words.unsafe_get(slot) += a * b
        else:
            var product = Int128(a) * Int128(b)
            self._add[wide](slot, UInt128(product), UInt64(0) - UInt64(Int(product < 0)))

    @always_inline
    def _add[wide: Bool](mut self, slot: Int, value: UInt128, extension: UInt64):
        """Add the two's complement value whose bits above 128 are all `extension`."""
        var old = self.lows.unsafe_get(slot)
        var new = old + value
        self.lows.unsafe_get(slot) = new
        comptime if wide:
            self.highs.unsafe_get(slot) += extension + UInt64(Int(new < old))

    def add_integer[wide: Bool](mut self, slot: Int, value: Integer) raises PolyError:
        """Add a ZZ coefficient below 2^127 in magnitude, or below 2^190 in a wide slot."""
        self.touch(slot)
        var bits = value.magnitude_bit_length()
        if bits <= 63:
            var word: Int64
            try:
                word = value.to_native_exact[DType.int64]()
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            self._add[wide](slot, UInt128(Int128(word)), UInt64(0) - UInt64(Int(word < 0)))
            return
        var bytes = value.to_bytes()
        if len(bytes) > (24 if wide else 16):
            raise PolyError(PolyErrorKind.CAPACITY, "packed accumulator coefficient overflow")
        var low = UInt128(0)
        var high = UInt64(0)
        for i in range(len(bytes)):
            if i < 16:
                low |= UInt128(bytes[i]) << UInt128(8 * i)
            else:
                high |= UInt64(bytes[i]) << UInt64(8 * (i - 16))
        if value.sign() < 0:
            # 192-bit two's complement negation.
            low = ~low + 1
            high = ~high + UInt64(Int(low == 0))
        self._add[wide](slot, low, high)

    @always_inline
    def add_modular[track: Bool](mut self, slot: Int, a: UInt64, b: UInt64, p: UInt64):
        """Add a * b for residues a, b < p < 2^62, reducing before the sum could overflow."""
        comptime if track:
            self.touch(slot)
        var acc = self.lows.unsafe_get(slot)
        if acc >= UInt128(1) << 127:
            acc %= UInt128(p)
        self.lows.unsafe_get(slot) = acc + UInt128(a) * UInt128(b)

    def modular(self, slot: Int, p: UInt64) -> UInt64:
        return UInt64(self.lows.unsafe_get(slot) % UInt128(p))

    def signed_words(self, slot: Int) -> Tuple[UInt64, UInt64, UInt64]:
        comptime if Self.narrow:
            var w = self.words.unsafe_get(slot)
            var extension = UInt64(0) - UInt64(Int(w < 0))
            return (UInt64(w), extension, extension)
        var low = self.lows.unsafe_get(slot)
        var hi = UInt64(low >> 64)
        var top = self.highs.unsafe_get(slot) if len(self.highs) else UInt64(0) - (hi >> 63)
        return (UInt64(low), hi, top)

    @always_inline
    def empty(self, slot: Int) -> Bool:
        """Whether the slot's sum is zero (before any modular reduction)."""
        comptime if Self.narrow:
            return self.words.unsafe_get(slot) == 0
        return self.lows.unsafe_get(slot) == 0 and (not len(self.highs) or self.highs.unsafe_get(slot) == 0)

    def clear(mut self, slot: Int):
        comptime if Self.narrow:
            self.words.unsafe_get(slot) = 0
            self.touched_flag.unsafe_get(slot) = 0
            return
        self.lows.unsafe_get(slot) = 0
        if len(self.highs):
            self.highs.unsafe_get(slot) = 0
        self.touched_flag.unsafe_get(slot) = 0

    def sorted_touched(mut self) -> List[Int]:
        var result = List[Int]()
        swap(result, self.touched)
        sort(result)
        return result^


def _is_zero(words: Tuple[UInt64, UInt64, UInt64]) -> Bool:
    return words[0] == 0 and words[1] == 0 and words[2] == 0


def _integer(words: Tuple[UInt64, UInt64, UInt64]) -> Integer:
    """The value of 192-bit two's complement words."""
    var negative = (words[2] >> 63) == 1
    var w0 = words[0]
    var w1 = words[1]
    var w2 = words[2]
    if negative:
        # Negate: invert and add one.
        w0 = ~w0
        w1 = ~w1
        w2 = ~w2
        w0 += 1
        if w0 == 0:
            w1 += 1
            if w1 == 0:
                w2 += 1
    if w1 == 0 and w2 == 0:
        if w0 <= UInt64(Int64.MAX):
            var v = Int64(w0)
            return Integer(-v if negative else v)
        if not negative:
            return Integer(w0)
    # One magnitude allocation, as _kronecker builds its coefficients.
    var count = 6 if w2 else 4 if w1 else 2
    var magnitude = List[UInt32](capacity=count)
    for w in [w0, w1, w2]:
        if len(magnitude) < count:
            magnitude.append(UInt32(w & 0xFFFFFFFF))
            magnitude.append(UInt32(w >> 32))
    return Integer._from_words(magnitude^, negative)


def _coefficient[C: CoefficientRing, narrow: Bool](acc: _Accumulator[narrow], slot: Int, template: C) -> C:
    comptime if C == ZZ and narrow:
        return rebind_var[C](ZZ(Integer(acc.words.unsafe_get(slot))))
    elif C == ZZ:
        return rebind_var[C](ZZ(_integer(acc.signed_words(slot))))
    elif C == ModInt:
        ref t = rebind[ModInt](template)
        return rebind_var[C](t._with_residue(acc.modular(slot, t.modulus)))
    else:
        return template.copy()


def _slot_zero[C: CoefficientRing, narrow: Bool](acc: _Accumulator[narrow], slot: Int, template: C) -> Bool:
    comptime if C == ZZ:
        return _is_zero(acc.signed_words(slot))
    elif C == ModInt:
        return acc.modular(slot, rebind[ModInt](template).modulus) == 0
    else:
        return True


# ---------------------------------------------------------------- groups

@fieldwise_init
struct _Packed(Movable):
    """Terms as increasing keys with word coefficients, and their groups by prefix."""
    var keys: List[UInt64]
    var lows: List[Int]
    """Each key's window part."""
    var values: List[Int64]
    var group_prefix: List[UInt64]
    var group_start: List[Int]
    """Group g holds terms group_start[g] ..< group_start[g + 1]."""


def _pack[C: CoefficientRing](p: SparsePoly[C], packing: _Packing, mut budget: Budget) raises PolyError -> _Packed:
    var count = len(p.terms)
    _allocation[UInt64](_product(3, count), budget)
    var keys = List[UInt64](capacity=count)
    var lows = List[Int](capacity=count)
    var values = List[Int64](capacity=count)
    var group_prefix = List[UInt64]()
    var group_start = List[Int]()
    for i in range(count):
        budget.charge(1)
        ref term = p.terms[i]
        var key = packing.key(term.monomial)
        var prefix = key // packing.window
        if not len(group_prefix) or group_prefix[len(group_prefix) - 1] != prefix:
            group_prefix.append(prefix)
            group_start.append(i)
        keys.append(key)
        lows.append(Int(key - prefix * packing.window))
        values.append(_signed(term.coefficient))
    group_start.append(count)
    return _Packed(keys^, lows^, values^, group_prefix^, group_start^)


def _pack_keys[C: CoefficientRing](var keys: List[UInt64], coefficients: List[C], packing: _Packing, mut budget: Budget) raises PolyError -> _Packed:
    """Increasing keys in the packing's radix with their coefficients."""
    _allocation[Int64](len(keys), budget)
    var values = List[Int64](capacity=len(keys))
    for c in coefficients:
        budget.charge(1)
        values.append(_signed(c))
    return _grouped(keys^, values^, packing, budget)


def _grouped(var keys: List[UInt64], var values: List[Int64], packing: _Packing, mut budget: Budget) raises PolyError -> _Packed:
    var count = len(keys)
    _allocation[UInt64](_product(2, count), budget)
    var lows = List[Int](capacity=count)
    var group_prefix = List[UInt64]()
    var group_start = List[Int]()
    for i in range(count):
        var key = keys[i]
        var prefix = key // packing.window
        if not len(group_prefix) or group_prefix[len(group_prefix) - 1] != prefix:
            group_prefix.append(prefix)
            group_start.append(i)
        lows.append(Int(key - prefix * packing.window))
    group_start.append(count)
    return _Packed(keys^, lows^, values^, group_prefix^, group_start^)


def _usable[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C]) -> Bool:
    comptime if not (C == ZZ or C == ModInt):
        return False
    if not len(a.order.variables):
        return False
    if not _lex_identity(a) or not _lex_identity(b) or a.order != b.order:
        return False
    for term in a.terms:
        if not _small(term.coefficient):
            return False
    for term in b.terms:
        if not _small(term.coefficient):
            return False
    return True


def _packed_product[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError -> Optional[SparsePoly[C]]:
    """The canonical product, or None when the operands do not fit the packed kernel.

    Both operands are nonzero, in the same ring and order.
    """
    if not _usable(a, b):
        return None
    var n = len(a.order.variables)
    var radix = _radices(_maxima(a, n), _maxima(b, n))
    if not radix:
        return None
    # A window larger than eight slots per term pair would cost more to clear than to use.
    var packing = _Packing(radix.take(), min(_PRODUCT_WINDOW, max(1024, 8 * len(a.terms) * len(b.terms))))
    var left = _pack(a, packing, budget)
    var right = _pack(b, packing, budget)
    # Squares visit each unordered pair of terms once, doubling the off-diagonal products.
    var square = a == b
    comptime if C == ZZ:
        for term in a.terms:
            if rebind[ZZ](term.coefficient).value.magnitude_bit_length() > 62:
                square = False
                break
    # Sums need 192 bits only when the coefficient bound reaches 2^127.
    var wide = False
    var narrow = False
    comptime if C == ZZ:
        var bits = _bits(a) + _bits(b) + _bit_length(min(len(a.terms), len(b.terms))) + 1
        wide = bits >= 127
        var top_a = UInt128(0)
        var top_b = UInt128(0)
        var norm_a = UInt128(0)
        var norm_b = UInt128(0)
        for term in a.terms:
            var v = UInt128(abs(_signed(term.coefficient)))
            top_a = max(top_a, v)
            norm_a += v
        for term in b.terms:
            var v = UInt128(abs(_signed(term.coefficient)))
            top_b = max(top_b, v)
            norm_b += v
        narrow = _fits_word(top_a, norm_a, top_b, norm_b, square)
    var sink = _Sink[C, False]()
    if not _product_core(left, right, packing, square, wide, narrow, a.terms[0].coefficient, sink, budget):
        return None
    var result = SparsePoly[C].zero(a.order, a.ring, budget)
    swap(result.terms, sink.terms)
    return result^


def _packed_poly_product[C: CoefficientRing](a: _PackedPoly[C], b: _PackedPoly[C], mut budget: Budget) raises PolyError -> Optional[_PackedPoly[C]]:
    """The product of packed operands, packed in its own radix, or None when the kernel does not apply."""
    comptime if not (C == ZZ or C == ModInt):
        return None
    if not len(a) or not len(b) or a.order != b.order or not len(a.radix):
        return None
    for v in a.values:
        if not _small(v):
            return None
    for v in b.values:
        if not _small(v):
            return None
    var top_a = List[UInt64](capacity=len(a.radix))
    var top_b = List[UInt64](capacity=len(b.radix))
    for r in a.radix:
        top_a.append(r - 1)
    for r in b.radix:
        top_b.append(r - 1)
    var radix = _radices(top_a, top_b)
    if not radix:
        return None
    var packing = _Packing(radix.value().copy(), min(_PRODUCT_WINDOW, max(1024, 8 * len(a) * len(b))))
    var left = _pack_keys(a.repacked(packing.radix), a.values, packing, budget)
    var right = _pack_keys(b.repacked(packing.radix), b.values, packing, budget)
    var square = a == b
    var narrow = False
    comptime if C == ZZ:
        for v in a.values:
            if rebind[ZZ](v).value.magnitude_bit_length() > 62:
                square = False
                break
    var wide = False
    comptime if C == ZZ:
        var bits_a = 0
        for v in a.values:
            bits_a = max(bits_a, rebind[ZZ](v).value.magnitude_bit_length())
        var bits_b = 0
        for v in b.values:
            bits_b = max(bits_b, rebind[ZZ](v).value.magnitude_bit_length())
        wide = bits_a + bits_b + _bit_length(min(len(a), len(b))) + 1 >= 127
        var top_a = UInt128(0)
        var norm_a = UInt128(0)
        for v in a.values:
            var m = UInt128(abs(_signed(v)))
            top_a = max(top_a, m)
            norm_a += m
        var top_b = UInt128(0)
        var norm_b = UInt128(0)
        for v in b.values:
            var m = UInt128(abs(_signed(v)))
            top_b = max(top_b, m)
            norm_b += m
        narrow = _fits_word(top_a, norm_a, top_b, norm_b, square)
    var sink = _Sink[C, True]()
    if not _product_core(left, right, packing, square, wide, narrow, a.values[0], sink, budget):
        return None
    var keys = List[UInt64]()
    var values = List[C]()
    swap(keys, sink.keys)
    swap(values, sink.values)
    return _PackedPoly[C](a.order.copy(), a.ring.copy(), radix.take(), keys^, values^)


def _product_core[C: CoefficientRing, packed: Bool](left: _Packed, right: _Packed, packing: _Packing, square: Bool, wide: Bool, narrow: Bool, template: C, mut sink: _Sink[C, packed], mut budget: Budget) raises PolyError -> Bool:
    """Accumulate the product of packed operands prefix by prefix into the sink; False when the kernel does not apply."""
    var groups_a = len(left.group_prefix)
    var groups_b = len(right.group_prefix)
    var pair_count = _product(groups_a, groups_b)
    # Fragmented supports would spend more on group pairs than on products.
    if pair_count > _product(len(left.keys), len(right.keys)) // 4 + 16:
        return False
    _allocation[UInt64](pair_count, budget)
    var pairs = List[UInt64](capacity=pair_count)
    var shift = 0
    while (1 << shift) < pair_count:
        shift += 1
    # Sort (output prefix, pair index) by packing both into one word when they fit.
    var prefix_bits = 0
    var prefix_range = left.group_prefix[groups_a - 1] + right.group_prefix[groups_b - 1] + 1
    while (UInt64(1) << UInt64(prefix_bits)) < prefix_range:
        prefix_bits += 1
    if prefix_bits + shift > 63:
        return False
    for i in range(groups_a):
        for j in range(i if square else 0, groups_b):
            budget.charge(1)
            var prefix = left.group_prefix[i] + right.group_prefix[j]
            pairs.append((prefix << UInt64(shift)) | UInt64(i * groups_b + j))
    pair_count = len(pairs)
    _charge_size(pair_count, budget)
    sort(pairs)
    _allocation[UInt64](_product(4, Int(packing.window)), budget)
    comptime if C == ZZ:
        if narrow:
            _product_windows[C, packed, True](pairs, shift, groups_b, left, right, packing, square, wide, template, sink, budget)
            return True
    _product_windows[C, packed, False](pairs, shift, groups_b, left, right, packing, square, wide, template, sink, budget)
    return True


def _product_windows[C: CoefficientRing, packed: Bool, narrow: Bool](pairs: List[UInt64], shift: Int, groups_b: Int, left: _Packed, right: _Packed, packing: _Packing, square: Bool, wide: Bool, template: C, mut sink: _Sink[C, packed], mut budget: Budget) raises PolyError:
    """Accumulate the sorted group pairs window by window into the sink."""
    var pair_count = len(pairs)
    var window = Int(packing.window)
    var acc = _Accumulator[narrow](window, wide)
    var p = UInt64(0)
    comptime if C == ModInt:
        p = rebind[ModInt](template).modulus
    var mask = (UInt64(1) << UInt64(shift)) - 1
    var position = 0
    while position < pair_count:
        var prefix = pairs[position] >> UInt64(shift)
        var end = position
        var work = 0
        while end < pair_count and (pairs[end] >> UInt64(shift)) == prefix:
            var index = Int(pairs[end] & mask)
            work += (left.group_start[index // groups_b + 1] - left.group_start[index // groups_b]) * (right.group_start[index % groups_b + 1] - right.group_start[index % groups_b])
            end += 1
        _charge_size(work, budget)
        # A prefix with work comparable to the window scans it instead of tracking slots.
        var dense = 4 * work >= window
        while position < end:
            var index = Int(pairs[position] & mask)
            position += 1
            var i = index // groups_b
            var j = index % groups_b
            var a0 = left.group_start[i]
            var a1 = left.group_start[i + 1]
            var b0 = right.group_start[j]
            var b1 = right.group_start[j + 1]
            var diagonal = square and i == j
            var doubled = square and i != j
            if dense:
                if wide:
                    _accumulate[C, True, False, narrow](acc, left, right, a0, a1, b0, b1, p, diagonal, doubled)
                else:
                    _accumulate[C, False, False, narrow](acc, left, right, a0, a1, b0, b1, p, diagonal, doubled)
            elif wide:
                _accumulate[C, True, True, narrow](acc, left, right, a0, a1, b0, b1, p, diagonal, doubled)
            else:
                _accumulate[C, False, True, narrow](acc, left, right, a0, a1, b0, b1, p, diagonal, doubled)
        if dense:
            _emit_window(acc, packing, prefix, template, sink, budget)
        else:
            _emit(acc, packing, prefix, template, sink, budget)


def _fits_word(top_a: UInt128, norm_a: UInt128, top_b: UInt128, norm_b: UInt128, square: Bool) -> Bool:
    """Whether every partial sum of a product stays below 2^63 in magnitude.

    A coefficient's partial sums are bounded by sums of |a_i| |b_j| over its
    pairs, so by max|a| * |b|_1 and by max|b| * |a|_1 (top_* the largest
    magnitudes, norm_* the 1-norms, all of 63-bit coefficients). Squares,
    which add doubled off-diagonal products, keep one bit in hand.
    """
    var bound = min(_bits128(top_a) + _bits128(norm_b), _bits128(top_b) + _bits128(norm_a))
    return bound + Int(square) <= 63


def _bits128(x: UInt128) -> Int:
    var bits = 0
    while bits < 128 and (x >> UInt128(bits)) != 0:
        bits += 1
    return bits


def _bits[C: CoefficientRing](p: SparsePoly[C]) -> Int:
    var bits = 0
    for term in p.terms:
        bits = max(bits, rebind[ZZ](term.coefficient).value.magnitude_bit_length())
    return bits


def _bit_length(n: Int) -> Int:
    var bits = 0
    while (1 << bits) <= n:
        bits += 1
    return bits


@always_inline
def _twice[C: CoefficientRing](v: Int64, p: UInt64) -> Int64:
    comptime if C == ZZ:
        return 2 * v
    else:
        var t = 2 * UInt64(v)
        return Int64(t - p if t >= p else t)


def _accumulate[C: CoefficientRing, wide: Bool, track: Bool, narrow: Bool](mut acc: _Accumulator[narrow], left: _Packed, right: _Packed, a0: Int, a1: Int, b0: Int, b1: Int, p: UInt64, diagonal: Bool, doubled: Bool):
    """Add the products of terms a0 ..< a1 and b0 ..< b1; for a diagonal group of a square, each unordered pair once."""
    for x in range(a0, a1):
        var ka = left.lows.unsafe_get(x)
        var va = left.values.unsafe_get(x)
        var start = b0
        if diagonal:
            comptime if C == ZZ:
                acc.add_signed[wide, track](ka + ka, va, va)
            else:
                acc.add_modular[track](ka + ka, UInt64(va), UInt64(va), p)
            start = x + 1
            va = _twice[C](va, p)
        elif doubled:
            va = _twice[C](va, p)
        for y in range(start, b1):
            var slot = ka + right.lows.unsafe_get(y)
            comptime if C == ZZ:
                acc.add_signed[wide, track](slot, va, right.values.unsafe_get(y))
            else:
                acc.add_modular[track](slot, UInt64(va), UInt64(right.values.unsafe_get(y)), p)


struct _Sink[C: CoefficientRing, packed: Bool](Movable):
    """Where products emit their terms: canonical sparse terms, or packed keys and coefficients in the product's radix."""
    var terms: List[Term[Self.C]]
    var keys: List[UInt64]
    var values: List[Self.C]

    def __init__(out self):
        self.terms = List[Term[Self.C]]()
        self.keys = List[UInt64]()
        self.values = List[Self.C]()

    @always_inline
    def add(mut self, packing: _Packing, key: UInt64, var coefficient: Self.C, mut budget: Budget) raises PolyError:
        comptime if Self.packed:
            self.keys.append(key)
            self.values.append(coefficient^)
        else:
            self.terms.append(Term[Self.C](packing.monomial(key, budget), coefficient^))


def _emit[C: CoefficientRing, packed: Bool, narrow: Bool](mut acc: _Accumulator[narrow], packing: _Packing, prefix: UInt64, template: C, mut sink: _Sink[C, packed], mut budget: Budget) raises PolyError:
    """Append the window's touched nonzero slots in increasing order and clear them."""
    var touched = acc.sorted_touched()
    _charge_size(len(touched), budget)
    _allocation[Term[C]](len(touched), budget)
    for slot in touched:
        if not _slot_zero(acc, slot, template):
            var key = prefix * packing.window + UInt64(slot)
            sink.add(packing, key, _coefficient(acc, slot, template), budget)
        acc.clear(slot)


def _emit_window[C: CoefficientRing, packed: Bool, narrow: Bool](mut acc: _Accumulator[narrow], packing: _Packing, prefix: UInt64, template: C, mut sink: _Sink[C, packed], mut budget: Budget) raises PolyError:
    """Append every nonzero slot of the window in increasing order and clear them."""
    var window = Int(packing.window)
    _charge_size(window, budget)
    for slot in range(window):
        if acc.empty(slot):
            continue
        if not _slot_zero(acc, slot, template):
            _allocation[Term[C]](1, budget)
            var key = prefix * packing.window + UInt64(slot)
            sink.add(packing, key, _coefficient(acc, slot, template), budget)
        acc.clear(slot)


# --------------------------------------------------------------- division

comptime _EXACT = 0
comptime _INEXACT = 1
"""The divisor does not divide the dividend."""
comptime _OVERFLOW = 2
"""A quotient coefficient needs more than 63 bits: the general path decides."""


def _divide_word(words: Tuple[UInt64, UInt64, UInt64], divisor: Int64, mut quotient: Int64) -> Int:
    """Divide a 192-bit two's complement value by a nonzero word: _EXACT, _INEXACT, or _OVERFLOW."""
    var sign = UInt64(0) - (words[0] >> 63)
    if words[1] == sign and words[2] == sign and not (divisor == -1 and words[0] == UInt64(1) << 63):
        # A value of one word divides natively; INT64_MIN / -1 alone would not fit.
        var v = Int64(words[0])
        var q = v // divisor
        if q * divisor != v:
            return _INEXACT
        quotient = q
        return _EXACT
    var value = _integer(words)
    try:
        var d = Integer(divisor)
        var q = value // d
        if q * d != value:
            return _INEXACT
        if q.magnitude_bit_length() > 63:
            return _OVERFLOW
        quotient = q.to_native_exact[DType.int64]()
        return _EXACT
    except:
        return _OVERFLOW


def _packed_divide[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], mut budget: Budget) raises PolyError -> Optional[Optional[SparsePoly[C]]]:
    """The exact quotient f / g (inner None when g does not divide f), or None when the kernel does not apply.

    Both operands are nonzero, in the same ring and order.
    """
    comptime if not (C == ZZ or C == ModInt):
        return None
    if not len(f.order.variables):
        return None
    if not _lex_identity(f) or not _lex_identity(g) or f.order != g.order:
        return None
    for term in g.terms:
        if not _small(term.coefficient):
            return None
    for term in f.terms:
        comptime if C == ZZ:
            if rebind[ZZ](term.coefficient).value.magnitude_bit_length() > 189:
                return None
    var n = len(f.order.variables)
    var top_f = _maxima(f, n)
    var top_g = _maxima(g, n)
    for i in range(n):
        if top_g[i] > top_f[i]:
            return Optional(Optional[SparsePoly[C]](None))
    # Quotient and divisor exponents both stay within the dividend's.
    var zeros = List[UInt64](length=n, fill=0)
    var radix = _radices(top_f, zeros)
    if not radix:
        return None
    var packing = _Packing(radix.take())
    var f_keys = List[UInt64](capacity=len(f.terms))
    for term in f.terms:
        budget.charge(1)
        f_keys.append(packing.key(term.monomial))
    var g_keys = List[UInt64](capacity=len(g.terms))
    var g_values = List[C](capacity=len(g.terms))
    for term in g.terms:
        budget.charge(1)
        g_keys.append(packing.key(term.monomial))
        g_values.append(term.coefficient.copy())
    var run = _division_core[C, True](f_keys, f.terms, List[C](), g_keys, g_values, top_f, top_g, packing.copy(), budget)
    if not run:
        return None
    if not run.value():
        return Optional(Optional[SparsePoly[C]](None))
    ref q = run.value().value()
    var template = g.terms[0].coefficient.copy()
    var result = SparsePoly[C].zero(f.order, f.ring, budget)
    _allocation[Term[C]](len(q[0]), budget)
    result.terms.reserve(len(q[0]))
    for x in range(len(q[0])):
        result.terms.append(Term[C](packing.monomial(q[0][x], budget), _word_coefficient(q[1][x], template)))
    return Optional(Optional(result^))


def _packed_poly_divide[C: CoefficientRing](f: _PackedPoly[C], g: _PackedPoly[C], mut budget: Budget) raises PolyError -> Optional[Optional[_PackedPoly[C]]]:
    """The exact quotient of packed operands, packed in its own radix (inner None when g does not divide f), or None when the kernel does not apply."""
    comptime if not (C == ZZ or C == ModInt):
        return None
    if not len(f) or not len(g) or f.order != g.order or not len(f.radix):
        return None
    for v in g.values:
        if not _small(v):
            return None
    comptime if C == ZZ:
        for v in f.values:
            if rebind[ZZ](v).value.magnitude_bit_length() > 189:
                return None
    var n = len(f.radix)
    var top_f = List[UInt64](capacity=n)
    var top_g = List[UInt64](capacity=n)
    for i in range(n):
        top_f.append(f.radix[i] - 1)
        top_g.append(g.radix[i] - 1)
        if top_g[i] > top_f[i]:
            return Optional(Optional[_PackedPoly[C]](None))
    var packing = _Packing(f.radix.copy())
    var run = _division_core[C, False](f.keys, List[Term[C]](), f.values, g.repacked(packing.radix), g.values, top_f, top_g, packing.copy(), budget)
    if not run:
        return None
    if not run.value():
        return Optional(Optional[_PackedPoly[C]](None))
    ref q = run.value().value()
    # The quotient's exponents are exactly top_f - top_g.
    var radix = List[UInt64](capacity=n)
    for i in range(n):
        radix.append(top_f[i] - top_g[i] + 1)
    var weights = List[UInt64](length=n, fill=1)
    var i = n - 2
    while i >= 0:
        weights[i] = weights[i + 1] * radix[i + 1]
        i -= 1
    var template = g.values[0].copy()
    _allocation[UInt64](len(q[0]), budget)
    _allocation[C](len(q[0]), budget)
    var keys = List[UInt64](capacity=len(q[0]))
    var values = List[C](capacity=len(q[0]))
    for x in range(len(q[0])):
        var rest = q[0][x]
        var key = UInt64(0)
        for v in range(n):
            var e = rest // packing.weight[v]
            rest -= e * packing.weight[v]
            key += e * weights[v]
        keys.append(key)
        values.append(_word_coefficient(q[1][x], template))
    return Optional(Optional(_PackedPoly[C](f.order.copy(), f.ring.copy(), radix^, keys^, values^)))


@always_inline
def _add_dividend[C: CoefficientRing](mut acc: _Accumulator, slot: Int, coefficient: C, wide: Bool, p: UInt64) raises PolyError:
    comptime if C == ZZ:
        if wide:
            acc.add_integer[True](slot, rebind[ZZ](coefficient).value)
        else:
            acc.add_integer[False](slot, rebind[ZZ](coefficient).value)
    else:
        acc.add_modular[True](slot, rebind[ModInt](coefficient).value, 1, p)


def _word_coefficient[C: CoefficientRing](value: Int64, template: C) -> C:
    comptime if C == ZZ:
        return rebind_var[C](ZZ(Integer(value)))
    else:
        return rebind_var[C](rebind[ModInt](template)._with_residue(UInt64(value)))


comptime _SCAN_DENSITY = 8
"""Window slots per dividend term in a chunk up to which division scans the
window rather than track every product's slot."""


@always_inline
def _subtract_products[C: CoefficientRing, track: Bool, wide: Bool](mut acc: _Accumulator, lows: List[Int], values: List[Int64], divisor: _Packed, g0: Int, g1: Int, p: UInt64):
    """Subtract a quotient chunk times divisor terms g0 ..< g1 from the window."""
    for x in range(len(lows)):
        var kq = lows.unsafe_get(x)
        var vq = values.unsafe_get(x)
        for y in range(g0, g1):
            var slot = kq + divisor.lows.unsafe_get(y)
            comptime if C == ZZ:
                acc.add_signed[wide, track](slot, -vq, divisor.values.unsafe_get(y))
            else:
                acc.add_modular[track](slot, UInt64(vq), p - UInt64(divisor.values.unsafe_get(y)), p)


def _division_core[C: CoefficientRing, terms: Bool](f_keys: List[UInt64], f_terms: List[Term[C]], f_values: List[C], g_keys: List[UInt64], g_values: List[C], top_f: List[UInt64], top_g: List[UInt64], var packing: _Packing, mut budget: Budget) raises PolyError -> Optional[Optional[Tuple[List[UInt64], List[Int64]]]]:
    """Chunked exact division on increasing keys in the dividend's radix: the quotient's keys (in that radix) and word coefficients, inner None when inexact, None when the kernel does not apply.

    Divisor coefficients fit a word and ZZ dividend coefficients 189 bits; the front ends check.
    The dividend's coefficients are f_terms' when terms, else f_values (the other list is unused).
    """
    var n = len(top_f)
    # Deepen the prefix until the divisor's leading chunk is its leading term
    # alone; with every variable in the prefix each chunk is a single term.
    var lead_key = g_keys[len(g_keys) - 1]
    var depth = max(packing.prefix_count, 1)
    while True:
        packing.split(depth)
        var lead_prefix = lead_key // packing.window
        var count = 0
        var i = len(g_keys) - 1
        while i >= 0 and count < 2:
            budget.charge(1)
            if g_keys[i] // packing.window != lead_prefix:
                break
            count += 1
            i -= 1
        if count == 1 or depth >= n:
            break
        depth += 1
    var divisor = _pack_keys(g_keys.copy(), g_values, packing, budget)
    var top = len(divisor.group_prefix) - 1
    var lead_value = divisor.values[len(divisor.values) - 1]
    var h = divisor.group_prefix[top]
    var template = g_values[0].copy()
    var p = UInt64(0)
    var lead_inverse = UInt64(0)
    comptime if C == ModInt:
        ref t = rebind[ModInt](template)
        p = t.modulus
        lead_inverse = rebind[ModInt](g_values[len(g_values) - 1]).inverse(budget).value
    # The dividend's chunks, by prefix.
    var dividend_prefix = List[UInt64]()
    var dividend_start = List[Int]()
    for i in range(len(f_keys)):
        budget.charge(1)
        var prefix = f_keys[i] // packing.window
        if not len(dividend_prefix) or dividend_prefix[len(dividend_prefix) - 1] != prefix:
            dividend_prefix.append(prefix)
            dividend_start.append(i)
    dividend_start.append(len(f_keys))
    # Every dividend prefix visits every divisor group. When that exceeds the
    # dividend's terms, chunks collide too little to pay; a heap division wins.
    if _product(len(dividend_prefix), len(divisor.group_prefix)) > 8 * len(f_keys):
        return None
    var window = Int(packing.window)
    _allocation[UInt64](_product(4, window), budget)
    # Each remainder slot sums one dividend coefficient and at most one quotient
    # product (of at most 63 bits) per divisor term.
    var wide = True
    comptime if C == ZZ:
        var dividend_bits = 0
        comptime if terms:
            for term in f_terms:
                dividend_bits = max(dividend_bits, rebind[ZZ](term.coefficient).value.magnitude_bit_length())
        else:
            for v in f_values:
                dividend_bits = max(dividend_bits, rebind[ZZ](v).value.magnitude_bit_length())
        var divisor_bits = 0
        for v in g_values:
            divisor_bits = max(divisor_bits, rebind[ZZ](v).value.magnitude_bit_length())
        wide = max(dividend_bits, 63 + divisor_bits + _bit_length(len(g_values))) + 1 >= 127
    var acc = _Accumulator(window, wide)
    # A window the dividend's chunks fill densely is scanned for its nonzero
    # slots, already in order, instead of tracking every product's slot.
    var scan = window <= _SCAN_DENSITY * max(1, len(f_keys) // len(dividend_prefix))
    # Quotient chunks in the order computed (decreasing prefix), as packed
    # terms, with each term's slot in the window.
    var q_prefix = List[UInt64]()
    var q_keys = List[List[UInt64]]()
    var q_lows = List[List[Int]]()
    var q_values = List[List[Int64]]()
    var q_index = Dict[UInt64, Int]()
    # Pending prefixes, processed largest first; duplicates are skipped.
    var pending = List[UInt64]()
    for prefix in dividend_prefix:
        pending.append(prefix)
    var next_dividend = len(dividend_prefix) - 1
    var last = UInt64.MAX
    var heap = _MaxHeap(pending^)
    while heap.size():
        var t = heap.pop()
        if t == last:
            continue
        last = t
        budget.charge(1)
        # The dividend's chunk at t, if any.
        while next_dividend >= 0 and dividend_prefix[next_dividend] > t:
            next_dividend -= 1
        if next_dividend >= 0 and dividend_prefix[next_dividend] == t:
            for i in range(dividend_start[next_dividend], dividend_start[next_dividend + 1]):
                var slot = Int(f_keys[i] % packing.window)
                comptime if terms:
                    _add_dividend(acc, slot, f_terms[i].coefficient, wide, p)
                else:
                    _add_dividend(acc, slot, f_values[i], wide, p)
        # Minus the known products: quotient chunks at t - e for divisor prefixes e below h.
        for e_index in range(top):
            var e = divisor.group_prefix[e_index]
            if e > t:
                break
            var d = t - e
            var found = q_index.get(d)
            if not found:
                continue
            var qi = found.value()
            ref ql = q_lows[qi]
            ref qv = q_values[qi]
            var g0 = divisor.group_start[e_index]
            var g1 = divisor.group_start[e_index + 1]
            _charge_size(_product(len(ql), g1 - g0), budget)
            if scan:
                if wide:
                    _subtract_products[C, False, True](acc, ql, qv, divisor, g0, g1, p)
                else:
                    _subtract_products[C, False, False](acc, ql, qv, divisor, g0, g1, p)
            elif wide:
                _subtract_products[C, True, True](acc, ql, qv, divisor, g0, g1, p)
            else:
                _subtract_products[C, True, False](acc, ql, qv, divisor, g0, g1, p)
        # The remainder chunk must be the next quotient chunk times the lead term.
        var touched: List[Int]
        if scan:
            touched = List[Int]()
            for slot in range(window):
                if acc.touched_flag.unsafe_get(slot) or not acc.empty(slot):
                    touched.append(slot)
            acc.touched.clear()
        else:
            touched = acc.sorted_touched()
        _charge_size(len(touched), budget)
        var keys = List[UInt64]()
        var values = List[Int64]()
        var status = _EXACT
        for slot in touched:
            if status != _EXACT or _slot_zero(acc, slot, template):
                acc.clear(slot)
                continue
            var full = t * packing.window + UInt64(slot)
            if t < h:
                status = _INEXACT
            else:
                # Exponentwise divisibility, and room for the divisor's exponents.
                for i in range(n):
                    var e_full = packing.exponent(full, i)
                    var e_lead = packing.exponent(lead_key, i)
                    if e_full < e_lead or e_full - e_lead + top_g[i] > top_f[i]:
                        status = _INEXACT
                        break
                if status == _EXACT:
                    comptime if C == ZZ:
                        var value = Int64(0)
                        status = _divide_word(acc.signed_words(slot), lead_value, value)
                        if status == _EXACT:
                            keys.append(full - lead_key)
                            values.append(value)
                    else:
                        keys.append(full - lead_key)
                        values.append(Int64(_mod_mul(acc.modular(slot, p), lead_inverse, p)))
            acc.clear(slot)
        if status == _OVERFLOW:
            return None
        if status == _INEXACT:
            return Optional(Optional[Tuple[List[UInt64], List[Int64]]](None))
        if len(keys):
            var d = t - h
            q_index[d] = len(q_prefix)
            q_prefix.append(d)
            for e_index in range(top):
                var e = divisor.group_prefix[e_index]
                heap.push(d + e)
            var lows = List[Int](capacity=len(keys))
            for key in keys:
                lows.append(Int(key % packing.window))
            q_keys.append(keys^)
            q_lows.append(lows^)
            q_values.append(values^)
    # Chunks were found from the top down; emit them in increasing order.
    var out_keys = List[UInt64]()
    var out_values = List[Int64]()
    var c = len(q_prefix) - 1
    while c >= 0:
        ref keys = q_keys[c]
        ref values = q_values[c]
        for x in range(len(keys)):
            out_keys.append(keys[x])
            out_values.append(values[x])
        c -= 1
    return Optional(Optional((out_keys^, out_values^)))


struct _MaxHeap(Movable):
    var items: List[UInt64]

    def __init__(out self, var items: List[UInt64]):
        self.items = items^
        var i = len(self.items) // 2 - 1
        while i >= 0:
            self._down(i)
            i -= 1

    def size(self) -> Int:
        return len(self.items)

    def push(mut self, value: UInt64):
        self.items.append(value)
        var i = len(self.items) - 1
        while i > 0:
            var parent = (i - 1) // 2
            if self.items[parent] >= self.items[i]:
                break
            var held = self.items[parent]
            self.items[parent] = self.items[i]
            self.items[i] = held
            i = parent

    def pop(mut self) -> UInt64:
        var top = self.items[0]
        var last = self.items.pop()
        if len(self.items):
            self.items[0] = last
            self._down(0)
        return top

    def _down(mut self, start: Int):
        var i = start
        var n = len(self.items)
        while True:
            var largest = i
            var l = 2 * i + 1
            var r = l + 1
            if l < n and self.items[l] > self.items[largest]:
                largest = l
            if r < n and self.items[r] > self.items[largest]:
                largest = r
            if largest == i:
                return
            var held = self.items[largest]
            self.items[largest] = self.items[i]
            self.items[i] = held
            i = largest


# ----------------------------------------------------------------- powers

def _exponent(m: Monomial, v: Int) -> Int:
    for pair in m.pairs:
        if Int(pair[0]) == v:
            return Int(pair[1])
    return 0


def _grading[C: CoefficientRing](f: SparsePoly[C], n: Int) -> Int:
    """A variable whose exponent is strictly largest at f's leading term, or -1."""
    ref lead = f.terms[len(f.terms) - 1].monomial
    for v in range(n):
        var top = _exponent(lead, v)
        var unique = top > 0
        for i in range(len(f.terms) - 1):
            if _exponent(f.terms[i].monomial, v) >= top:
                unique = False
                break
        if unique:
            return v
    return -1


def _term_estimate(t: Int, j: Int, top: List[UInt64]) -> Float64:
    """Terms of f^j for t terms, at most C(t + j - 1, j) and the exponent box (j top_i + 1)."""
    var choose = Float64(1)
    for i in range(1, t):
        choose = choose * Float64(j + i) / Float64(i)
    var box = Float64(1)
    for e in top:
        box *= Float64(j) * Float64(e) + 1
    return min(choose, box)


def _packed_power[C: CoefficientRing](f: SparsePoly[C], k: Int, mut budget: Budget) raises PolyError -> Optional[SparsePoly[C]]:
    """f^k by the power recurrence, or None when repeated squaring is expected to be faster or the kernel does not apply.

    The recurrence takes about (terms of f - 1) * (terms of f^k) heap steps and
    the last square of repeated squaring about (terms of f^(k/2))^2 / 2 products;
    both counts are estimated from the exponents alone.

    With D the derivation multiplying a monomial by its exponent of a variable v,
    g = f^k satisfies f D(g) = k D(f) g. If f's leading term f_0 has the
    strictly largest exponent of v, the coefficient of f_0 g_j in that identity
    gives g_j from the larger terms of g:
        g_j (d_j - k e_0) f_0 = - sum over i > 0, f_i g_l = f_0 g_j of f_i g_l (d_l - k e_i),
    where e_i and d_l are exponents of v. A heap merges the products f_i g as
    g grows, from the top down (Monagan and Pearce), so the work is about
    (terms of f - 1) * (terms of g) instead of the products of the last square.
    """
    if not _prefers_recurrence(f, k):
        return None
    return _power_recurrence(f, k, budget)


def _packed_poly_power[C: CoefficientRing](f: _PackedPoly[C], k: Int, mut budget: Budget) raises PolyError -> Optional[_PackedPoly[C]]:
    """f^k for k >= 1, packed, or None when f does not fit the packed kernels.

    Takes the power recurrence where SparsePoly.pow would, and otherwise
    squares from the exponent's leading bit down with packed products; a step
    whose operands outgrow the word kernel multiplies as SparsePoly does.
    """
    comptime if not (C == ZZ or C == ModInt):
        return None
    if not len(f) or k < 1:
        return None
    var base = f.to_sparse(budget)
    comptime if C == ZZ:
        if _prefers_recurrence(base, k):
            var run = _packed_recurrence(base, k, budget)
            if run:
                return run.take()
    var result = f.copy()
    var bit = 63
    while not (UInt64(k) >> UInt64(bit)) & 1:
        bit -= 1
    while bit > 0:
        budget.charge(1)
        bit -= 1
        var squared = _packed_step(result, result, budget)
        if not squared:
            return None
        result = squared.take()
        if (UInt64(k) >> UInt64(bit)) & 1:
            var times = _packed_step(result, f, budget)
            if not times:
                return None
            result = times.take()
    return result^


def _packed_step[C: CoefficientRing](a: _PackedPoly[C], b: _PackedPoly[C], mut budget: Budget) raises PolyError -> Optional[_PackedPoly[C]]:
    """a * b packed: the packed kernel, or SparsePoly.mul when the operands outgrow it."""
    var direct = _packed_poly_product(a, b, budget)
    if direct:
        return direct^
    var product = a.to_sparse(budget).mul(b.to_sparse(budget), budget)
    return _PackedPoly[C].from_sparse(product, budget)


def _prefers_recurrence[C: CoefficientRing](f: SparsePoly[C], k: Int) -> Bool:
    """Whether the power recurrence is expected to beat repeated squaring for f^k."""
    comptime if not (C == ZZ):
        return False
    var t = len(f.terms)
    if t < 2 or k < 3 or not _usable(f, f):
        return False
    var n = len(f.order.variables)
    var top = _maxima(f, n)
    # Bases with a wide exponent spread (Monagan and Pearce's u^5) give sparse
    # windows that repeated squaring fills cheaply; the recurrence measured
    # faster only for exponents of at most _POWER_SPREAD.
    for e in top:
        if e > _POWER_SPREAD:
            return False
    return Float64(t - 1) * _term_estimate(t, k, top) < Float64(_POWER_GATE) * _term_estimate(t, k // 2, top) ** 2 / 2


@fieldwise_init
struct _Recurrence(Movable):
    """f^k from the recurrence, from its trailing term up: exponent vectors and coefficients."""
    var n: Int
    var exponents: List[Int]
    var coefficients: List[Int128]


def _power_recurrence[C: CoefficientRing](f: SparsePoly[C], k: Int, mut budget: Budget) raises PolyError -> Optional[SparsePoly[C]]:
    """f^k by the power recurrence for any k >= 2, or None when the kernel does not apply."""
    var run = _recurrence(f, k, budget)
    if not run:
        return None
    var g = run.take()
    var n = g.n
    ref ge = g.exponents
    ref gc = g.coefficients
    var result = SparsePoly[C].zero(f.order, f.ring, budget)
    var count = len(gc)
    _allocation[Term[C]](count, budget)
    result.terms.reserve(count)
    var j = count - 1
    while j >= 0:
        budget.charge(1)
        var used = 0
        for i in range(n):
            used += Int(ge[j * n + i] != 0)
        _allocation[Tuple[Int64, UInt32]](used, budget)
        var monomial = Monomial()
        monomial.pairs.reserve(used)
        for i in range(n):
            if ge[j * n + i]:
                monomial.pairs.append((Int64(i), UInt32(ge[j * n + i])))
        var c = gc[j]
        var words = UInt128(c)
        result.terms.append(Term[C](monomial^, rebind_var[C](ZZ(_integer((UInt64(words), UInt64(words >> 64), UInt64(0) - UInt64(Int(c < 0))))))))
        j -= 1
    return result^


def _packed_recurrence[C: CoefficientRing](f: SparsePoly[C], k: Int, mut budget: Budget) raises PolyError -> Optional[_PackedPoly[C]]:
    """f^k by the power recurrence, packed in the radix k * top + 1, or None when the kernel does not apply."""
    var run = _recurrence(f, k, budget)
    if not run:
        return None
    ref g = run.value()
    var n = g.n
    var top = _maxima(f, n)
    var radix = List[UInt64](capacity=n)
    for e in top:
        radix.append(UInt64(k) * e + 1)
    var weights = List[UInt64](length=n, fill=1)
    var i = n - 2
    while i >= 0:
        weights[i] = weights[i + 1] * radix[i + 1]
        i -= 1
    var count = len(g.coefficients)
    _allocation[UInt64](count, budget)
    _allocation[C](count, budget)
    var keys = List[UInt64](capacity=count)
    var values = List[C](capacity=count)
    var j = count - 1
    while j >= 0:
        budget.charge(1)
        var key = UInt64(0)
        for v in range(n):
            key += UInt64(g.exponents[j * n + v]) * weights[v]
        keys.append(key)
        var c = g.coefficients[j]
        var words = UInt128(c)
        values.append(rebind_var[C](ZZ(_integer((UInt64(words), UInt64(words >> 64), UInt64(0) - UInt64(Int(c < 0)))))))
        j -= 1
    return _PackedPoly[C](f.order.copy(), f.ring.copy(), radix^, keys^, values^)


def _recurrence[C: CoefficientRing](f: SparsePoly[C], k: Int, mut budget: Budget) raises PolyError -> Optional[_Recurrence]:
    """The power recurrence's exponent vectors and coefficients for f^k, or None when the kernel does not apply."""
    comptime if not (C == ZZ):
        return None
    var t = len(f.terms)
    if t < 2 or k < 2 or not _usable(f, f):
        return None
    var n = len(f.order.variables)
    var v = _grading(f, n)
    if v < 0:
        return None
    var top = _maxima(f, n)
    # Products f_i g_l have exponents up to (k + 1) top_i: keys carry no digit.
    var scaled = List[UInt64](capacity=n)
    var zeros = List[UInt64](length=n, fill=0)
    for i in range(n):
        scaled.append(UInt64(k + 1) * top[i])
    var radix = _radices(scaled, zeros)
    if not radix:
        return None
    var packing = _Packing(radix.take())
    # |g| <= |f|_1^k; each recurrence sum adds t - 1 products of at most
    # max|f_i| * (k + 1) top_v * |g| and must stay below 2^127.
    var norm = Float64(0)
    var largest = Float64(1)
    for term in f.terms:
        var c = abs(Float64(_signed(term.coefficient)))
        norm += c
        largest = max(largest, c)
    var bits = Float64(k) * log2(norm) + log2(Float64(t) * largest * Float64(k + 1) * (Float64(top[v]) + 1)) + 2
    if bits > 124:
        return None
    var key_bits = _bit_length(Int(packing.weight[0] * packing.radix[0]))
    var stream_bits = _bit_length(t - 1)
    if key_bits + stream_bits > 63:
        return None
    # f's terms from the leading one down.
    var fk = List[UInt64](capacity=t)
    var fc = List[Int128](capacity=t)
    var fd = List[Int](capacity=t)
    var fe = List[Int](capacity=t * n)
    for x in range(t):
        ref term = f.terms[t - 1 - x]
        budget.charge(1)
        fk.append(packing.key(term.monomial))
        fc.append(Int128(_signed(term.coefficient)))
        fd.append(_exponent(term.monomial, v))
        for i in range(n):
            fe.append(_exponent(term.monomial, i))
    _allocation[UInt64](_product(n + 3, t + 1), budget)
    var expected = Int(min(_term_estimate(t, k, top), Float64(1 << 20)))
    var gk = List[UInt64](capacity=expected)
    var gc = List[Int128](capacity=expected)
    var gd = List[Int](capacity=expected)
    var ge = List[Int](capacity=expected * n)
    var lead = Int128(1)
    for _ in range(k):
        lead *= fc[0]
    gk.append(UInt64(k) * fk[0])
    gc.append(lead)
    gd.append(k * fd[0])
    for i in range(n):
        ge.append(k * fe[i])
    var shift = UInt64(stream_bits)
    var mask = (UInt64(1) << shift) - 1
    var entries = List[UInt64](capacity=t - 1)
    for x in range(1, t):
        entries.append(((fk[x] + gk[0]) << shift) | UInt64(x))
    var heap = _StreamHeap(entries^)
    var position = List[Int](length=t, fill=0)
    var stalled = List[Int](capacity=t)
    var exponents = List[Int](length=n, fill=0)
    while heap.size:
        var key = heap.top() >> shift
        var sum = Int128(0)
        var x = 0
        var l = 0
        while heap.size and (heap.top() >> shift) == key:
            budget.charge(1)
            x = Int(heap.top() & mask)
            l = position[x]
            sum += fc[x] * Int128(gd[l] - k * fd[x]) * gc[l]
            position[x] = l + 1
            if l + 1 < len(gk):
                heap.replace_top(((fk[x] + gk[l + 1]) << shift) | UInt64(x))
            else:
                heap.pop()
                stalled.append(x)
        if sum == 0:
            continue
        # The new term is (f_x g_l) / f_0, which the identity makes exact.
        for i in range(n):
            exponents[i] = fe[x * n + i] + ge[l * n + i] - fe[i]
            if exponents[i] < 0:
                return None
        var d = fd[x] + gd[l] - fd[0]
        var denominator = fc[0] * Int128(d - k * fd[0])
        var quotient: Int128
        # Most sums fit a word, where division is a machine instruction.
        if Int128(Int64.MIN) < sum and sum <= Int128(Int64.MAX) and Int128(Int64.MIN) < denominator and denominator <= Int128(Int64.MAX):
            var word = Int64(sum)
            var by = Int64(denominator)
            if word % by != 0:
                return None
            quotient = Int128(word // by)
        else:
            quotient = sum // denominator
            if quotient * denominator != sum:
                return None
        _allocation[UInt64](n + 3, budget)
        gk.append(key - fk[0])
        gc.append(-quotient)
        gd.append(d)
        for i in range(n):
            ge.append(exponents[i])
        for s in stalled:
            heap.push(((fk[s] + gk[position[s]]) << shift) | UInt64(s))
        stalled.clear()
    return _Recurrence(n, ge^, gc^)


struct _StreamHeap(Movable):
    """A max-heap of at most its initial number of words, without reallocation."""
    var items: List[UInt64]
    var size: Int

    def __init__(out self, var items: List[UInt64]):
        self.size = len(items)
        self.items = items^
        var i = self.size // 2 - 1
        while i >= 0:
            self._down(i, self.items.unsafe_get(i))
            i -= 1

    @always_inline
    def top(self) -> UInt64:
        return self.items.unsafe_get(0)

    @always_inline
    def replace_top(mut self, value: UInt64):
        self._down(0, value)

    @always_inline
    def pop(mut self):
        self.size -= 1
        if self.size:
            self._down(0, self.items.unsafe_get(self.size))

    def push(mut self, value: UInt64):
        var i = self.size
        self.size += 1
        while i > 0:
            var parent = (i - 1) // 2
            var above = self.items.unsafe_get(parent)
            if above >= value:
                break
            self.items.unsafe_get(i) = above
            i = parent
        self.items.unsafe_get(i) = value

    @always_inline
    def _down(mut self, start: Int, value: UInt64):
        """Place value at start and sift it down."""
        var i = start
        while True:
            var child = 2 * i + 1
            if child >= self.size:
                break
            var c = self.items.unsafe_get(child)
            if child + 1 < self.size:
                var other = self.items.unsafe_get(child + 1)
                if other > c:
                    child += 1
                    c = other
            if c <= value:
                break
            self.items.unsafe_get(i) = c
            i = child
        self.items.unsafe_get(i) = value
