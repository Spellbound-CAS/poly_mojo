"""Number-theoretic transforms in double precision, for products that need few primes.

Residues modulo primes 2^49.99 < p < 2^50 are exact integers in a Float64, and
every operation below is exact. The bounds hold for any residues below p and
any length, so no input needs the integer kernel instead. The argument needs
IEEE double arithmetic with round-to-nearest and a correctly rounded fma
(hardware FMA, or a correctly rounded software fma; both qualify). With
u = 2^-53:

- _mulmod(a, b) for |a b| <= M p^2, M <= 2. h = fl(a b) has |a b - h| <= u M p^2
  < 2^48, and l = fma(a, b, -h) is that error exactly. fl(h fl(1/p)) is within
  3u M p < 0.375 M of a b / p, so the rounded quotient Q is within one of the
  nearest integer to a b / p and leaves |a b - Q p| <= (0.5 + 0.375 M) p. Q is
  not corrected: the result is a b - Q p, a correct residue in that wider range
  (within 1.25 p), and the stage invariants below absorb it. fma(-Q, p, h) =
  h - Q p is an integer below 2^52, so it is exact, and adding l is exact.
- _reduce(x) for |x| < 4 p: x - Q p, within p / 2 + 1, exactly by one fma.
- _float_digits, for any K of the eight primes (largest first, so q_j - q_i
  < 0.007 q_i for j < i), computes digit i modulo p = q_i from residues and
  earlier digits in [0, 1.007 p). For i = 1, r_1 - v_0 lies within 1.007 p
  and its product with the prefix inverse (below p) has M <= 1.007. For
  i >= 2, each Horner step multiplies acc (within 1.007 p) by q_j mod q_i =
  q_j - q_i, so M < 0.01 and the product is within 0.504 p; adding v_j keeps
  it below 1.52 p, and _reduce brings it within p / 2 + 1. Then r_i - acc
  lies within 1.5 p + 1, and its product with the prefix inverse has
  M <= 1.51, within 1.07 p. _canonical reduces that below 4 p to within
  p / 2 + 1 and adds p to a negative result, giving the digit in [0, q_i).
- Twiddles are reduced, within p / 2 + 1.
- Forward (decimation in frequency): inputs lie in [0, p). With values within p,
  u + v is reduced and (u - v) w has M <= 1, so outputs stay within 0.875 p.
- Radix-3 layer of 3 m points, forward: from inputs in [0, p), t = omega (x1 - x2)
  has M <= 0.5 (within 0.69 p); x0 + x1 + x2 < 3 p and x0 - x2 + t, x0 - x1 - t
  (below 1.7 p) are reduced to within p / 2 + 1, and the twiddle products then
  have M <= 0.26: within 0.6 p, inside the power-of-two stages' input range.
  Inverse, from values within p / 2 + 1: the twiddle products are within
  0.6 p, t = omega (y2 - y1) has M <= 0.6 (within 0.73 p), and the three
  sums, below 1.9 p, are reduced to within p / 2 + 1.
- Pointwise products of values within p have M <= 1: within 0.875 p.
- Inverse (decimation in time): t = v w has M <= 0.5, within 0.69 p, and both
  u + t and u - t (below 1.6 p) are reduced to within p / 2 + 1.
- The final scaling by 1/n < p has M <= 0.5. Every sum is of integers below
  2^52, hence exact.

Eight-lane vectors run every butterfly, the last three stages through lane
shuffles: a convolution runs 2.4-2.7 times faster than the scalar 61-bit
kernel, about twice as fast per bit of modulus. The tables follow _ntt's
layout, entry m + j holding w_2m^j or w_2m^-j.
"""
from std.math import fma

comptime _FPRIME_COUNT = 8
comptime _FPRIMES: Array[UInt64, _FPRIME_COUNT] = [
    1125844072267777, 1125818302464001, 1125625028935681, 1125122517762049,
    1124903474429953, 1124877704626177, 1124658661294081, 1124581351882753,
]
"""Primes c * 2^32 + 1 just below 2^50 with 3 | c, largest first: lengths 2^k and 3 * 2^k both divide p - 1."""
comptime _FROOTS: Array[UInt64, _FPRIME_COUNT] = [
    786008014450235, 147641925747491, 417876965932711, 616716876639881,
    601605979588170, 511210619017435, 288846415555182, 641149848381095,
]
"""Roots of unity of order 2^32: the cubes of _FCUBE_ROOTS."""
comptime _FCUBE_ROOTS: Array[UInt64, _FPRIME_COUNT] = [
    909181527941361, 235975324430315, 880865589704237, 283701703289718,
    819158334494372, 442533454869943, 65778916453337, 425046294173167,
]
"""Roots of unity of order 3 * 2^32."""
comptime _FPRIME_BITS = 49
"""Each prime exceeds 2^49.99: k primes hold a symmetric range of 49 k bits."""
comptime _FMAX_LOG = 32
comptime _W = 8


def _fpower(a: UInt64, n: UInt64, q: UInt64) -> UInt64:
    """a^n mod a float prime q, by exact float products."""
    var p = SIMD[DType.float64, 1](Float64(q))
    var inverse = 1.0 / p
    var base = SIMD[DType.float64, 1](Float64(a % q))
    var result = SIMD[DType.float64, 1](1.0)
    var exponent = n
    while exponent:
        if exponent & 1:
            result = _reduce[1](_mulmod[1](result, base, p, inverse), p, inverse)
        exponent >>= 1
        if exponent:
            base = _reduce[1](_mulmod[1](base, base, p, inverse), p, inverse)
    var r = result[0]
    return UInt64(r + Float64(q)) if r < 0 else UInt64(r)


@always_inline
def _mulmod[w: Int](a: SIMD[DType.float64, w], b: SIMD[DType.float64, w], p: SIMD[DType.float64, w], inverse: SIMD[DType.float64, w]) -> SIMD[DType.float64, w]:
    """a * b mod p within (0.5 + 0.375 M) p, for |a b| <= M p^2 with M <= 2."""
    var h = a * b
    var l = fma(a, b, -h)
    var q = (h * inverse).__round__()
    return fma(-q, p, h) + l


@always_inline
def _reduce[w: Int](x: SIMD[DType.float64, w], p: SIMD[DType.float64, w], inverse: SIMD[DType.float64, w]) -> SIMD[DType.float64, w]:
    """x mod p within p / 2 + 1, for |x| < 4 p."""
    return fma(-(x * inverse).__round__(), p, x)


@always_inline
def _fmul(x: UInt64, y: UInt64, q: UInt64) -> UInt64:
    """x y mod a float prime q in [0, q), for x, y < q: a product with M <= 1, then a reduction."""
    var p = SIMD[DType.float64, 1](Float64(q))
    var inverse = SIMD[DType.float64, 1](1.0 / Float64(q))
    var r = _reduce[1](_mulmod[1](SIMD[DType.float64, 1](Float64(x)), SIMD[DType.float64, 1](Float64(y)), p, inverse), p, inverse)[0]
    return UInt64(r + Float64(q)) if r < 0 else UInt64(r)


struct _FTables(Movable):
    """Twiddles as signed floats within p / 2: w^-j is -w^(n/2 - j).

    A length n = 3 m runs one radix-3 layer over three power-of-two
    transforms of length m: `radix` holds w^j, w^2j, w^-j and w^-2j for
    j < m, w a primitive n-th root, and `omega` the cube root w^m.
    """
    var forward: List[Float64]
    var inverse: List[Float64]
    var third: Int
    """m when the length is 3 m, else 0."""
    var radix: List[Float64]
    var omega: Float64

    def __init__(out self, q: UInt64, root: UInt64, length: Int, cube_root: UInt64 = 0):
        """Tables for `length` points; a length 3 * 2^k needs the root of order 3 * 2^32, a power of two the root of order 2^32."""
        var n = length
        self.third = 0
        self.radix = List[Float64]()
        self.omega = 0
        var log = 0
        while (1 << log) < length:
            log += 1
        # A primitive n-th root for the power-of-two tables.
        var w = UInt64(1)
        if length >= 3 and length % 3 == 0:
            n = length // 3
            self.third = n
            self.radix = List[Float64](length=4 * n, fill=0)
            # w_3m = R^(2^32 / m) for R of order 3 * 2^32; from w^(m - 1),
            # omega = w^m and w^-1 = w^(m - 1) omega^2 follow by products.
            var w3 = _fpower(cube_root, UInt64((1 << _FMAX_LOG) // n), q)
            var rp = self.radix.unsafe_ptr()
            _powers(rp, n, w3, q)
            var top = _canonical_word(rp.unsafe_offset(n - 1).unsafe_load(), q)
            var omega = _fmul(top, w3, q)
            var w_inverse = _fmul(top, _fmul(omega, omega, q), q)
            _powers(rp.unsafe_offset(n), n, _fmul(w3, w3, q), q)
            _powers(rp.unsafe_offset(2 * n), n, w_inverse, q)
            _powers(rp.unsafe_offset(3 * n), n, _fmul(w_inverse, w_inverse, q), q)
            self.omega = _signed(omega, q)
            w = _fmul(w3, _fmul(w3, w3, q), q)
        elif length >= 2:
            w = _fpower(root, UInt64(1) << UInt64(_FMAX_LOG - log), q)
        # Every entry from index 1 up is written below; index 0 is unused.
        self.forward = List[Float64](unsafe_uninit_length=max(n, 1))
        self.inverse = List[Float64](unsafe_uninit_length=max(n, 1))
        self.forward[0] = 0
        self.inverse[0] = 0
        if n < 2:
            return
        var half = n // 2
        var fw = self.forward.unsafe_ptr()
        # Signed representatives keep every twiddle within p / 2.
        _powers(fw.unsafe_offset(half), half, w, q)
        var iv = self.inverse.unsafe_ptr()
        iv.unsafe_offset(half).unsafe_store(1.0)
        for j in range(1, half):
            iv.unsafe_offset(half + j).unsafe_store(-fw.unsafe_offset(n - j).unsafe_load())
        # Level m holds w_2m^j = w_4m^2j: every other entry of level 2m, read
        # in order from the level above rather than strided from the top.
        var m = half // 2
        while m >= 1:
            var j = 0
            if m >= _W:
                while j < m:
                    var f = fw.unsafe_offset(2 * m + 2 * j).unsafe_load[width=_W]().join(fw.unsafe_offset(2 * m + 2 * j + _W).unsafe_load[width=_W]())
                    var g = iv.unsafe_offset(2 * m + 2 * j).unsafe_load[width=_W]().join(iv.unsafe_offset(2 * m + 2 * j + _W).unsafe_load[width=_W]())
                    fw.unsafe_offset(m + j).unsafe_store(f.deinterleave()[0])
                    iv.unsafe_offset(m + j).unsafe_store(g.deinterleave()[0])
                    j += _W
            while j < m:
                fw.unsafe_offset(m + j).unsafe_store(fw.unsafe_offset(2 * m + 2 * j).unsafe_load())
                iv.unsafe_offset(m + j).unsafe_store(iv.unsafe_offset(2 * m + 2 * j).unsafe_load())
                j += 1
            m //= 2




@always_inline
def _canonical_word(x: Float64, q: UInt64) -> UInt64:
    """The residue in [0, q) of an integer float within q."""
    return UInt64(x + Float64(q)) if x < 0 else UInt64(x)


def _inverse_length(n: Int, q: UInt64) -> UInt64:
    """1 / n mod a float prime q = c 2^32 + 1, for n = 2^k or 3 * 2^k: 2^-k = q - c 2^(32 - k) and 3^-1 = (2 q + 1) / 3."""
    var m = n // 3 if n % 3 == 0 else n
    var k = 0
    while (1 << k) < m:
        k += 1
    var inverse = q - (((q - 1) >> 32) << UInt64(_FMAX_LOG - k))
    if n % 3 == 0:
        inverse = _fmul(inverse, (2 * q + 1) // 3, q)
    return inverse


@always_inline
def _signed(x: UInt64, q: UInt64) -> Float64:
    """The representative of x < q within q / 2."""
    return Float64(x) - Float64(q) if x > q // 2 else Float64(x)


def _powers(dst: Pointer[mut=True, Float64, _], count: Int, w: UInt64, q: UInt64):
    """w^0 ..< w^count (count < 8 or a multiple of 8) as signed floats within p / 2 + 1: eight exact entries, one chain of w^8 to 32, then 32 chains of w^32, four vectors abreast."""
    comptime span = 4 * _W
    var x = UInt64(1)
    for j in range(min(count, _W)):
        dst.unsafe_offset(j).unsafe_store(_signed(x, q))
        x = _fmul(x, w, q)
    if count > _W:
        var prime = Float64(q)
        var p = SIMD[DType.float64, _W](prime)
        var inv = SIMD[DType.float64, _W](1.0 / prime)
        # Entries 8 ..< 32 by one chain of w^8, then four chains of w^32.
        var step8 = SIMD[DType.float64, _W](_signed(x, q))
        var k = _W
        while k < min(count, span):
            dst.unsafe_offset(k).unsafe_store(_reduce[_W](_mulmod[_W](dst.unsafe_offset(k - _W).unsafe_load[width=_W](), step8, p, inv), p, inv))
            k += _W
        var x16 = _fmul(x, x, q)
        x = _fmul(x16, x16, q)
        var step = SIMD[DType.float64, _W](_signed(x, q))
        var j = span
        while j + span <= count:
            comptime for v in range(4):
                var previous = dst.unsafe_offset(j + v * _W - span).unsafe_load[width=_W]()
                dst.unsafe_offset(j + v * _W).unsafe_store(_reduce[_W](_mulmod[_W](previous, step, p, inv), p, inv))
            j += span
        while j < count:
            var previous = dst.unsafe_offset(j - span).unsafe_load[width=1]()
            dst.unsafe_offset(j).unsafe_store(_reduce[1](_mulmod[1](previous, SIMD[DType.float64, 1](_signed(x, q)), SIMD[DType.float64, 1](prime), SIMD[DType.float64, 1](1.0 / prime)), SIMD[DType.float64, 1](prime), SIMD[DType.float64, 1](1.0 / prime)))
            j += 1


def _ftables(index: Int, n: Int) -> _FTables:
    """The tables of float prime `index` for n = 2^k or 3 * 2^k points."""
    var primes = materialize[_FPRIMES]()
    var roots = materialize[_FROOTS]()
    var cube_roots = materialize[_FCUBE_ROOTS]()
    return _FTables(primes[index], roots[index], n, cube_roots[index])


def _ftransform_length(count: Int) -> Int:
    """The least 2^k or 3 * 2^k (with 2^k >= 16) holding count points."""
    var n = 1
    while n < count:
        n *= 2
    var third = 16
    while 3 * third < count:
        third *= 2
    return min(n, 3 * third)


comptime _BLOCK = 1 << 14
"""Stages whose butterflies stay inside blocks of this many values run block by block, in cache."""


@always_inline
def _forward_stages(ptr: Pointer[mut=True, Float64, _], length: Int, top: Int, stop: Int, wp: Pointer[mut=False, Float64, _], prime: Float64):
    """Decimation-in-frequency stages m = top, top / 2, ... down to just above `stop` (and at least _W) over `length` values."""
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var m = top
    # Stages m and m / 2 together while both are vector-wide: one pass loads
    # four vectors and applies the same four butterflies in the same order.
    while m >= 2 * _W and m // 2 > stop:
        var h = m // 2
        var s = 0
        while s < length:
            var j = 0
            while j < h:
                var a0 = ptr.unsafe_offset(s + j).unsafe_load[width=_W]()
                var a1 = ptr.unsafe_offset(s + j + h).unsafe_load[width=_W]()
                var a2 = ptr.unsafe_offset(s + j + m).unsafe_load[width=_W]()
                var a3 = ptr.unsafe_offset(s + j + m + h).unsafe_load[width=_W]()
                var b0 = _reduce[_W](a0 + a2, p, inv)
                var b2 = _mulmod[_W](a0 - a2, wp.unsafe_offset(m + j).unsafe_load[width=_W](), p, inv)
                var b1 = _reduce[_W](a1 + a3, p, inv)
                var b3 = _mulmod[_W](a1 - a3, wp.unsafe_offset(m + j + h).unsafe_load[width=_W](), p, inv)
                var w = wp.unsafe_offset(h + j).unsafe_load[width=_W]()
                ptr.unsafe_offset(s + j).unsafe_store(_reduce[_W](b0 + b1, p, inv))
                ptr.unsafe_offset(s + j + h).unsafe_store(_mulmod[_W](b0 - b1, w, p, inv))
                ptr.unsafe_offset(s + j + m).unsafe_store(_reduce[_W](b2 + b3, p, inv))
                ptr.unsafe_offset(s + j + m + h).unsafe_store(_mulmod[_W](b2 - b3, w, p, inv))
                j += _W
            s += 2 * m
        m //= 4
    while m >= _W and m > stop:
        var s = 0
        while s < length:
            var j = 0
            while j < m:
                var u = ptr.unsafe_offset(s + j).unsafe_load[width=_W]()
                var v = ptr.unsafe_offset(s + j + m).unsafe_load[width=_W]()
                var w = wp.unsafe_offset(m + j).unsafe_load[width=_W]()
                ptr.unsafe_offset(s + j).unsafe_store(_reduce[_W](u + v, p, inv))
                ptr.unsafe_offset(s + j + m).unsafe_store(_mulmod[_W](u - v, w, p, inv))
                j += _W
            s += 2 * m
        m //= 2


@always_inline
def _forward_tail(ptr: Pointer[mut=True, Float64, _], length: Int, wp: Pointer[mut=False, Float64, _], prime: Float64):
    """The last three stages inside each vector of eight."""
    var p4 = SIMD[DType.float64, 4](prime)
    var inv4 = SIMD[DType.float64, 4](1.0 / prime)
    var w4 = wp.unsafe_offset(4).unsafe_load[width=4]()
    var w2 = wp.unsafe_offset(2).unsafe_load[width=2]()
    var w22 = w2.join(w2)
    var s = 0
    while s < length:
        var x = ptr.unsafe_offset(s).unsafe_load[width=_W]()
        # m = 4: lanes 0-3 against 4-7.
        var u = x.slice[4, offset=0]()
        var v = x.slice[4, offset=4]()
        x = _reduce[4](u + v, p4, inv4).join(_mulmod[4](u - v, w4, p4, inv4))
        # m = 2: pairs 0-1 against 2-3 in each half.
        var y = x.shuffle[0, 1, 4, 5, 2, 3, 6, 7]()
        u = y.slice[4, offset=0]()
        v = y.slice[4, offset=4]()
        x = _reduce[4](u + v, p4, inv4).join(_mulmod[4](u - v, w22, p4, inv4)).shuffle[0, 1, 4, 5, 2, 3, 6, 7]()
        # m = 1: even lanes against odd lanes; the twiddle is one.
        y = x.shuffle[0, 2, 4, 6, 1, 3, 5, 7]()
        u = y.slice[4, offset=0]()
        v = y.slice[4, offset=4]()
        x = _reduce[4](u + v, p4, inv4).join(_reduce[4](u - v, p4, inv4)).shuffle[0, 4, 1, 5, 2, 6, 3, 7]()
        ptr.unsafe_offset(s).unsafe_store(x)
        s += _W


def _fforward(mut a: List[Float64], t: _FTables, prime: Float64):
    """Decimation in frequency; bit-reversed output (within each third for 3 * 2^k points), values within 0.875 p.

    Every butterfly is the same in any grouping of stages, so blocking only
    reorders independent work.
    """
    var ptr = a.unsafe_ptr()
    var wp = t.forward.unsafe_ptr()
    var m = t.third
    if m:
        _forward_radix3(ptr, m, t.radix.unsafe_ptr(), t.omega, prime)
        for s in range(3):
            _forward_power(ptr.unsafe_offset(s * m), m, wp, prime)
        return
    _forward_power(ptr, len(a), wp, prime)


def _forward_radix3(ptr: Pointer[mut=True, Float64, _], m: Int, rp: Pointer[mut=False, Float64, _], omega: Float64, prime: Float64):
    """y_t[j] = w^(t j) sum_s omega^(s t) x[j + s m], for thirds t of a 3 m-point transform; m is a multiple of _W."""
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var o = SIMD[DType.float64, _W](omega)
    var j = 0
    while j < m:
        var x0 = ptr.unsafe_offset(j).unsafe_load[width=_W]()
        var x1 = ptr.unsafe_offset(j + m).unsafe_load[width=_W]()
        var x2 = ptr.unsafe_offset(j + 2 * m).unsafe_load[width=_W]()
        # omega^2 = -1 - omega: y_1 = (x0 - x2) + omega (x1 - x2), y_2 = (x0 - x1) - omega (x1 - x2).
        var t = _mulmod[_W](x1 - x2, o, p, inv)
        ptr.unsafe_offset(j).unsafe_store(_reduce[_W](x0 + x1 + x2, p, inv))
        var y1 = _reduce[_W](x0 - x2 + t, p, inv)
        var y2 = _reduce[_W](x0 - x1 - t, p, inv)
        ptr.unsafe_offset(j + m).unsafe_store(_mulmod[_W](y1, rp.unsafe_offset(j).unsafe_load[width=_W](), p, inv))
        ptr.unsafe_offset(j + 2 * m).unsafe_store(_mulmod[_W](y2, rp.unsafe_offset(m + j).unsafe_load[width=_W](), p, inv))
        j += _W


def _inverse_radix3(ptr: Pointer[mut=True, Float64, _], m: Int, rp: Pointer[mut=False, Float64, _], omega: Float64, prime: Float64):
    """x[j + s m] = sum_t omega^(-s t) w^(-t j) z_t[j], undoing _forward_radix3 up to the factor 3."""
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var o = SIMD[DType.float64, _W](omega)
    var j = 0
    while j < m:
        var y0 = ptr.unsafe_offset(j).unsafe_load[width=_W]()
        var y1 = _mulmod[_W](ptr.unsafe_offset(j + m).unsafe_load[width=_W](), rp.unsafe_offset(2 * m + j).unsafe_load[width=_W](), p, inv)
        var y2 = _mulmod[_W](ptr.unsafe_offset(j + 2 * m).unsafe_load[width=_W](), rp.unsafe_offset(3 * m + j).unsafe_load[width=_W](), p, inv)
        # omega^-1 = omega^2: x_1 = (y0 - y1) + omega (y2 - y1), x_2 = (y0 - y2) - omega (y2 - y1).
        var t = _mulmod[_W](y2 - y1, o, p, inv)
        ptr.unsafe_offset(j).unsafe_store(_reduce[_W](y0 + y1 + y2, p, inv))
        ptr.unsafe_offset(j + m).unsafe_store(_reduce[_W](y0 - y1 + t, p, inv))
        ptr.unsafe_offset(j + 2 * m).unsafe_store(_reduce[_W](y0 - y2 - t, p, inv))
        j += _W


def _forward_power(ptr: Pointer[mut=True, Float64, _], n: Int, wp: Pointer[mut=False, Float64, _], prime: Float64):
    """The decimation-in-frequency transform of n = 2^k values at ptr."""
    if n >= 2 * _W:
        var block = min(n, _BLOCK)
        _forward_stages(ptr, n, n // 2, block // 2, wp, prime)
        var b = 0
        while b < n:
            _forward_stages(ptr.unsafe_offset(b), block, block // 2, 0, wp, prime)
            _forward_tail(ptr.unsafe_offset(b), block, wp, prime)
            b += block
        return
    var p1 = SIMD[DType.float64, 1](prime)
    var inv1 = SIMD[DType.float64, 1](1.0 / prime)
    var m = n // 2
    while m >= 1:
        var s = 0
        while s < n:
            for j in range(m):
                var u = ptr.unsafe_offset(s + j).unsafe_load[width=1]()
                var v = ptr.unsafe_offset(s + j + m).unsafe_load[width=1]()
                var w = wp.unsafe_offset(m + j).unsafe_load[width=1]()
                ptr.unsafe_offset(s + j).unsafe_store(_reduce[1](u + v, p1, inv1))
                ptr.unsafe_offset(s + j + m).unsafe_store(_mulmod[1](u - v, w, p1, inv1))
            s += 2 * m
        m //= 2


@always_inline
def _inverse_tail(ptr: Pointer[mut=True, Float64, _], length: Int, wp: Pointer[mut=False, Float64, _], prime: Float64):
    """The first three stages inside each vector of eight."""
    var p4 = SIMD[DType.float64, 4](prime)
    var inv4 = SIMD[DType.float64, 4](1.0 / prime)
    var w4 = wp.unsafe_offset(4).unsafe_load[width=4]()
    var w2 = wp.unsafe_offset(2).unsafe_load[width=2]()
    var w22 = w2.join(w2)
    var s = 0
    while s < length:
        var x = ptr.unsafe_offset(s).unsafe_load[width=_W]()
        # m = 1: even lanes against odd lanes; the twiddle is one.
        var y = x.shuffle[0, 2, 4, 6, 1, 3, 5, 7]()
        var u = y.slice[4, offset=0]()
        var v = y.slice[4, offset=4]()
        x = _reduce[4](u + v, p4, inv4).join(_reduce[4](u - v, p4, inv4)).shuffle[0, 4, 1, 5, 2, 6, 3, 7]()
        # m = 2.
        y = x.shuffle[0, 1, 4, 5, 2, 3, 6, 7]()
        u = y.slice[4, offset=0]()
        var tw = _mulmod[4](y.slice[4, offset=4](), w22, p4, inv4)
        x = _reduce[4](u + tw, p4, inv4).join(_reduce[4](u - tw, p4, inv4)).shuffle[0, 1, 4, 5, 2, 3, 6, 7]()
        # m = 4.
        u = x.slice[4, offset=0]()
        tw = _mulmod[4](x.slice[4, offset=4](), w4, p4, inv4)
        ptr.unsafe_offset(s).unsafe_store(_reduce[4](u + tw, p4, inv4).join(_reduce[4](u - tw, p4, inv4)))
        s += _W


@always_inline
def _inverse_stages(ptr: Pointer[mut=True, Float64, _], length: Int, bottom: Int, stop: Int, wp: Pointer[mut=False, Float64, _], prime: Float64):
    """Decimation-in-time stages m = bottom, 2 bottom, ... below `stop` over `length` values; bottom >= _W."""
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var m = bottom
    # Stages m and 2 m together: four vectors per pass, butterflies unchanged.
    while 4 * m <= stop:
        var d = 2 * m
        var s = 0
        while s < length:
            var j = 0
            while j < m:
                var a0 = ptr.unsafe_offset(s + j).unsafe_load[width=_W]()
                var a1 = ptr.unsafe_offset(s + j + m).unsafe_load[width=_W]()
                var a2 = ptr.unsafe_offset(s + j + d).unsafe_load[width=_W]()
                var a3 = ptr.unsafe_offset(s + j + d + m).unsafe_load[width=_W]()
                var w = wp.unsafe_offset(m + j).unsafe_load[width=_W]()
                var t1 = _mulmod[_W](a1, w, p, inv)
                var b0 = _reduce[_W](a0 + t1, p, inv)
                var b1 = _reduce[_W](a0 - t1, p, inv)
                var t3 = _mulmod[_W](a3, w, p, inv)
                var b2 = _reduce[_W](a2 + t3, p, inv)
                var b3 = _reduce[_W](a2 - t3, p, inv)
                var u2 = _mulmod[_W](b2, wp.unsafe_offset(d + j).unsafe_load[width=_W](), p, inv)
                ptr.unsafe_offset(s + j).unsafe_store(_reduce[_W](b0 + u2, p, inv))
                ptr.unsafe_offset(s + j + d).unsafe_store(_reduce[_W](b0 - u2, p, inv))
                var u3 = _mulmod[_W](b3, wp.unsafe_offset(d + j + m).unsafe_load[width=_W](), p, inv)
                ptr.unsafe_offset(s + j + m).unsafe_store(_reduce[_W](b1 + u3, p, inv))
                ptr.unsafe_offset(s + j + d + m).unsafe_store(_reduce[_W](b1 - u3, p, inv))
                j += _W
            s += 2 * d
        m *= 4
    while m < stop:
        var s = 0
        while s < length:
            var j = 0
            while j < m:
                var u = ptr.unsafe_offset(s + j).unsafe_load[width=_W]()
                var t1 = _mulmod[_W](ptr.unsafe_offset(s + j + m).unsafe_load[width=_W](), wp.unsafe_offset(m + j).unsafe_load[width=_W](), p, inv)
                ptr.unsafe_offset(s + j).unsafe_store(_reduce[_W](u + t1, p, inv))
                ptr.unsafe_offset(s + j + m).unsafe_store(_reduce[_W](u - t1, p, inv))
                j += _W
            s += 2 * m
        m *= 2


def _finverse(mut a: List[Float64], t: _FTables, prime: Float64):
    """Decimation in time without the 1/n factor; bit-reversed input (within each third for 3 * 2^k points), values within p / 2 + 1."""
    var ptr = a.unsafe_ptr()
    var wp = t.inverse.unsafe_ptr()
    var m = t.third
    if m:
        for s in range(3):
            _inverse_power(ptr.unsafe_offset(s * m), m, wp, prime)
        _inverse_radix3(ptr, m, t.radix.unsafe_ptr(), t.omega, prime)
        return
    _inverse_power(ptr, len(a), wp, prime)


def _inverse_power(ptr: Pointer[mut=True, Float64, _], n: Int, wp: Pointer[mut=False, Float64, _], prime: Float64):
    """The decimation-in-time transform of n = 2^k values at ptr, without the 1/n factor."""
    if n >= 2 * _W:
        var block = min(n, _BLOCK)
        var b = 0
        while b < n:
            _inverse_tail(ptr.unsafe_offset(b), block, wp, prime)
            _inverse_stages(ptr.unsafe_offset(b), block, _W, block, wp, prime)
            b += block
        if n > block:
            _inverse_stages(ptr, n, block, n, wp, prime)
        return
    var p1 = SIMD[DType.float64, 1](prime)
    var inv1 = SIMD[DType.float64, 1](1.0 / prime)
    var m = 1
    while m < n:
        var s = 0
        while s < n:
            for j in range(m):
                var u = ptr.unsafe_offset(s + j).unsafe_load[width=1]()
                var t1 = _mulmod[1](ptr.unsafe_offset(s + j + m).unsafe_load[width=1](), wp.unsafe_offset(m + j).unsafe_load[width=1](), p1, inv1)
                ptr.unsafe_offset(s + j).unsafe_store(_reduce[1](u + t1, p1, inv1))
                ptr.unsafe_offset(s + j + m).unsafe_store(_reduce[1](u - t1, p1, inv1))
            s += 2 * m
        m *= 2


def _floats(values: List[UInt64], n: Int) -> List[Float64]:
    """Residues below the prime as floats, zero-padded to n."""
    var out = List[Float64](length=n, fill=0)
    var count = min(n, len(values))
    var src = values.unsafe_ptr()
    var dst = out.unsafe_ptr()
    var i = 0
    while i + _W <= count:
        dst.unsafe_offset(i).unsafe_store(src.unsafe_offset(i).unsafe_load[width=_W]().cast[DType.float64]())
        i += _W
    while i < count:
        dst.unsafe_offset(i).unsafe_store(Float64(src.unsafe_offset(i).unsafe_load()))
        i += 1
    return out^


def _ftransform(values: List[UInt64], tables: _FTables, prime: Float64, n: Int) -> List[Float64]:
    """The forward transform of residues below the prime, zero-padded to n."""
    var x = _floats(values, n)
    _fforward(x, tables, prime)
    return x^


def _fmultiply(x: List[Float64], y: List[Float64], tables: _FTables, prime: Float64, count: Int) -> List[UInt64]:
    """The first `count` cyclic convolution coefficients of two forward transforms, reduced below the prime."""
    var n = len(x)
    var z = List[Float64](length=n, fill=0)
    var px = x.unsafe_ptr()
    var py = y.unsafe_ptr()
    var pz = z.unsafe_ptr()
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var p1 = SIMD[DType.float64, 1](prime)
    var inv1 = SIMD[DType.float64, 1](1.0 / prime)
    var i = 0
    while i + _W <= n:
        pz.unsafe_offset(i).unsafe_store(_mulmod[_W](px.unsafe_offset(i).unsafe_load[width=_W](), py.unsafe_offset(i).unsafe_load[width=_W](), p, inv))
        i += _W
    while i < n:
        pz.unsafe_offset(i).unsafe_store(_mulmod[1](px.unsafe_offset(i).unsafe_load[width=1](), py.unsafe_offset(i).unsafe_load[width=1](), p1, inv1))
        i += 1
    _finverse(z, tables, prime)
    return _fscaled(z, 0, count, prime)


def _faccumulate(mut acc: List[Float64], x: List[Float64], y: List[Float64], prime: Float64):
    """acc += x * y pointwise for forward transforms, each entry kept within
    prime / 2 + 1, as _fmultiply's products are before their inverse."""
    var n = len(acc)
    var pa = acc.unsafe_ptr()
    var px = x.unsafe_ptr()
    var py = y.unsafe_ptr()
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var p1 = SIMD[DType.float64, 1](prime)
    var inv1 = SIMD[DType.float64, 1](1.0 / prime)
    var i = 0
    while i + _W <= n:
        var product = _mulmod[_W](px.unsafe_offset(i).unsafe_load[width=_W](), py.unsafe_offset(i).unsafe_load[width=_W](), p, inv)
        pa.unsafe_offset(i).unsafe_store(_reduce[_W](pa.unsafe_offset(i).unsafe_load[width=_W]() + product, p, inv))
        i += _W
    while i < n:
        var product = _mulmod[1](px.unsafe_offset(i).unsafe_load[width=1](), py.unsafe_offset(i).unsafe_load[width=1](), p1, inv1)
        pa.unsafe_offset(i).unsafe_store(_reduce[1](pa.unsafe_offset(i).unsafe_load[width=1]() + product, p1, inv1))
        i += 1


def _fcolumn(x: List[Float64], var y: List[Float64], tables: _FTables, prime: Float64, lo: Int, hi: Int) -> List[UInt64]:
    """Entries lo ..< hi of the cyclic convolution of a forward transform x with the points y, reduced below the prime."""
    _fforward(y, tables, prime)
    var n = len(x)
    var px = x.unsafe_ptr()
    var py = y.unsafe_ptr()
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var i = 0
    while i + _W <= n:
        py.unsafe_offset(i).unsafe_store(_mulmod[_W](px.unsafe_offset(i).unsafe_load[width=_W](), py.unsafe_offset(i).unsafe_load[width=_W](), p, inv))
        i += _W
    while i < n:
        py.unsafe_offset(i).unsafe_store(_mulmod[1](px.unsafe_offset(i).unsafe_load[width=1](), py.unsafe_offset(i).unsafe_load[width=1](), SIMD[DType.float64, 1](prime), SIMD[DType.float64, 1](1.0 / prime)))
        i += 1
    _finverse(y, tables, prime)
    return _fscaled(y, lo, hi, prime)


def _fscaled(z: List[Float64], lo: Int, hi: Int, prime: Float64) -> List[UInt64]:
    """Entries lo ..< hi of an unscaled inverse transform, divided by its length and moved into [0, q)."""
    var n = len(z)
    var pz = z.unsafe_ptr()
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var p1 = SIMD[DType.float64, 1](prime)
    var inv1 = SIMD[DType.float64, 1](1.0 / prime)
    var q = UInt64(prime)
    var scale = Float64(_inverse_length(n, q))
    var sv = SIMD[DType.float64, _W](scale)
    var result = List[UInt64](length=hi - lo, fill=0)
    var out = result.unsafe_ptr()
    var zero = SIMD[DType.float64, _W](0)
    var i = lo
    while i + _W <= hi:
        var r = _reduce[_W](_mulmod[_W](pz.unsafe_offset(i).unsafe_load[width=_W](), sv, p, inv), p, inv)
        r = r.lt(zero).select(r + p, r)
        out.unsafe_offset(i - lo).unsafe_store(r.cast[DType.uint64]())
        i += _W
    var s1 = SIMD[DType.float64, 1](scale)
    while i < hi:
        var r = _reduce[1](_mulmod[1](pz.unsafe_offset(i).unsafe_load[width=1](), s1, p1, inv1), p1, inv1)[0]
        if r < 0:
            r += prime
        out.unsafe_offset(i - lo).unsafe_store(UInt64(r))
        i += 1
    return result^


def _fconvolve(a: List[UInt64], b: List[UInt64], count: Int, square: Bool, index: Int, n: Int) -> List[UInt64]:
    """The first `count` cyclic convolution coefficients of residues below float prime `index`, reduced below it.

    The transform length n is 2^k or 3 * 2^k (2^k >= 16); operands are read up to n.
    """
    var q = materialize[_FPRIMES]()[index]
    var prime = Float64(q)
    var tables = _ftables(index, n)
    var x = _ftransform(a, tables, prime, n)
    if square:
        return _fmultiply(x, x, tables, prime, count)
    return _fmultiply(x, _ftransform(b, tables, prime, n), tables, prime, count)


def _fpoints(low: List[UInt64], high: List[UInt64], q: UInt64, n: Int) -> List[Float64]:
    """Signed digits as residues in [0, q), zero-padded to n points: low holds bits 0-63, high bits 64-119 with the sign in bit 63.

    With low = l1 2^32 + l0 and high = h1 2^32 + h0, the residue is
    l1 2^32 + l0 + h0 (2^64 mod q) + h1 (2^96 mod q): products of factors
    below q (M < 1), sums below 4 q reduced, then moved into [0, q) and
    negated when the sign is set.
    """
    var out = List[Float64](length=n, fill=0)
    var count = len(low)
    var lp = low.unsafe_ptr()
    var hp = high.unsafe_ptr()
    var op = out.unsafe_ptr()
    var prime = Float64(q)
    var p = SIMD[DType.float64, _W](prime)
    var inverse = SIMD[DType.float64, _W](1.0 / prime)
    var c32 = SIMD[DType.float64, _W](Float64(UInt64(1) << 32))
    var c64 = SIMD[DType.float64, _W](Float64(UInt64((UInt128(1) << 64) % UInt128(q))))
    var c96 = SIMD[DType.float64, _W](Float64(UInt64((UInt128(1) << 96) % UInt128(q))))
    var half = SIMD[DType.uint64, _W]((UInt64(1) << 32) - 1)
    var magnitude = SIMD[DType.uint64, _W]((UInt64(1) << 63) - 1)
    var zero = SIMD[DType.float64, _W](0)
    var i = 0
    while i < count:
        var l = SIMD[DType.uint64, _W](0)
        var h = SIMD[DType.uint64, _W](0)
        if i + _W <= count:
            l = lp.unsafe_offset(i).unsafe_load[width=_W]()
            h = hp.unsafe_offset(i).unsafe_load[width=_W]()
        else:
            for k in range(count - i):
                l[k] = lp.unsafe_offset(i + k).unsafe_load()
                h[k] = hp.unsafe_offset(i + k).unsafe_load()
        var negative = (h >> 63).ne(SIMD[DType.uint64, _W](0))
        h = h & magnitude
        var r = _reduce[_W](_mulmod[_W]((l >> 32).cast[DType.float64](), c32, p, inverse) + (l & half).cast[DType.float64](), p, inverse)
        var t = _reduce[_W](_mulmod[_W]((h & half).cast[DType.float64](), c64, p, inverse) + _mulmod[_W]((h >> 32).cast[DType.float64](), c96, p, inverse), p, inverse)
        r = _canonical[_W](r + t, p, inverse)
        r = (negative & r.ne(zero)).select(p - r, r)
        if i + _W <= count:
            op.unsafe_offset(i).unsafe_store(r)
        else:
            for k in range(count - i):
                op.unsafe_offset(i + k).unsafe_store(r[k])
        i += _W
    return out^


def _fconvolve_floats(var x: List[Float64], var y: List[Float64], lo: Int, hi: Int, square: Bool, index: Int) -> List[UInt64]:
    """Entries lo ..< hi of the cyclic convolution of residues below float prime `index`, as floats, reduced below it.

    The operands hold the transform length of values, 2^k or 3 * 2^k; the
    products run in place, with the same arithmetic as _fconvolve.
    """
    var q = materialize[_FPRIMES]()[index]
    var prime = Float64(q)
    var n = len(x)
    var tables = _ftables(index, n)
    _fforward(x, tables, prime)
    if not square:
        _fforward(y, tables, prime)
    var px = x.unsafe_ptr()
    var py = x.unsafe_ptr() if square else y.unsafe_ptr()
    var p = SIMD[DType.float64, _W](prime)
    var inv = SIMD[DType.float64, _W](1.0 / prime)
    var i = 0
    while i + _W <= n:
        px.unsafe_offset(i).unsafe_store(_mulmod[_W](px.unsafe_offset(i).unsafe_load[width=_W](), py.unsafe_offset(i).unsafe_load[width=_W](), p, inv))
        i += _W
    while i < n:
        px.unsafe_offset(i).unsafe_store(_mulmod[1](px.unsafe_offset(i).unsafe_load[width=1](), py.unsafe_offset(i).unsafe_load[width=1](), SIMD[DType.float64, 1](prime), SIMD[DType.float64, 1](1.0 / prime)))
        i += 1
    _finverse(x, tables, prime)
    return _fscaled(x, lo, hi, prime)


@always_inline
def _canonical[w: Int](x: SIMD[DType.float64, w], p: SIMD[DType.float64, w], inverse: SIMD[DType.float64, w]) -> SIMD[DType.float64, w]:
    """x mod p in [0, p), for |x| < 4 p."""
    var r = _reduce[w](x, p, inverse)
    return r.lt(SIMD[DType.float64, w](0)).select(r + p, r)


@always_inline
def _digit_column[w: Int, i: Int, K: Int](r: Pointer[mut=False, UInt64, _], d: Pointer[mut=True, UInt64, _], at: Int, count: Int, p: SIMD[DType.float64, w], inverse: SIMD[DType.float64, w], prefix: SIMD[DType.float64, w], remainder: Array[SIMD[DType.float64, w], K]):
    """Digit i of the w points from `at`, from their residue modulo q_i and their digits below i."""
    # acc = v_0 + q_0 v_1 + ... + q_0 ... q_(i-2) v_(i-1) mod q_i, by Horner.
    var acc = d.unsafe_offset((i - 1) * count + at).unsafe_load[width=w]().cast[DType.float64]()
    comptime for step in range(i - 1):
        comptime j = i - 2 - step
        acc = _reduce[w](_mulmod[w](acc, remainder[j], p, inverse) + d.unsafe_offset(j * count + at).unsafe_load[width=w]().cast[DType.float64](), p, inverse)
    var x = r.unsafe_offset(at).unsafe_load[width=w]().cast[DType.float64]()
    d.unsafe_offset(i * count + at).unsafe_store(_canonical[w](_mulmod[w](x - acc, prefix, p, inverse), p, inverse).cast[DType.uint64]())


@always_inline
def _digit_columns[i: Int, K: Int](residues: List[List[UInt64]], mut digits: List[UInt64], count: Int, q: Float64, prefix: Float64, remainder: Array[Float64, K * K]):
    """Digit i of every point, eight lanes wide: one pass per digit, with its constants and pointers hoisted."""
    var r = residues[i].unsafe_ptr()
    var d = digits.unsafe_ptr()
    var p = SIMD[DType.float64, _W](q)
    var inverse = SIMD[DType.float64, _W](1.0 / q)
    var c = SIMD[DType.float64, _W](prefix)
    var lanes = Array[SIMD[DType.float64, _W], K](fill=0)
    var ones = Array[SIMD[DType.float64, 1], K](fill=0)
    comptime for j in range(K):
        lanes[j] = SIMD[DType.float64, _W](remainder[i * K + j])
        ones[j] = SIMD[DType.float64, 1](remainder[i * K + j])
    var t = 0
    while t + _W <= count:
        _digit_column[_W, i, K](r, d, t, count, p, inverse, c, lanes)
        t += _W
    while t < count:
        _digit_column[1, i, K](r, d, t, count, SIMD[DType.float64, 1](q), SIMD[DType.float64, 1](1.0 / q), SIMD[DType.float64, 1](prefix), ones)
        t += 1


def _float_digits[K: Int](residues: List[List[UInt64]], moduli: List[UInt64], mut digits: List[UInt64]):
    """The digits v_i < q_i of x = v_0 + q_0 v_1 + ... + q_0 ... q_(K-2) v_(K-1), for residues modulo K float primes.

    Digit i of point t lands at digits[i * count + t], for count points. The
    moduli are the first K float primes, which descend, so q_j mod q_i =
    q_j - q_i is tiny for j < i, and the Horner products have M < 0.01;
    r_i - acc lies within 1.5 p + 1, so the last product has M <= 1.51 (the
    module docstring has the bounds). Every reduction is of a value below 3 p,
    so all of it is exact, and each digit is the unique representative in
    [0, q_i): the same digits _digits gives.
    """
    var q = Array[Float64, K](fill=0)
    var prefix = Array[Float64, K](fill=0)
    var remainder = Array[Float64, K * K](fill=0)
    for i in range(K):
        q[i] = Float64(moduli[i])
        var product = UInt64(1)
        for j in range(i):
            var r = moduli[j] % moduli[i]
            remainder[i * K + j] = Float64(r)
            product = _fmul(product, r, moduli[i])
        prefix[i] = Float64(_fpower(product, moduli[i] - 2, moduli[i]))
    var count = len(residues[0])
    var d = digits.unsafe_ptr()
    var r0 = residues[0].unsafe_ptr()
    for t in range(count):
        d.unsafe_offset(t).unsafe_store(r0.unsafe_offset(t).unsafe_load())
    comptime for i in range(1, K):
        _digit_columns[i, K](residues, digits, count, q[i], prefix[i], remainder)
