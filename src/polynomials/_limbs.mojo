"""Montgomery arithmetic on three or four 64-bit limbs, for odd moduli below 2^254.

Residues are SIMD[uint64, 4] limbs, lowest first, in Montgomery form with
R = 2^(64 N). Products multiply by coarsely integrated operand scanning;
sums of products can instead accumulate unreduced and be reduced once
while they stay below m R. Shared by Hensel lifting and Wang's sparse
interpolation, whose moduli outgrow 128 bits.
"""
from apn_mojo import Integer, div_rem_euclid
from .coefficients import ZZ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _product, _charge_size

comptime _LOW = UInt128(UInt64.MAX)


comptime _Limbs = SIMD[DType.uint64, 4]
"""Up to four 64-bit limbs, lowest first; limbs at and above the width stay zero."""


@fieldwise_init
struct _Wide[N: Int](Copyable, Movable):
    """Montgomery arithmetic modulo an odd m < 2^(64 N - 2) with R = 2^(64 N)."""
    var m: _Limbs
    var minv: UInt64
    """-m^-1 modulo 2^64."""
    var r2: _Limbs
    """R^2 modulo m."""
    var one: _Limbs
    """R modulo m."""
    var bits: Int
    """The bit length of m."""


@always_inline
def _limbs_less[N: Int](a: _Limbs, b: _Limbs) -> Bool:
    comptime for k in range(N):
        comptime i = N - 1 - k
        if a[i] != b[i]:
            return a[i] < b[i]
    return False


@always_inline
def _limbs_sub[N: Int](a: _Limbs, b: _Limbs) -> _Limbs:
    """a - b modulo 2^(64 N)."""
    var out = _Limbs(0)
    var borrow = UInt64(0)
    comptime for i in range(N):
        var x = a[i]
        var y = b[i]
        var d = x - y - borrow
        borrow = UInt64(1) if (x < y or (x == y and borrow == 1)) else UInt64(0)
        out[i] = d
    return out


@always_inline
def _limbs_add[N: Int](a: _Limbs, b: _Limbs) -> Tuple[_Limbs, UInt64]:
    var out = _Limbs(0)
    var carry = UInt64(0)
    comptime for i in range(N):
        var t = UInt128(a[i]) + UInt128(b[i]) + UInt128(carry)
        out[i] = UInt64(t & _LOW)
        carry = UInt64(t >> 64)
    return (out, carry)


@always_inline
def _wadd_mod[N: Int](a: _Limbs, b: _Limbs, ctx: _Wide[N]) -> _Limbs:
    var s = _limbs_add[N](a, b)
    if s[1] != 0 or not _limbs_less[N](s[0], ctx.m):
        return _limbs_sub[N](s[0], ctx.m)
    return s[0]


@always_inline
def _wsub_mod[N: Int](a: _Limbs, b: _Limbs, ctx: _Wide[N]) -> _Limbs:
    if _limbs_less[N](a, b):
        return _limbs_add[N](_limbs_sub[N](a, b), ctx.m)[0]
    return _limbs_sub[N](a, b)


@always_inline
def _wmont[N: Int](a: _Limbs, b: _Limbs, ctx: _Wide[N]) -> _Limbs:
    """a b R^-1 modulo m by coarsely integrated operand scanning."""
    var t = SIMD[DType.uint64, 8](0)
    comptime for i in range(N):
        var carry = UInt128(0)
        comptime for j in range(N):
            var x = UInt128(t[j]) + UInt128(a[j]) * UInt128(b[i]) + carry
            t[j] = UInt64(x & _LOW)
            carry = x >> 64
        var x = UInt128(t[N]) + carry
        t[N] = UInt64(x & _LOW)
        t[N + 1] = UInt64(x >> 64)
        var u = t[0] * ctx.minv
        carry = (UInt128(t[0]) + UInt128(u) * UInt128(ctx.m[0])) >> 64
        comptime for j in range(1, N):
            var y = UInt128(t[j]) + UInt128(u) * UInt128(ctx.m[j]) + carry
            t[j - 1] = UInt64(y & _LOW)
            carry = y >> 64
        var z = UInt128(t[N]) + carry
        t[N - 1] = UInt64(z & _LOW)
        t[N] = t[N + 1] + UInt64(z >> 64)
    var out = _Limbs(0)
    comptime for i in range(N):
        out[i] = t[i]
    if t[N] != 0 or not _limbs_less[N](out, ctx.m):
        return _limbs_sub[N](out, ctx.m)
    return out


def _to_limbs(a: Integer) -> _Limbs:
    """A nonnegative integer below 2^256."""
    var out = _Limbs(0)
    for i in range(min(a._word_count(), 8)):
        out[i // 2] |= UInt64(a._word(i)) << UInt64(32 * (i % 2))
    return out


def _from_limbs[N: Int](a: _Limbs, negative: Bool) -> Integer:
    var words = List[UInt32](length=2 * N, fill=0)
    for i in range(N):
        words[2 * i] = UInt32(a[i] & 0xFFFFFFFF)
        words[2 * i + 1] = UInt32(a[i] >> 32)
    return Integer._from_words(words^, negative)


def _wide_context[N: Int](modulus: Integer) raises PolyError -> _Wide[N]:
    var m = _to_limbs(modulus)
    var inverse = m[0]
    for _ in range(6):
        inverse = inverse * (UInt64(2) - m[0] * inverse)
    var r2: Integer
    try:
        r2 = div_rem_euclid(Integer(1) << (128 * N), modulus)[1].copy()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var ctx = _Wide[N](m, UInt64(0) - inverse, _to_limbs(r2), _Limbs(0), modulus.magnitude_bit_length())
    var unit = _Limbs(0)
    unit[0] = 1
    ctx.one = _wmont[N](unit, ctx.r2, ctx)
    return ctx^


comptime _Acc = SIMD[DType.uint64, 16]
"""An unreduced sum of limb products, lowest limb first."""


@always_inline
def _acc_mac[N: Int](mut acc: _Acc, a: _Limbs, b: _Limbs):
    """acc += a b."""
    comptime for i in range(N):
        var carry = UInt128(0)
        comptime for j in range(N):
            var x = UInt128(acc[i + j]) + UInt128(a[i]) * UInt128(b[j]) + carry
            acc[i + j] = UInt64(x & _LOW)
            carry = x >> 64
        comptime for k in range(i + N, 2 * N + 2):
            var y = UInt128(acc[k]) + carry
            acc[k] = UInt64(y & _LOW)
            carry = y >> 64


@always_inline
def _acc_redc[N: Int](var acc: _Acc, ctx: _Wide[N]) -> _Limbs:
    """acc R^-1 modulo m for acc < m R: the sum over R plus less than m, below 2 m."""
    comptime for i in range(N):
        var u = acc[i] * ctx.minv
        var carry = UInt128(0)
        comptime for j in range(N):
            var x = UInt128(acc[i + j]) + UInt128(u) * UInt128(ctx.m[j]) + carry
            acc[i + j] = UInt64(x & _LOW)
            carry = x >> 64
        comptime for k in range(i + N, 2 * N + 2):
            var y = UInt128(acc[k]) + carry
            acc[k] = UInt64(y & _LOW)
            carry = y >> 64
    var out = _Limbs(0)
    comptime for i in range(N):
        out[i] = acc[N + i]
    if acc[2 * N] != 0 or not _limbs_less[N](out, ctx.m):
        return _limbs_sub[N](out, ctx.m)
    return out


@always_inline
def _lazy_terms[N: Int](terms: Int, ctx: _Wide[N]) -> Bool:
    """Whether sums of `terms` products of residues stay below m R."""
    var length = 0
    var t = terms + 1
    while t > 0:
        length += 1
        t >>= 1
    return ctx.bits + length < 64 * N


def _ltrim(mut a: List[_Limbs]):
    while len(a) > 0 and a[len(a) - 1].reduce_or() == 0:
        _ = a.pop()


def _ladd[N: Int](a: List[_Limbs], b: List[_Limbs], ctx: _Wide[N], mut budget: Budget) raises PolyError -> List[_Limbs]:
    var n = max(len(a), len(b))
    _allocation[_Limbs](n, budget)
    var c = List[_Limbs](capacity=n)
    for i in range(n):
        c.append(_wadd_mod[N](a[i] if i < len(a) else _Limbs(0), b[i] if i < len(b) else _Limbs(0), ctx))
    _ltrim(c)
    return c^


def _lsub[N: Int](a: List[_Limbs], b: List[_Limbs], ctx: _Wide[N], mut budget: Budget) raises PolyError -> List[_Limbs]:
    var n = max(len(a), len(b))
    _allocation[_Limbs](n, budget)
    var c = List[_Limbs](capacity=n)
    for i in range(n):
        c.append(_wsub_mod[N](a[i] if i < len(a) else _Limbs(0), b[i] if i < len(b) else _Limbs(0), ctx))
    _ltrim(c)
    return c^


def _lmul[N: Int](a: List[_Limbs], b: List[_Limbs], ctx: _Wide[N], mut budget: Budget) raises PolyError -> List[_Limbs]:
    if len(a) == 0 or len(b) == 0:
        return List[_Limbs]()
    var n = len(a) + len(b) - 1
    _allocation[_Limbs](n, budget)
    _charge_size(_product(len(a), len(b)), budget)
    if _lazy_terms[N](min(len(a), len(b)), ctx):
        # Each coefficient sums its products unreduced, then one reduction.
        var lazy = List[_Limbs](capacity=n)
        for k in range(n):
            var acc = _Acc(0)
            for i in range(max(0, k - len(b) + 1), min(k + 1, len(a))):
                _acc_mac[N](acc, a[i], b[k - i])
            lazy.append(_acc_redc[N](acc, ctx))
        _ltrim(lazy)
        return lazy^
    var c = List[_Limbs](length=n, fill=_Limbs(0))
    for i in range(len(a)):
        if a[i].reduce_or() == 0:
            continue
        for j in range(len(b)):
            c[i + j] = _wadd_mod[N](c[i + j], _wmont[N](a[i], b[j], ctx), ctx)
    _ltrim(c)
    return c^


def _ldivmod[N: Int](a: List[_Limbs], h: List[_Limbs], ctx: _Wide[N], mut budget: Budget) raises PolyError -> Tuple[List[_Limbs], List[_Limbs]]:
    """Quotient and remainder by a monic h (leading residue R mod m)."""
    var n = len(h) - 1
    var r = a.copy()
    if len(r) <= n:
        return (List[_Limbs](), r^)
    var q = List[_Limbs](length=len(r) - n, fill=_Limbs(0))
    _charge_size(_product(len(q), len(h)), budget)
    if _lazy_terms[N](min(len(q), n) + 1, ctx):
        # Entries stay unreduced at scale R^2: (x R)(R mod m) to start, then
        # (c R)((m - h_i) R) for each elimination; each coefficient is
        # reduced once, when it becomes the top or the remainder.
        var negated = List[_Limbs](capacity=n)
        for i in range(n):
            negated.append(_limbs_sub[N](ctx.m, h[i]) if h[i].reduce_or() != 0 else _Limbs(0))
        var acc = List[_Acc](capacity=len(r))
        for i in range(len(r)):
            var start = _Acc(0)
            _acc_mac[N](start, r[i], ctx.one)
            acc.append(start)
        var top_lazy = len(r) - 1
        while top_lazy >= n:
            var c = _acc_redc[N](acc[top_lazy], ctx)
            if c.reduce_or() != 0:
                var shift = top_lazy - n
                q[shift] = c
                for i in range(n):
                    _acc_mac[N](acc[shift + i], c, negated[i])
            top_lazy -= 1
        var rest = List[_Limbs](capacity=n)
        for i in range(n):
            rest.append(_acc_redc[N](acc[i], ctx))
        _ltrim(rest)
        _ltrim(q)
        return (q^, rest^)
    var top = len(r) - 1
    while top >= n:
        var c = r[top]
        if c.reduce_or() != 0:
            var shift = top - n
            q[shift] = c
            for i in range(n):
                r[shift + i] = _wsub_mod[N](r[shift + i], _wmont[N](c, h[i], ctx), ctx)
        r[top] = _Limbs(0)
        top -= 1
    _ltrim(r)
    _ltrim(q)
    return (q^, r^)


def _limbs_of[N: Int](f: DensePoly[ZZ], ctx: _Wide[N], modulus: Integer, mut budget: Budget) raises PolyError -> List[_Limbs]:
    _allocation[_Limbs](len(f.coefficients), budget)
    var result = List[_Limbs](capacity=len(f.coefficients))
    for c in f.coefficients:
        budget.charge(1)
        var r: Integer
        try:
            r = div_rem_euclid(c.value, modulus)[1].copy()
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        result.append(_wmont[N](_to_limbs(r), ctx.r2, ctx))
    _ltrim(result)
    return result^


def _symmetric_limbs[N: Int](a: List[_Limbs], ctx: _Wide[N], variable: Int64, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(a), budget)
    var unit = _Limbs(0)
    unit[0] = 1
    var values = List[ZZ](capacity=len(a))
    for x in a:
        budget.charge(1)
        var v = _wmont[N](x, unit, ctx)
        var w = _limbs_sub[N](ctx.m, v)
        if _limbs_less[N](w, v):
            values.append(ZZ(_from_limbs[N](w, True)))
        else:
            values.append(ZZ(_from_limbs[N](v, False)))
    return DensePoly[ZZ].from_coefficients(variable, NoContext(), values^, budget)
