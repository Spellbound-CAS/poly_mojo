"""Quadratic multifactor Hensel lifting over a binary factor tree (Z-4)."""
from apn_mojo import Integer, inverse_mod, div_rem_euclid
from ._integer_ops import _primitive, _divide_scalar
from .coefficients import ZZ, ModInt, NoContext, _work, _pow, _submul
from .dense import DensePoly, _slice, _trim
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _product, _sum, _charge_size
from ._coefficient_ops import _copy as _copy_coefficient
from ._coefficient_ops import _quotient
from ._dense_ops import _constant, _remainder, _exact
from ._euclidean import _extended
from ._images import _mod_image
from .primes import _symmetric
from ._ntt import _ntt_low_product
from ._kronecker import _kronecker_mul
from ._limbs import _Limbs, _Wide, _wide_context, _wmont, _ladd, _lsub, _lmul, _ldivmod, _limbs_of, _symmetric_limbs

comptime _KRONECKER_BITS = 384
"""Coefficient sizes up to which short lift products pack into one integer product."""


def _hmul(a: DensePoly[ZZ], b: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """a * b, by Kronecker packing for short operands with mid-size coefficients."""
    var shorter = min(len(a.coefficients), len(b.coefficients))
    if shorter >= 3 and shorter < 32:
        var bits = 0
        for c in a.coefficients:
            bits = max(bits, c.value.magnitude_bit_length())
        for c in b.coefficients:
            bits = max(bits, c.value.magnitude_bit_length())
        if bits <= _KRONECKER_BITS:
            return _kronecker_mul(a, b, budget)
    return a.mul(b, budget)


def _lift_image(f: DensePoly[ModInt], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(f.coefficients), budget)
    var values = List[ZZ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        values.append(ZZ(Integer(Int(c.value))))
    return DensePoly[ZZ]._from_canonical(f.variable, NoContext(), values^)


def _truncate(f: DensePoly[ZZ], modulus: ZZ, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(f.coefficients), budget)
    var values = List[ZZ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        values.append(ZZ(_symmetric(c.value, modulus.value, budget)))
    return DensePoly[ZZ].from_coefficients(f.variable, NoContext(), values^, budget)


def _trial_divide(f: DensePoly[ZZ], g: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Optional[DensePoly[ZZ]]:
    f._same(g)
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "zero integer polynomial divisor")
    var n = max(f.degree() - g.degree() + 1, 0)
    _allocation[ZZ](n, budget)
    var quotient = List[ZZ](length=n, fill=ZZ.zero(f.ring))
    var remainder = _slice(_view(f.coefficients), 0, len(f.coefficients), budget)
    while len(remainder) >= len(g.coefficients):
        budget.charge(1)
        var shift = len(remainder) - len(g.coefficients)
        var scalar = _quotient(remainder[len(remainder) - 1], g.coefficients[g.degree()], budget)
        if not scalar:
            return None
        var q = scalar.take()
        for i in range(len(g.coefficients)):
            budget.charge(1)
            _submul(remainder[shift + i], q, g.coefficients[i], budget)
        quotient[shift] = q^
        _trim(remainder, budget)
    if len(remainder):
        return None
    return DensePoly[ZZ].from_coefficients(f.variable, f.ring, quotient^, budget)


@fieldwise_init
struct _HenselNode(Copyable, Movable):
    # Internal nodes hold Bezout cofactors s*left + t*right = 1 for their
    # children; leaves name their position in the lifted factor list.
    var poly: DensePoly[ZZ]
    var s: DensePoly[ZZ]
    var t: DensePoly[ZZ]
    var left: Int
    var right: Int
    var leaf: Int
    # For a right child of degree _CARRIED_DIVISION or more: the inverse of
    # its reversal modulo x^degree and inverse_modulus (zero when none). It
    # holds for every later lift of the child, which agrees modulo that.
    var inverse: DensePoly[ZZ]
    var inverse_modulus: ZZ


@fieldwise_init
struct _PendingCofactors(Copyable, Movable):
    """What a step that skipped an internal node's cofactors needs to lift
    them later: the children before the step and their corrections."""
    var g: DensePoly[ZZ]
    var h: DensePoly[ZZ]
    var x: DensePoly[ZZ]
    var r: DensePoly[ZZ]


@fieldwise_init
struct _HenselState(Movable):
    var lifted: List[DensePoly[ZZ]]
    var modulus: ZZ
    # modulus = p^exponent for the tree below; exponent zero means the
    # factor tree has not been built and lifting starts from the images.
    var exponent: Int
    var nodes: List[_HenselNode]
    # When the last step lifted only the factors, the cofactors still hold
    # modulo `pending_modulus` (zero when none are pending), and `pending`
    # keeps, by node, what lifting them the rest of the way needs.
    var pending_modulus: ZZ
    var pending: List[_PendingCofactors]


def _hensel_start(factors: List[DensePoly[ModInt]], mut budget: Budget) raises PolyError -> _HenselState:
    if len(factors) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "Hensel lifting requires modular factors")
    _allocation[DensePoly[ZZ]](len(factors), budget)
    var lifted = List[DensePoly[ZZ]](capacity=len(factors))
    for g in factors:
        budget.charge(1)
        lifted.append(_lift_image(g, budget))
    return _HenselState(lifted^, ZZ(Integer(Int(factors[0].ring.value))), 0, List[_HenselNode](), ZZ.zero(NoContext()), List[_PendingCofactors]())


def _hensel_resume(var lifted: List[DensePoly[ZZ]], var modulus: ZZ) -> _HenselState:
    """A state holding lifted factors without their tree; extending rebuilds it."""
    return _HenselState(lifted^, modulus^, 0, List[_HenselNode](), ZZ.zero(NoContext()), List[_PendingCofactors]())


comptime _CARRIED_DIVISION = 48
"""Divisor degree from which a lift step divides by a node's carried series
inverse: measured, block division with the inverse at hand beats long
division 1.3-1.8x at degree 50, 2.5-3x at 100 and 3-12x from 200."""


def _image_inverse(h: DensePoly[ModInt], mut budget: Budget) raises PolyError -> DensePoly[ModInt]:
    """The inverse of rev(h) modulo x^n for a monic h of degree n >= 1, by the
    recurrence g_i = -sum_(j=1..i) rev(h)_j g_(i-j)."""
    var n = h.degree()
    var p = h.ring.value
    _charge_size(_product(n, n) // 2 + 1, budget)
    _allocation[UInt64](_product(2, n), budget)
    var f = List[UInt64](capacity=n + 1)
    for j in range(n + 1):
        f.append(h.coefficients[n - j].value)
    var g = List[UInt64](length=n, fill=0)
    g[0] = 1
    for i in range(1, n):
        var acc = UInt128(0)
        for j in range(1, i + 1):
            acc += UInt128(f[j]) * UInt128(g[i - j])
        var r = UInt64(acc % UInt128(p))
        g[i] = 0 if r == 0 else p - r
    _allocation[ModInt](n, budget)
    var values = List[ModInt](capacity=n)
    for v in g:
        values.append(ModInt(v, h.ring))
    return DensePoly[ModInt].from_coefficients(h.variable, h.ring, values^, budget)


def _raise_inverse(h: DensePoly[ZZ], mut inverse: DensePoly[ZZ], mut inverse_modulus: ZZ, target: ZZ, mut budget: Budget) raises PolyError:
    """Newton's p-adic steps on the inverse of rev(h) modulo x^deg(h), until it
    holds modulo the target, a power of p: with rev(h) g = 1 + e modulo
    the current modulus m, g - g e holds modulo m^2."""
    if not inverse_modulus.value < target.value:
        return
    var n = h.degree()
    var flipped = _low(_reverse(h, n + 1, budget), n, budget)
    var one = _constant(h, ZZ.one(NoContext()), budget)
    while inverse_modulus.value < target.value:
        budget.charge(1)
        var next = inverse_modulus.mul(inverse_modulus, budget)
        if next.value > target.value:
            next = _copy_coefficient(target, budget)
        var e = _truncate(_low_product(flipped, inverse, n, budget).sub(one, budget), next, budget)
        inverse = _truncate(inverse.sub(_low_product(inverse, e, n, budget), budget), next, budget)
        inverse_modulus = next^


def _carried_divmod(a: DensePoly[ZZ], divisor: DensePoly[ZZ], current: DensePoly[ZZ], mut inverse: DensePoly[ZZ], mut inverse_modulus: ZZ, modulus: ZZ, mut budget: Budget) raises PolyError -> Tuple[DensePoly[ZZ], DensePoly[ZZ]]:
    """_divmod_monic by divisor, which is current modulo the modulus, with the
    carried inverse raised from current as far as the modulus."""
    if len(a.coefficients) <= divisor.degree():
        return _divmod_monic(a, divisor, modulus, budget)
    _raise_inverse(current, inverse, inverse_modulus, modulus, budget)
    return _block_divmod_with(a, divisor, _truncate(inverse, modulus, budget), modulus, budget)


def _build_tree(factors: List[DensePoly[ModInt]], mut state: _HenselState, mut budget: Budget) raises PolyError:
    # Merge least-degree nodes first (Huffman), so high-degree work sits near
    # the root and the total degree lifted below it is least; the root is the
    # last node.
    var ctx = factors[0].ring
    var x = factors[0].variable
    _allocation[_HenselNode](_product(2, len(factors)), budget)
    var nodes = List[_HenselNode]()
    var images = List[DensePoly[ModInt]]()
    var level = List[Int]()
    for i in range(len(factors)):
        budget.charge(1)
        nodes.append(_HenselNode(_lift_image(factors[i], budget), DensePoly[ZZ](x, NoContext()), DensePoly[ZZ](x, NoContext()), -1, -1, i, DensePoly[ZZ](x, NoContext()), ZZ.zero(NoContext())))
        images.append(factors[i].copy())
        level.append(i)
    # Merge the two pending nodes of least degree (ties: the earlier node).
    while len(level) > 1:
        budget.charge(1)
        var first = 0
        for k in range(1, len(level)):
            if images[level[k]].degree() < images[level[first]].degree():
                first = k
        var left = level[first]
        _ = level.pop(first)
        var second = 0
        for k in range(1, len(level)):
            if images[level[k]].degree() < images[level[second]].degree():
                second = k
        var right = level[second]
        _ = level.pop(second)
        var egcd = _extended(images[left], images[right], budget)
        if egcd[0].degree() != 0:
            raise PolyError(PolyErrorKind.DOMAIN, "Hensel lifting requires coprime modular factors")
        var product = images[left].mul(images[right], budget)
        var inverse = DensePoly[ZZ](x, NoContext())
        var inverse_modulus = ZZ.zero(NoContext())
        if images[right].degree() >= _CARRIED_DIVISION:
            inverse = _lift_image(_image_inverse(images[right], budget), budget)
            inverse_modulus = ZZ(Integer(Int(ctx.value)))
        nodes.append(_HenselNode(_lift_image(product, budget), _lift_image(egcd[1], budget), _lift_image(egcd[2], budget), left, right, -1, inverse^, inverse_modulus^))
        images.append(product^)
        level.append(len(nodes) - 1)
    for i in range(len(factors)):
        state.lifted[i] = nodes[i].poly.copy()
    state.nodes = nodes^
    state.modulus = ZZ(Integer(Int(ctx.value)))
    state.exponent = 1


comptime _FAST_DIVISION = 64
"""Divisor length from which division multiplies by a power series inverse
instead of eliminating one coefficient at a time."""
comptime _FAST_QUOTIENT = 128
"""Quotient length from which the switch pays at any modulus size: measured
at 300 to 4000 bits. Shorter quotients, down to the divisor threshold, pay
only for moduli up to _FAST_BITS."""
comptime _FAST_BITS = 1000
comptime _WIDE_DIVISION = 32
"""Divisor and quotient length from which the series inverse pays once the
modulus has _WIDE_DIVISION_BITS bits: measured on products of 256-bit
factors, 1.1-1.2x at degree 100 to 800 and neutral at degree 80; from
length 16 it lost at degree 80."""
comptime _WIDE_DIVISION_BITS = 2000


def _low_product(a: DensePoly[ZZ], b: DensePoly[ZZ], length: Int, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """a b modulo x^length."""
    if a.is_zero() or b.is_zero() or length <= 0:
        return DensePoly[ZZ](a.variable, NoContext())
    var values = _ntt_low_product(_view(a.coefficients), _view(b.coefficients), length, budget)
    if values:
        return DensePoly[ZZ].from_coefficients(a.variable, NoContext(), values.take(), budget)
    return _low(a.mul(b, budget), length, budget)


def _low(p: DensePoly[ZZ], length: Int, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """p modulo x^length."""
    var values = _slice(_view(p.coefficients), 0, min(len(p.coefficients), length), budget)
    return DensePoly[ZZ].from_coefficients(p.variable, NoContext(), values^, budget)


def _high(p: DensePoly[ZZ], lo: Int, hi: Int, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """(p modulo x^hi) div x^lo."""
    var top = min(len(p.coefficients), hi)
    var values = _slice(_view(p.coefficients), min(lo, top), top, budget)
    return DensePoly[ZZ].from_coefficients(p.variable, NoContext(), values^, budget)


def _reverse(p: DensePoly[ZZ], length: Int, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """x^(length - 1) p(1/x) for deg p < length."""
    _allocation[ZZ](length, budget)
    var values = List[ZZ](capacity=length)
    for i in range(length):
        budget.charge(1)
        var j = length - 1 - i
        values.append(_copy_coefficient(p.coefficients[j], budget) if j < len(p.coefficients) else ZZ.zero(NoContext()))
    return DensePoly[ZZ].from_coefficients(p.variable, NoContext(), values^, budget)


def _series_inverse(h: DensePoly[ZZ], length: Int, modulus: ZZ, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """g with rev(h) g = 1 modulo x^length and the modulus, for monic h.

    Newton's step from precision k to 2k: rev(h) g = 1 + x^k t, and
    g - x^k (g t) is correct to x^(2k), so only g t modulo x^k is needed."""
    var flipped = _reverse(h, h.degree() + 1, budget)
    var g = _constant(h, ZZ.one(NoContext()), budget)
    var precision = 1
    while precision < length:
        budget.charge(1)
        var next = min(_sum(precision, precision), length)
        var t = _truncate(_high(_low_product(_low(flipped, next, budget), g, next, budget), precision, next, budget), modulus, budget)
        var step = _low_product(g, t, next - precision, budget)
        _allocation[ZZ](next, budget)
        var values = List[ZZ](capacity=next)
        for i in range(next):
            budget.charge(1)
            if i < precision:
                values.append(_copy_coefficient(g.coefficients[i], budget) if i < len(g.coefficients) else ZZ.zero(NoContext()))
            else:
                var j = i - precision
                values.append(step.coefficients[j].neg(budget) if j < len(step.coefficients) else ZZ.zero(NoContext()))
        g = _truncate(DensePoly[ZZ].from_coefficients(h.variable, NoContext(), values^, budget), modulus, budget)
        precision = next
    return g^


def _divmod_monic(a: DensePoly[ZZ], h: DensePoly[ZZ], modulus: ZZ, mut budget: Budget) raises PolyError -> Tuple[DensePoly[ZZ], DensePoly[ZZ]]:
    var n = h.degree()
    var m = len(a.coefficients) - n
    var bits = modulus.value.magnitude_bit_length()
    if n >= _FAST_DIVISION and m >= _FAST_DIVISION and (m >= _FAST_QUOTIENT or bits <= _FAST_BITS):
        return _block_divmod_monic(a, h, modulus, budget)
    if n >= _WIDE_DIVISION and m >= _WIDE_DIVISION and bits >= _WIDE_DIVISION_BITS:
        return _block_divmod_monic(a, h, modulus, budget)
    return _long_divmod_monic(a, h, modulus, budget)


def _block_divmod_monic(a: DensePoly[ZZ], h: DensePoly[ZZ], modulus: ZZ, mut budget: Budget) raises PolyError -> Tuple[DensePoly[ZZ], DensePoly[ZZ]]:
    # The quotient comes from the top in blocks of at most n coefficients:
    # with the top k coefficients of the current remainder reversed, times
    # the inverse of rev(h) modulo x^k, a block clears them. Quotient and
    # remainder are unique modulo the modulus, as long division's are.
    var n = h.degree()
    var m = len(a.coefficients) - n
    return _block_divmod_with(a, h, _series_inverse(h, min(n, m), modulus, budget), modulus, budget)


def _block_divmod_with(a: DensePoly[ZZ], h: DensePoly[ZZ], inverse: DensePoly[ZZ], modulus: ZZ, mut budget: Budget) raises PolyError -> Tuple[DensePoly[ZZ], DensePoly[ZZ]]:
    """_block_divmod_monic with the inverse of rev(h) modulo x^min(n, m) and
    the modulus given, or to a higher precision in x."""
    var n = h.degree()
    var m = len(a.coefficients) - n
    var r = _slice(_view(a.coefficients), 0, len(a.coefficients), budget)
    _allocation[ZZ](m, budget)
    var q = List[ZZ](length=m, fill=ZZ.zero(NoContext()))
    var top = len(r)
    while top > n:
        budget.charge(1)
        var k = min(n, top - n)
        var shift = top - n - k
        _allocation[ZZ](k, budget)
        var head = List[ZZ](capacity=k)
        for i in range(k):
            budget.charge(1)
            head.append(_copy_coefficient(r[top - 1 - i], budget))
        var reversed_head = DensePoly[ZZ].from_coefficients(a.variable, NoContext(), head^, budget)
        var block = _truncate(_reverse(_low_product(reversed_head, _low(inverse, k, budget), k, budget), k, budget), modulus, budget)
        # The top k coefficients of block * h cancel those of r by construction.
        var product = _low_product(block, h, n, budget)
        for i in range(len(product.coefficients)):
            budget.charge(1)
            r[shift + i] = r[shift + i].sub(product.coefficients[i], budget)
        for i in range(len(block.coefficients)):
            budget.charge(1)
            q[shift + i] = _copy_coefficient(block.coefficients[i], budget)
        for i in range(shift, top - k):
            budget.charge(1)
            r[i] = ZZ(_symmetric(r[i].value, modulus.value, budget))
        top -= k
    var quotient = _truncate(DensePoly[ZZ].from_coefficients(a.variable, NoContext(), q^, budget), modulus, budget)
    var remainder = _truncate(_low(DensePoly[ZZ].from_coefficients(a.variable, NoContext(), r^, budget), n, budget), modulus, budget)
    return (quotient^, remainder^)


def _long_divmod_monic(a: DensePoly[ZZ], h: DensePoly[ZZ], modulus: ZZ, mut budget: Budget) raises PolyError -> Tuple[DensePoly[ZZ], DensePoly[ZZ]]:
    # Long division by a monic divisor over Z/modulus.
    var n = h.degree()
    var r = _slice(_view(a.coefficients), 0, len(a.coefficients), budget)
    var length = max(len(r) - n, 0)
    _allocation[ZZ](length, budget)
    var q = List[ZZ](length=length, fill=ZZ.zero(NoContext()))
    var top = len(r) - 1
    while top >= n:
        budget.charge(1)
        var c = ZZ(_symmetric(r[top].value, modulus.value, budget))
        if not c.is_zero():
            var shift = top - n
            for i in range(n):
                budget.charge(1)
                _submul(r[shift + i], c, h.coefficients[i], budget)
            q[shift] = c^
        r[top] = ZZ.zero(NoContext())
        top -= 1
    var quotient = _truncate(DensePoly[ZZ].from_coefficients(a.variable, NoContext(), q^, budget), modulus, budget)
    var remainder = _truncate(DensePoly[ZZ].from_coefficients(a.variable, NoContext(), r^, budget), modulus, budget)
    return (quotient^, remainder^)


def _quotient_exact(a: ZZ, b: ZZ, mut budget: Budget) raises PolyError -> ZZ:
    var q = _quotient(a, b, budget)
    if not q:
        raise PolyError(PolyErrorKind.DOMAIN, "the new lifting modulus is not a multiple of the old one")
    return q.take()


def _lift_tree(f: DensePoly[ZZ], mut state: _HenselState, modulus: ZZ, cofactors: Bool, mut budget: Budget) raises PolyError:
    # Without cofactors, a step lifts only the factors and keeps what the
    # cofactors' own step needs, for _finish_cofactors.
    # One quadratic step (von zur Gathen and Gerhard, Algorithm 15.10) at
    # every internal node, from the root's monic image of f downwards. The
    # new modulus divides the square of the old one.
    var inverse: Integer
    try:
        inverse = inverse_mod(f.coefficients[f.degree()].value, modulus.value)
    except e:
        raise PolyError(PolyErrorKind.DOMAIN, "leading coefficient is not invertible modulo the lifting prime")
    var target = _truncate(f.scale(ZZ(inverse^), budget), modulus, budget)
    var one = _constant(f, ZZ.one(NoContext()), budget)
    var old = _copy_coefficient(state.modulus, budget)
    var step = _quotient_exact(modulus, old, budget)
    if not cofactors:
        _allocation[_PendingCofactors](len(state.nodes), budget)
        var empty = DensePoly[ZZ](f.variable, NoContext())
        state.pending = List[_PendingCofactors](length=len(state.nodes), fill=_PendingCofactors(empty.copy(), empty.copy(), empty.copy(), empty.copy()))
        state.pending_modulus = _copy_coefficient(old, budget)
    var pending: List[Tuple[Int, DensePoly[ZZ]]] = [(len(state.nodes) - 1, target^)]
    while len(pending):
        budget.charge(1)
        var item = pending.pop()
        var v = item[0]
        if state.nodes[v].leaf >= 0:
            state.nodes[v].poly = item[1].copy()
            state.lifted[state.nodes[v].leaf] = item[1].copy()
            continue
        var left = state.nodes[v].left
        var right = state.nodes[v].right
        ref g = state.nodes[left].poly
        ref h = state.nodes[right].poly
        ref s = state.nodes[v].s
        ref t = state.nodes[v].t
        # Every correction is a multiple of the old modulus m, so with the
        # new one m * k, each is computed as m times a residue modulo k,
        # from half-size operands (von zur Gathen and Gerhard 15.10 scaled).
        var e = _truncate(_divide_scalar(item[1].sub(_hmul(g, h, budget), budget), old, budget), step, budget)
        var gs = _truncate(g, step, budget)
        var hs = _truncate(h, step, budget)
        var ss = _truncate(s, step, budget)
        var ts = _truncate(t, step, budget)
        var carried = not state.nodes[v].inverse_modulus.is_zero()
        var inverse = state.nodes[v].inverse.copy()
        var inverse_modulus = _copy_coefficient(state.nodes[v].inverse_modulus, budget)
        var dividend = _truncate(_hmul(ss, e, budget), step, budget)
        var qr = _carried_divmod(dividend, hs, h, inverse, inverse_modulus, step, budget) if carried else _divmod_monic(dividend, hs, step, budget)
        var x = _truncate(_hmul(ts, e, budget).add(_hmul(qr[0], gs, budget), budget), step, budget)
        var g2 = _truncate(g.add(x.scale(old, budget), budget), modulus, budget)
        var h2 = _truncate(h.add(qr[1].scale(old, budget), budget), modulus, budget)
        if not cofactors:
            state.pending[v] = _PendingCofactors(g.copy(), h.copy(), x^, qr[1].copy())
            state.nodes[v].poly = item[1].copy()
            state.nodes[v].inverse = inverse^
            state.nodes[v].inverse_modulus = inverse_modulus^
            pending.append((left, g2^))
            pending.append((right, h2^))
            continue
        # (s g2 + t h2 - 1) / m = (s g + t h - 1) / m + s x + t r.
        var drift = _divide_scalar(_hmul(s, g, budget).add(_hmul(t, h, budget), budget).sub(one, budget), old, budget)
        var b = _truncate(drift.add(_hmul(ss, x, budget), budget).add(_hmul(ts, qr[1], budget), budget), step, budget)
        var cofactor_dividend = _truncate(_hmul(ss, b, budget), step, budget)
        var h2s = _truncate(h2, step, budget)
        var cd = _carried_divmod(cofactor_dividend, h2s, h2, inverse, inverse_modulus, step, budget) if carried else _divmod_monic(cofactor_dividend, h2s, step, budget)
        var s2 = _truncate(s.sub(cd[1].scale(old, budget), budget), modulus, budget)
        var y = _truncate(_hmul(ts, b, budget).add(_hmul(cd[0], _truncate(g2, step, budget), budget), budget), step, budget)
        var t2 = _truncate(t.sub(y.scale(old, budget), budget), modulus, budget)
        state.nodes[v].s = s2^
        state.nodes[v].t = t2^
        state.nodes[v].poly = item[1].copy()
        state.nodes[v].inverse = inverse^
        state.nodes[v].inverse_modulus = inverse_modulus^
        pending.append((left, g2^))
        pending.append((right, h2^))


def _finish_cofactors(mut state: _HenselState, mut budget: Budget) raises PolyError:
    """Lift the cofactors a factor-only step left behind, by the same
    expressions _lift_tree would have used, to the state's modulus."""
    if state.pending_modulus.is_zero():
        return
    var old = _copy_coefficient(state.pending_modulus, budget)
    ref modulus = state.modulus
    var step = _quotient_exact(modulus, old, budget)
    for v in range(len(state.nodes)):
        budget.charge(1)
        if state.nodes[v].leaf >= 0:
            continue
        ref saved = state.pending[v]
        ref g = saved.g
        ref h = saved.h
        ref g2 = state.nodes[state.nodes[v].left].poly
        ref h2 = state.nodes[state.nodes[v].right].poly
        ref s = state.nodes[v].s
        ref t = state.nodes[v].t
        var one = _constant(g2, ZZ.one(NoContext()), budget)
        var ss = _truncate(s, step, budget)
        var ts = _truncate(t, step, budget)
        # (s g2 + t h2 - 1) / m = (s g + t h - 1) / m + s x + t r.
        var drift = _divide_scalar(_hmul(s, g, budget).add(_hmul(t, h, budget), budget).sub(one, budget), old, budget)
        var b = _truncate(drift.add(_hmul(ss, saved.x, budget), budget).add(_hmul(ts, saved.r, budget), budget), step, budget)
        var carried = not state.nodes[v].inverse_modulus.is_zero()
        var inverse = state.nodes[v].inverse.copy()
        var inverse_modulus = _copy_coefficient(state.nodes[v].inverse_modulus, budget)
        var cofactor_dividend = _truncate(_hmul(ss, b, budget), step, budget)
        var h2s = _truncate(h2, step, budget)
        var cd = _carried_divmod(cofactor_dividend, h2s, h2, inverse, inverse_modulus, step, budget) if carried else _divmod_monic(cofactor_dividend, h2s, step, budget)
        var s2 = _truncate(s.sub(cd[1].scale(old, budget), budget), modulus, budget)
        var y = _truncate(_hmul(ts, b, budget).add(_hmul(cd[0], _truncate(g2, step, budget), budget), budget), step, budget)
        var t2 = _truncate(t.sub(y.scale(old, budget), budget), modulus, budget)
        state.nodes[v].s = s2^
        state.nodes[v].t = t2^
        state.nodes[v].inverse = inverse^
        state.nodes[v].inverse_modulus = inverse_modulus^
    state.pending = List[_PendingCofactors]()
    state.pending_modulus = ZZ.zero(NoContext())


def _hensel_extend(f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], bound: ZZ, mut state: _HenselState, mut budget: Budget) raises PolyError:
    # Lift the monic factors until the modulus exceeds the bound, to the
    # least such power of p, by quadratic steps over a binary factor tree.
    # The leading coefficient of f is carried separately, so every lifted
    # factor stays monic with its modular degree.
    var p = ZZ(Integer(Int(factors[0].ring.value)))
    _work(state.modulus.size(), bound.size(), False, budget)
    if state.modulus.value > bound.value:
        return
    if state.exponent == 0:
        # Without a tree, a lift below 2^254 runs on fixed-width residues,
        # whose lifted monic factors are the same. It builds no tree here,
        # so a later extension builds one and lifts from the images again.
        var fast = _word_hensel(f, factors, bound, budget)
        if not fast:
            fast = _double_hensel(f, factors, bound, budget)
        if not fast:
            fast = _wide_hensel(f, factors, bound, budget)
        if fast:
            var lifted = fast.take()
            state.lifted = lifted[0].copy()
            state.modulus = _copy_coefficient(lifted[1], budget)
            return
        _build_tree(factors, state, budget)
    _finish_cofactors(state, budget)
    var exponent = state.exponent
    var power = _copy_coefficient(state.modulus, budget)
    while not power.value > bound.value:
        budget.charge(1)
        _work(power.size(), bound.size(), False, budget)
        power = power.mul(p, budget)
        exponent = _sum(exponent, 1)
    var schedule = List[Int]()
    var e = exponent
    while e > state.exponent:
        budget.charge(1)
        _allocation[Int](1, budget)
        schedule.append(e)
        e = (e + 1) // 2
    var i = len(schedule) - 1
    while i >= 0:
        budget.charge(1)
        var modulus = _pow(p, UInt64(schedule[i]), budget)
        # The last step leaves the cofactors to a later extension, if any.
        _lift_tree(f, state, modulus, i > 0, budget)
        state.modulus = modulus^
        state.exponent = schedule[i]
        i -= 1


# ------------------------------------------------ word-size moduli

comptime _WORD_BITS = 62
"""Moduli below 2^62 keep products of residues below 2^124 and sums of
two residues below 2^63."""


def _wtrim(mut a: List[UInt64]):
    while len(a) > 0 and a[len(a) - 1] == 0:
        _ = a.pop()


def _wadd(a: List[UInt64], b: List[UInt64], m: UInt64, mut budget: Budget) raises PolyError -> List[UInt64]:
    var n = max(len(a), len(b))
    _allocation[UInt64](n, budget)
    var c = List[UInt64](capacity=n)
    for i in range(n):
        var x = a[i] if i < len(a) else 0
        var y = b[i] if i < len(b) else 0
        var z = x + y
        c.append(z - m if z >= m else z)
    _wtrim(c)
    return c^


def _wsub(a: List[UInt64], b: List[UInt64], m: UInt64, mut budget: Budget) raises PolyError -> List[UInt64]:
    var n = max(len(a), len(b))
    _allocation[UInt64](n, budget)
    var c = List[UInt64](capacity=n)
    for i in range(n):
        var x = a[i] if i < len(a) else 0
        var y = b[i] if i < len(b) else 0
        c.append(x - y if x >= y else x + (m - y))
    _wtrim(c)
    return c^


def _wmul(a: List[UInt64], b: List[UInt64], m: UInt64, mut budget: Budget) raises PolyError -> List[UInt64]:
    if len(a) == 0 or len(b) == 0:
        return List[UInt64]()
    var n = len(a) + len(b) - 1
    _allocation[UInt64](n, budget)
    _charge_size(_product(len(a), len(b)), budget)
    var c = List[UInt64](length=n, fill=0)
    var mm = UInt128(m)
    for i in range(len(a)):
        if a[i] == 0:
            continue
        var x = UInt128(a[i])
        for j in range(len(b)):
            var z = UInt64((x * UInt128(b[j])) % mm) + c[i + j]
            c[i + j] = z - m if z >= m else z
    _wtrim(c)
    return c^


def _wdivmod(a: List[UInt64], h: List[UInt64], m: UInt64, mut budget: Budget) raises PolyError -> Tuple[List[UInt64], List[UInt64]]:
    """Quotient and remainder by a monic h modulo m."""
    var n = len(h) - 1
    var r = a.copy()
    if len(r) <= n:
        return (List[UInt64](), r^)
    var q = List[UInt64](length=len(r) - n, fill=0)
    _charge_size(_product(len(q), len(h)), budget)
    var mm = UInt128(m)
    var top = len(r) - 1
    while top >= n:
        var c = r[top]
        if c != 0:
            var shift = top - n
            q[shift] = c
            var x = UInt128(c)
            for i in range(n):
                var y = UInt64((x * UInt128(h[i])) % mm)
                var z = r[shift + i]
                r[shift + i] = z - y if z >= y else z + (m - y)
        r[top] = 0
        top -= 1
    _wtrim(r)
    _wtrim(q)
    return (q^, r^)


def _word(f: DensePoly[ZZ], m: UInt64, mut budget: Budget) raises PolyError -> List[UInt64]:
    var modulus = Integer(Int(m))
    _allocation[UInt64](len(f.coefficients), budget)
    var result = List[UInt64](capacity=len(f.coefficients))
    for c in f.coefficients:
        budget.charge(1)
        try:
            result.append(UInt64(Int(div_rem_euclid(c.value, modulus)[1])))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    _wtrim(result)
    return result^


def _symmetric_word(a: List[UInt64], m: UInt64, variable: Int64, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(a), budget)
    var values = List[ZZ](capacity=len(a))
    for x in a:
        budget.charge(1)
        if x > m - x:
            values.append(ZZ(Integer(-Int(m - x))))
        else:
            values.append(ZZ(Integer(Int(x))))
    return DensePoly[ZZ].from_coefficients(variable, NoContext(), values^, budget)


def _word_hensel(f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], bound: ZZ, mut budget: Budget) raises PolyError -> Optional[Tuple[List[DensePoly[ZZ]], ZZ]]:
    """_hensel on machine words when the final modulus is below 2^62: the
    same quadratic tree steps on residues in [0, m), whose lifted monic
    factors are unique modulo m and so equal the big-integer lift's."""
    var p = factors[0].ring.value
    var power = Integer(Int(p))
    var exponent = 1
    while not power > bound.value:
        budget.charge(1)
        try:
            power = power * Integer(Int(p))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        exponent += 1
        if power.magnitude_bit_length() >= _WORD_BITS:
            return None
    # No step at all: the big-integer path returns the images as they are.
    if exponent == 1 or power.magnitude_bit_length() >= _WORD_BITS:
        return None
    var state = _hensel_start(factors, budget)
    _build_tree(factors, state, budget)
    var count = len(state.nodes)
    var polys = List[List[UInt64]]()
    var ss = List[List[UInt64]]()
    var ts = List[List[UInt64]]()
    for node in state.nodes:
        budget.charge(1)
        polys.append(_word(node.poly, p, budget))
        ss.append(_word(node.s, p, budget))
        ts.append(_word(node.t, p, budget))
    var schedule = List[Int]()
    var e = exponent
    while e > 1:
        budget.charge(1)
        schedule.append(e)
        e = (e + 1) // 2
    var m = p
    var i = len(schedule) - 1
    while i >= 0:
        budget.charge(1)
        m = 1
        for _ in range(schedule[i]):
            m *= p
        var inverse: Integer
        try:
            inverse = inverse_mod(f.coefficients[f.degree()].value, Integer(Int(m)))
        except err:
            raise PolyError(PolyErrorKind.DOMAIN, "leading coefficient is not invertible modulo the lifting prime")
        var monic = _word(f.scale(ZZ(inverse^), budget), m, budget)
        var one = List[UInt64]()
        one.append(1)
        var pending: List[Tuple[Int, List[UInt64]]] = [(count - 1, monic^)]
        while len(pending):
            budget.charge(1)
            var item = pending.pop()
            var v = item[0]
            if state.nodes[v].leaf >= 0:
                polys[v] = item[1].copy()
                continue
            var left = state.nodes[v].left
            var right = state.nodes[v].right
            var g = polys[left].copy()
            var h = polys[right].copy()
            var err = _wsub(item[1], _wmul(g, h, m, budget), m, budget)
            var qr = _wdivmod(_wmul(ss[v], err, m, budget), h, m, budget)
            var g2 = _wadd(_wadd(g, _wmul(ts[v], err, m, budget), m, budget), _wmul(qr[0], g, m, budget), m, budget)
            var h2 = _wadd(h, qr[1], m, budget)
            var b = _wsub(_wadd(_wmul(ss[v], g2, m, budget), _wmul(ts[v], h2, m, budget), m, budget), one, m, budget)
            var cd = _wdivmod(_wmul(ss[v], b, m, budget), h2, m, budget)
            var s2 = _wsub(ss[v], cd[1], m, budget)
            var t2 = _wsub(_wsub(ts[v], _wmul(ts[v], b, m, budget), m, budget), _wmul(cd[0], g2, m, budget), m, budget)
            ss[v] = s2^
            ts[v] = t2^
            polys[v] = item[1].copy()
            pending.append((left, g2^))
            pending.append((right, h2^))
        i -= 1
    var lifted = List[DensePoly[ZZ]](capacity=len(factors))
    for _ in range(len(factors)):
        lifted.append(DensePoly[ZZ](f.variable, NoContext()))
    for v in range(count):
        budget.charge(1)
        if state.nodes[v].leaf >= 0:
            lifted[state.nodes[v].leaf] = _symmetric_word(polys[v], m, f.variable, budget)
    return (lifted^, ZZ(Integer(Int(m))))


# ------------------------------------------------ two-word moduli

comptime _DOUBLE_BITS = 124
"""Odd moduli below 2^124: Montgomery residues with R = 2^128, products of
residues below 2^248."""
comptime _LOW = UInt128(UInt64.MAX)


@fieldwise_init
struct _Montgomery(Copyable, Movable):
    var m: UInt128
    var minv: UInt128
    """-m^-1 modulo 2^128."""
    var r2: UInt128
    """R^2 modulo m."""
    var one: UInt128
    """R modulo m."""


def _wide(a: UInt128, b: UInt128) -> Tuple[UInt128, UInt128]:
    """The 256-bit product as (high, low) 128-bit halves."""
    var a0 = a & _LOW
    var a1 = a >> 64
    var b0 = b & _LOW
    var b1 = b >> 64
    var p00 = a0 * b0
    var p01 = a0 * b1
    var p10 = a1 * b0
    var p11 = a1 * b1
    var mid = (p00 >> 64) + (p01 & _LOW) + (p10 & _LOW)
    var low = (p00 & _LOW) | (mid << 64)
    var high = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64)
    return (high, low)


def _redc(high: UInt128, low: UInt128, ctx: _Montgomery) -> UInt128:
    var u = low * ctx.minv
    var um = _wide(u, ctx.m)
    # low + um.low is 0 or 2^128 by the choice of u.
    var carry = UInt128(1) if low != 0 else UInt128(0)
    var t = high + um[0] + carry
    return t - ctx.m if t >= ctx.m else t


def _mmul(a: UInt128, b: UInt128, ctx: _Montgomery) -> UInt128:
    var w = _wide(a, b)
    return _redc(w[0], w[1], ctx)


def _montgomery(m: UInt128, r2: UInt128) -> _Montgomery:
    var inverse = m
    for _ in range(7):
        inverse = inverse * (UInt128(2) - m * inverse)
    var ctx = _Montgomery(m, UInt128(0) - inverse, r2, 0)
    ctx.one = _redc(0, r2, ctx)
    return ctx^


def _dtrim(mut a: List[UInt128]):
    while len(a) > 0 and a[len(a) - 1] == 0:
        _ = a.pop()


def _dadd(a: List[UInt128], b: List[UInt128], ctx: _Montgomery, mut budget: Budget) raises PolyError -> List[UInt128]:
    var n = max(len(a), len(b))
    _allocation[UInt128](n, budget)
    var c = List[UInt128](capacity=n)
    for i in range(n):
        var x = a[i] if i < len(a) else 0
        var y = b[i] if i < len(b) else 0
        var z = x + y
        c.append(z - ctx.m if z >= ctx.m else z)
    _dtrim(c)
    return c^


def _dsub(a: List[UInt128], b: List[UInt128], ctx: _Montgomery, mut budget: Budget) raises PolyError -> List[UInt128]:
    var n = max(len(a), len(b))
    _allocation[UInt128](n, budget)
    var c = List[UInt128](capacity=n)
    for i in range(n):
        var x = a[i] if i < len(a) else 0
        var y = b[i] if i < len(b) else 0
        c.append(x - y if x >= y else x + (ctx.m - y))
    _dtrim(c)
    return c^


def _dmul(a: List[UInt128], b: List[UInt128], ctx: _Montgomery, mut budget: Budget) raises PolyError -> List[UInt128]:
    if len(a) == 0 or len(b) == 0:
        return List[UInt128]()
    var n = len(a) + len(b) - 1
    _allocation[UInt128](n, budget)
    _charge_size(_product(len(a), len(b)), budget)
    var c = List[UInt128](length=n, fill=0)
    for i in range(len(a)):
        if a[i] == 0:
            continue
        for j in range(len(b)):
            var z = _mmul(a[i], b[j], ctx) + c[i + j]
            c[i + j] = z - ctx.m if z >= ctx.m else z
    _dtrim(c)
    return c^


def _ddivmod(a: List[UInt128], h: List[UInt128], ctx: _Montgomery, mut budget: Budget) raises PolyError -> Tuple[List[UInt128], List[UInt128]]:
    """Quotient and remainder by a monic h (leading residue R mod m)."""
    var n = len(h) - 1
    var r = a.copy()
    if len(r) <= n:
        return (List[UInt128](), r^)
    var q = List[UInt128](length=len(r) - n, fill=0)
    _charge_size(_product(len(q), len(h)), budget)
    var top = len(r) - 1
    while top >= n:
        var c = r[top]
        if c != 0:
            var shift = top - n
            q[shift] = c
            for i in range(n):
                var y = _mmul(c, h[i], ctx)
                var z = r[shift + i]
                r[shift + i] = z - y if z >= y else z + (ctx.m - y)
        r[top] = 0
        top -= 1
    _dtrim(r)
    _dtrim(q)
    return (q^, r^)


def _to_u128(a: Integer) raises PolyError -> UInt128:
    """A nonnegative integer below 2^128, by 32-bit chunks."""
    var result = UInt128(0)
    try:
        for k in range(4):
            var chunk = Int((a >> (32 * k)) & Integer(0xFFFFFFFF))
            result |= UInt128(UInt64(chunk)) << UInt128(32 * k)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return result


def _from_u128(a: UInt128) raises PolyError -> Integer:
    var result = Integer(0)
    try:
        var k = 3
        while k >= 0:
            result = (result << 32) + Integer(Int(UInt64((a >> UInt128(32 * k)) & UInt128(0xFFFFFFFF))))
            k -= 1
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return result^


def _double(f: DensePoly[ZZ], ctx: _Montgomery, modulus: Integer, mut budget: Budget) raises PolyError -> List[UInt128]:
    _allocation[UInt128](len(f.coefficients), budget)
    var result = List[UInt128](capacity=len(f.coefficients))
    for c in f.coefficients:
        budget.charge(1)
        var r: Integer
        try:
            r = div_rem_euclid(c.value, modulus)[1].copy()
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        result.append(_mmul(_to_u128(r), ctx.r2, ctx))
    _dtrim(result)
    return result^


def _symmetric_double(a: List[UInt128], ctx: _Montgomery, variable: Int64, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(a), budget)
    var values = List[ZZ](capacity=len(a))
    for x in a:
        budget.charge(1)
        var v = _redc(0, x, ctx)
        if v > ctx.m - v:
            values.append(ZZ(_from_u128(ctx.m - v)).neg(budget))
        else:
            values.append(ZZ(_from_u128(v)))
    return DensePoly[ZZ].from_coefficients(variable, NoContext(), values^, budget)


def _double_hensel(f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], bound: ZZ, mut budget: Budget) raises PolyError -> Optional[Tuple[List[DensePoly[ZZ]], ZZ]]:
    """_word_hensel for moduli in [2^62, 2^124), on Montgomery residues."""
    var p = factors[0].ring.value
    var power = Integer(Int(p))
    var exponent = 1
    while not power > bound.value:
        budget.charge(1)
        try:
            power = power * Integer(Int(p))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        exponent += 1
        if power.magnitude_bit_length() >= _DOUBLE_BITS:
            return None
    if exponent == 1 or power.magnitude_bit_length() >= _DOUBLE_BITS:
        return None
    var state = _hensel_start(factors, budget)
    _build_tree(factors, state, budget)
    var count = len(state.nodes)
    var schedule = List[Int]()
    var e = exponent
    while e > 1:
        budget.charge(1)
        schedule.append(e)
        e = (e + 1) // 2
    # Residues mod p are valid representatives mod every power of p, so the
    # tree enters each modulus by conversion from its previous values.
    var polys = List[DensePoly[ZZ]]()
    var ss = List[DensePoly[ZZ]]()
    var ts = List[DensePoly[ZZ]]()
    for node in state.nodes:
        polys.append(node.poly.copy())
        ss.append(node.s.copy())
        ts.append(node.t.copy())
    var modulus = Integer(1)
    var i = len(schedule) - 1
    while i >= 0:
        budget.charge(1)
        modulus = _pow(ZZ(Integer(Int(p))), UInt64(schedule[i]), budget).value.copy()
        var r2: Integer
        try:
            r2 = div_rem_euclid(Integer(1) << 256, modulus)[1].copy()
        except err:
            raise PolyError(PolyErrorKind.CAPACITY, String(err))
        var ctx = _montgomery(_to_u128(modulus), _to_u128(r2))
        var wp = List[List[UInt128]]()
        var ws = List[List[UInt128]]()
        var wt = List[List[UInt128]]()
        for v in range(count):
            wp.append(_double(polys[v], ctx, modulus, budget))
            ws.append(_double(ss[v], ctx, modulus, budget))
            wt.append(_double(ts[v], ctx, modulus, budget))
        var inverse: Integer
        try:
            inverse = inverse_mod(f.coefficients[f.degree()].value, modulus)
        except err:
            raise PolyError(PolyErrorKind.DOMAIN, "leading coefficient is not invertible modulo the lifting prime")
        var monic = _double(f.scale(ZZ(inverse^), budget), ctx, modulus, budget)
        var one = List[UInt128]()
        one.append(ctx.one)
        var pending: List[Tuple[Int, List[UInt128]]] = [(count - 1, monic^)]
        while len(pending):
            budget.charge(1)
            var item = pending.pop()
            var v = item[0]
            if state.nodes[v].leaf >= 0:
                wp[v] = item[1].copy()
                continue
            var left = state.nodes[v].left
            var right = state.nodes[v].right
            var g = wp[left].copy()
            var h = wp[right].copy()
            var err = _dsub(item[1], _dmul(g, h, ctx, budget), ctx, budget)
            var qr = _ddivmod(_dmul(ws[v], err, ctx, budget), h, ctx, budget)
            var g2 = _dadd(_dadd(g, _dmul(wt[v], err, ctx, budget), ctx, budget), _dmul(qr[0], g, ctx, budget), ctx, budget)
            var h2 = _dadd(h, qr[1], ctx, budget)
            var b = _dsub(_dadd(_dmul(ws[v], g2, ctx, budget), _dmul(wt[v], h2, ctx, budget), ctx, budget), one, ctx, budget)
            var cd = _ddivmod(_dmul(ws[v], b, ctx, budget), h2, ctx, budget)
            ws[v] = _dsub(ws[v], cd[1], ctx, budget)
            wt[v] = _dsub(_dsub(wt[v], _dmul(wt[v], b, ctx, budget), ctx, budget), _dmul(cd[0], g2, ctx, budget), ctx, budget)
            wp[v] = item[1].copy()
            pending.append((left, g2^))
            pending.append((right, h2^))
        for v in range(count):
            polys[v] = _symmetric_double(wp[v], ctx, f.variable, budget)
            ss[v] = _symmetric_double(ws[v], ctx, f.variable, budget)
            ts[v] = _symmetric_double(wt[v], ctx, f.variable, budget)
        i -= 1
    var lifted = List[DensePoly[ZZ]](capacity=len(factors))
    for _ in range(len(factors)):
        lifted.append(DensePoly[ZZ](f.variable, NoContext()))
    for v in range(count):
        budget.charge(1)
        if state.nodes[v].leaf >= 0:
            lifted[state.nodes[v].leaf] = polys[v].copy()
    return (lifted^, ZZ(modulus^))


comptime _WIDE_BITS = 254
"""Moduli below 2^254 lift on three or four 64-bit limbs."""
comptime _WIDE_DEGREE = 128
"""Degrees up to which quadratic limb products beat the transform products of the integer tree."""


def _wide_hensel(f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], bound: ZZ, mut budget: Budget) raises PolyError -> Optional[Tuple[List[DensePoly[ZZ]], ZZ]]:
    """_double_hensel for moduli in [2^124, 2^254), on three or four limbs."""
    var p = factors[0].ring.value
    if p % 2 == 0 or f.degree() > _WIDE_DEGREE:
        return None
    var power = Integer(Int(p))
    var exponent = 1
    while not power > bound.value:
        budget.charge(1)
        try:
            power = power * Integer(Int(p))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        exponent += 1
        if power.magnitude_bit_length() >= _WIDE_BITS:
            return None
    if exponent == 1 or power.magnitude_bit_length() < _DOUBLE_BITS or power.magnitude_bit_length() >= _WIDE_BITS:
        return None
    if power.magnitude_bit_length() <= 64 * 3 - 2:
        return _wide_lift[3](f, factors, p, exponent, budget)
    return _wide_lift[4](f, factors, p, exponent, budget)


def _wide_lift[N: Int](f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], p: UInt64, exponent: Int, mut budget: Budget) raises PolyError -> Tuple[List[DensePoly[ZZ]], ZZ]:
    var state = _hensel_start(factors, budget)
    _build_tree(factors, state, budget)
    var count = len(state.nodes)
    var schedule = List[Int]()
    var e = exponent
    while e > 1:
        budget.charge(1)
        schedule.append(e)
        e = (e + 1) // 2
    # Residues mod p are valid representatives mod every power of p, so the
    # tree enters each modulus by conversion from its previous values.
    var polys = List[DensePoly[ZZ]]()
    var ss = List[DensePoly[ZZ]]()
    var ts = List[DensePoly[ZZ]]()
    for node in state.nodes:
        polys.append(node.poly.copy())
        ss.append(node.s.copy())
        ts.append(node.t.copy())
    var modulus = Integer(1)
    var i = len(schedule) - 1
    while i >= 0:
        budget.charge(1)
        modulus = _pow(ZZ(Integer(Int(p))), UInt64(schedule[i]), budget).value.copy()
        var ctx = _wide_context[N](modulus)
        var wp = List[List[_Limbs]]()
        var ws = List[List[_Limbs]]()
        var wt = List[List[_Limbs]]()
        for v in range(count):
            wp.append(_limbs_of[N](polys[v], ctx, modulus, budget))
            ws.append(_limbs_of[N](ss[v], ctx, modulus, budget))
            wt.append(_limbs_of[N](ts[v], ctx, modulus, budget))
        var inverse: Integer
        try:
            inverse = inverse_mod(f.coefficients[f.degree()].value, modulus)
        except err:
            raise PolyError(PolyErrorKind.DOMAIN, "leading coefficient is not invertible modulo the lifting prime")
        var monic = _limbs_of[N](f.scale(ZZ(inverse^), budget), ctx, modulus, budget)
        var one = List[_Limbs]()
        one.append(ctx.one)
        var pending: List[Tuple[Int, List[_Limbs]]] = [(count - 1, monic^)]
        while len(pending):
            budget.charge(1)
            var item = pending.pop()
            var v = item[0]
            if state.nodes[v].leaf >= 0:
                wp[v] = item[1].copy()
                continue
            var left = state.nodes[v].left
            var right = state.nodes[v].right
            var g = wp[left].copy()
            var h = wp[right].copy()
            var err = _lsub[N](item[1], _lmul[N](g, h, ctx, budget), ctx, budget)
            var qr = _ldivmod[N](_lmul[N](ws[v], err, ctx, budget), h, ctx, budget)
            var g2 = _ladd[N](_ladd[N](g, _lmul[N](wt[v], err, ctx, budget), ctx, budget), _lmul[N](qr[0], g, ctx, budget), ctx, budget)
            var h2 = _ladd[N](h, qr[1], ctx, budget)
            var b = _lsub[N](_ladd[N](_lmul[N](ws[v], g2, ctx, budget), _lmul[N](wt[v], h2, ctx, budget), ctx, budget), one, ctx, budget)
            var cd = _ldivmod[N](_lmul[N](ws[v], b, ctx, budget), h2, ctx, budget)
            ws[v] = _lsub[N](ws[v], cd[1], ctx, budget)
            wt[v] = _lsub[N](_lsub[N](wt[v], _lmul[N](wt[v], b, ctx, budget), ctx, budget), _lmul[N](cd[0], g2, ctx, budget), ctx, budget)
            wp[v] = item[1].copy()
            pending.append((left, g2^))
            pending.append((right, h2^))
        for v in range(count):
            polys[v] = _symmetric_limbs[N](wp[v], ctx, f.variable, budget)
            ss[v] = _symmetric_limbs[N](ws[v], ctx, f.variable, budget)
            ts[v] = _symmetric_limbs[N](wt[v], ctx, f.variable, budget)
        i -= 1
    var lifted = List[DensePoly[ZZ]](capacity=len(factors))
    for _ in range(len(factors)):
        lifted.append(DensePoly[ZZ](f.variable, NoContext()))
    for v in range(count):
        budget.charge(1)
        if state.nodes[v].leaf >= 0:
            lifted[state.nodes[v].leaf] = polys[v].copy()
    return (lifted^, ZZ(modulus^))


def _hensel(f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], bound: ZZ, mut budget: Budget) raises PolyError -> Tuple[List[DensePoly[ZZ]], ZZ]:
    var word = _word_hensel(f, factors, bound, budget)
    if word:
        return word.take()
    var double = _double_hensel(f, factors, bound, budget)
    if double:
        return double.take()
    var wide = _wide_hensel(f, factors, bound, budget)
    if wide:
        return wide.take()
    var state = _hensel_start(factors, budget)
    _hensel_extend(f, factors, bound, state, budget)
    var lifted = List[DensePoly[ZZ]]()
    var modulus = ZZ.zero(NoContext())
    swap(lifted, state.lifted)
    swap(modulus, state.modulus)
    return (lifted^, modulus^)
