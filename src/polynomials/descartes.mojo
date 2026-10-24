"""Real root isolation by Descartes' rule of signs, with checkable certificates.

Every real root of the squarefree part g of f lies in (-2^k, 2^k). Dyadic
bisection in the manner of Vincent, Collins and Akritas splits each half of
that interval until each piece is resolved: Descartes' bound on it is zero
or one, both exact counts, or its v roots show as v sign changes. Newton
steps shorten the descent to clusters, and a close real pair is split at a
point between its roots. The certificate holds the pieces, the rational
roots met at split points, and a bracket for every other root. Checking it
recomputes one bound per piece and two signs per bracket; f = g gcd(f, f')
shows that f and g have the same roots.
"""
from apn_mojo import Integer, Rational
from .coefficients import ZZ, QQ, NoContext, CoefficientRing
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from .numeric_bounds import NumericBox
from .roots import RootInterval, RootIsolation, _rational_root
from .numeric_options import NumericOptions, NumericOutcome, _accurate
from .numeric_evaluate import _exact_dense_boxes
from ._images import _integer_dense
from ._integer_ops import _primitive
from ._univariate_modular import _univariate_gcd
from ._hensel import _trial_divide
from .certificates import _fujiwara_exponent
from ._root_ops import _sign_at, _midpoint, _compare, _order, _twos, _bits, _power_of_two, _error_bits, _fixed, _exact_sign, _sign, _add_shifted, _Words
from ._real_refine import _RootOracle, _Enclosure, _enclosure, _check_sign, _qir_step, _side, _ratio_bits, _wider, _derivative
from ._numeric_roots import _rational_polynomial, _require_success
from .root_stream import RootStreamReceipt
from .division import divmod


@fieldwise_init
struct DescartesLeaf(Copyable, Movable):
    """An open interval of a Descartes certificate with its root bound."""
    var low: QQ
    """The lower endpoint."""
    var high: QQ
    """The upper endpoint, above low."""
    var variations: Int
    """Descartes' bound for the roots of g in (low, high): at least their number, with the same parity; zero and one are exact counts."""
    var resolved: Bool
    """Whether `variations` brackets follow for this leaf; always true for bounds below two.

    Disjoint brackets with sign changes each hold an odd number of roots, so
    v of them under a bound of v hold exactly one root each.
    """


@fieldwise_init
struct DescartesCertificate(Copyable, Movable):
    """Leaves partitioning (low, high), the roots of g at their shared ends, and a bracket for every other root."""
    var low: QQ
    """The lower end of the partitioned interval, not a root of g."""
    var high: QQ
    """The upper end of the partitioned interval, not a root of g."""
    var leaves: List[DescartesLeaf]
    """Increasing leaves; each leaf's high is the next leaf's low."""
    var exact_roots: List[QQ]
    """Increasing: exactly the shared leaf endpoints where g is zero."""
    var brackets: List[RootInterval]
    """As many per resolved leaf as its bound, increasing and disjoint inside it: an open interval where g changes sign, or a root low == high of g."""


@fieldwise_init
struct DescartesRootResult(Copyable, Movable):
    """The distinct real roots of an exact polynomial with a Descartes certificate and outcome flags."""
    var outcome: NumericOutcome
    """SUCCESS only when every requested guarantee holds."""
    var complete: Bool
    """Every leaf is resolved, so every real root is counted."""
    var accurate: Bool
    """Every bracket meets the requested widths."""
    var isolated: Bool
    """Every counted root has its own bracket or exact value; equals complete."""
    var squarefree: DensePoly[ZZ]
    """The primitive positive-leading g = f / gcd(f, f'): the same roots, each simple."""
    var cofactor: DensePoly[ZZ]
    """The cofactor gcd(f, f'), so that the primitive part of f is g times it."""
    var bound_exponent: Int
    """Every root lies in (-2^k, 2^k), the certificate's interval."""
    var certificate: DescartesCertificate
    """The leaves, exact roots and brackets for g."""


@fieldwise_init
struct _Node(Movable):
    # The interval (low, low + 2^-width) for low = index / 2^depth, depth >=
    # width, in the coordinates of one scaled side. A search node holds its
    # polynomial; a root node is the exact root low; an empty node is a
    # piece without roots.
    var poly: List[Integer]
    var index: Integer
    var depth: Int
    var width: Int
    var kind: Int
    var newton: Int
    """log2 of the Newton grid: 2^newton parts of the interval."""
    var parent: Int
    """The parent's bound, so that a node whose bound did not drop can try a Newton step."""
    var paired: Int
    """The width at which a search for a point between a close pair last failed, or -1."""
    var seen: Int
    """Sign changes found by the last failed grid search above this node, or -1."""
    var seen_width: Int
    var seen_bound: Int
    """The exact bound where that search ran, which bounds this node's."""
    var end: Integer
    """An empty node's upper end, end / 2^depth, or a known node's bound."""


comptime _SEARCH = 0
comptime _ROOT = 1
comptime _EMPTY = 2
comptime _KNOWN = 3
"""A search node whose bound is end."""
comptime _MAX_NEWTON = 1 << 16


def _error(e: Error) -> PolyError:
    return PolyError(PolyErrorKind.CAPACITY, String(e))


def _shift_one(mut c: List[Integer]) raises:
    # c(x) <- c(x + 1), ascending.
    var n = len(c) - 1
    var w = _shift_width(c)
    var data = _load(c, w, False)
    _ = _taylor_words(data, n, w, 0)
    c = _store(data, n + 1, w)


# Taylor shifts by one on fixed-width two's complement limbs: every
# coefficient takes w 64-bit limbs, lowest first, so the n^2 / 2 additions
# are carry chains over plain words rather than allocating integer sums.


def _shift_width(c: List[Integer]) -> Int:
    # Coefficients of c(x + 1) stay below 2^(n + 1) max |c_i| in magnitude.
    var top = 0
    for i in range(len(c)):
        top = max(top, c[i].magnitude_bit_length())
    return (top + len(c) + 2) // 64 + 1


def _load(c: List[Integer], w: Int, reverse: Bool) -> List[UInt64]:
    var count = len(c)
    var data = List[UInt64](length=count * w, fill=0)
    for i in range(count):
        ref x = c[count - 1 - i if reverse else i]
        var base = i * w
        var words = x._word_count()
        for j in range(words):
            data[base + j // 2] |= UInt64(x._word(j)) << UInt64(32 * (j % 2))
        if x.sign() < 0:
            var carry = UInt64(1)
            for t in range(w):
                var v = ~data[base + t] + carry
                carry = UInt64(1) if carry == 1 and v == 0 else UInt64(0)
                data[base + t] = v
    return data^


def _store(data: List[UInt64], count: Int, w: Int) -> List[Integer]:
    var out = List[Integer](capacity=count)
    for i in range(count):
        var base = i * w
        var negative = (data[base + w - 1] >> 63) != 0
        var words = List[UInt32](length=2 * w, fill=0)
        var carry = UInt64(1)
        for t in range(w):
            var v = data[base + t]
            if negative:
                v = ~v + carry
                carry = UInt64(1) if carry == 1 and v == 0 else UInt64(0)
            words[2 * t] = UInt32(v & 0xFFFFFFFF)
            words[2 * t + 1] = UInt32(v >> 32)
        out.append(Integer._from_words(words^, negative))
    return out^


@always_inline
def _limb_sign(p: Pointer[mut=True, UInt64, _], base: Int, w: Int) -> Int:
    if (p.unsafe_offset(base + w - 1).unsafe_load() >> 63) != 0:
        return -1
    for t in range(w):
        if p.unsafe_offset(base + t).unsafe_load() != 0:
            return 1
    return 0


def _taylor_words(mut data: List[UInt64], n: Int, w: Int, stop: Int) -> Int:
    """Shift in place by one, round by round; with stop > 0, end once the fixed prefix shows stop variations.

    Returns the variations counted, capped at stop, when stop is positive.
    """
    var p = data.unsafe_ptr()
    var top = _limb_sign(p, n * w, w)
    var count = 0
    var last = 0
    for i in range(n):
        if w == 1:
            for j in range(n - 1, i - 1, -1):
                p.unsafe_offset(j).unsafe_store(p.unsafe_offset(j).unsafe_load() + p.unsafe_offset(j + 1).unsafe_load())
        else:
            for j in range(n - 1, i - 1, -1):
                var a = j * w
                var carry = UInt64(0)
                for t in range(w):
                    var x = p.unsafe_offset(a + t).unsafe_load()
                    var sum = x + p.unsafe_offset(a + w + t).unsafe_load()
                    var low = sum + carry
                    carry = UInt64(1) if sum < x or low < sum else UInt64(0)
                    p.unsafe_offset(a + t).unsafe_store(low)
        if stop > 0:
            var sign = _limb_sign(p, i * w, w)
            if sign != 0:
                if last != 0 and sign != last:
                    count += 1
                last = sign
            if count + (1 if top != 0 and last != 0 and top != last else 0) >= stop:
                return stop
    if stop > 0:
        if top != 0 and last != 0 and top != last:
            count += 1
        return min(count, stop)
    return 0


def _shift(mut c: List[Integer], p: Integer) raises:
    # c(x) <- c(x + p).
    var n = len(c) - 1
    for i in range(n):
        for j in range(n - 1, i - 1, -1):
            c[j] += c[j + 1] * p


def _descartes_bound(c: List[Integer], cap: Int) raises -> Int:
    """Sign variations of (x + 1)^n c(1 / (x + 1)), capped: a bound on the roots of c in (0, 1).

    Round i of the Taylor shift fixes coefficient i. Dropping terms never adds
    variations, so the fixed prefix with the leading term c(0) already bounds
    the count from below, and the shift stops once that reaches the cap.
    """
    var n = len(c) - 1
    if n <= 0:
        return 0
    var w = _shift_width(c)
    var data = _load(c, w, True)
    return _taylor_words(data, n, w, cap)


def _root_endpoint(c: List[Integer]) raises -> Bool:
    # Whether c(0) = 0 or c(1) = 0.
    if not c[0]:
        return True
    var total = Integer(0)
    for i in range(len(c)):
        total += c[i].copy()
    return not total


def _halve(c: List[Integer]) raises -> List[Integer]:
    # 2^n c(x / 2) without the powers of two common to every coefficient.
    var n = len(c) - 1
    var out = List[Integer](capacity=n + 1)
    var common = Int.MAX
    for i in range(n + 1):
        out.append(c[i] << (n - i))
        if out[i]:
            common = min(common, _twos(out[i]))
    if common != Int.MAX and common > 0:
        for i in range(n + 1):
            out[i] = out[i] >> common
    return out^


def _point(index: Integer, depth: Int, k: Int) raises -> Rational:
    # index 2^(k - depth), in lowest terms.
    if not index:
        return Rational(0)
    var e = k - depth
    if e >= 0:
        return Rational(index << e)
    var t = min(_twos(index), -e)
    return Rational(_numerator=index >> t, _denominator=Integer(1) << (-e - t))


def _scaled(g: List[Integer], k: Int, negative: Bool) raises -> List[Integer]:
    # A positive multiple of g(+-2^k x): its roots in (0, 1) are those of g in
    # (0, 2^k) or (-2^k, 0), scaled.
    var n = len(g) - 1
    var b = List[Integer](capacity=n + 1)
    for i in range(n + 1):
        var v = g[i] << (k * i if k >= 0 else -k * (n - i))
        b.append(-v if negative and i % 2 == 1 else v^)
    return b^


def _newton_step(p: List[Integer], v: Int, m: Int, end: Int) raises -> Optional[Integer]:
    """The grid point of 2^(m+1) parts nearest to the Newton estimate from end 0 or 1.

    A cluster of v roots near x attracts end - v p(end) / p'(end).
    """
    var value = p[0].copy()
    var slope = p[1].copy() if len(p) > 1 else Integer(0)
    if end == 1:
        value = Integer(0)
        slope = Integer(0)
        for i in range(len(p)):
            value += p[i].copy()
            if i:
                slope += p[i] * Integer(i)
    if not slope:
        return None
    # round((end p' - v p) 2^(m+1) / p').
    var numerator = (slope * Integer(end) - value * Integer(v)) << (m + 2)
    numerator += slope.copy()
    var denominator = slope << 1
    if denominator.sign() < 0:
        numerator = -numerator
        denominator = -denominator
    if numerator.sign() < 0:
        return None
    var grid = numerator // denominator
    if grid > (Integer(1) << (m + 1)):
        return None
    return grid^


def _window(p: List[Integer], start: Integer, m: Int) raises -> List[Integer]:
    """A positive multiple of p((start + 4 x) / 2^(m+1)), without common powers of two."""
    var n = len(p) - 1
    var c = List[Integer](capacity=n + 1)
    for i in range(n + 1):
        c.append(p[i] << ((m + 1) * (n - i)))
    if start:
        _shift(c, start)
    var common = Int.MAX
    for i in range(n + 1):
        if i:
            c[i] = c[i] << (2 * i)
        if c[i]:
            common = min(common, _twos(c[i]))
    if common != Int.MAX and common > 0:
        for i in range(n + 1):
            c[i] = c[i] >> common
    return c^


comptime _GUARD = 64
"""Fractional bits of the fixed-point grid evaluation."""


def _local_signs(c: List[Integer], log: Int, points: List[UInt64]) -> List[Int]:
    """Exact signs of c at a / 2^(log + 1) for each a in points, with log + 1 <= 62.

    Fixed-point Horner, t <- floor(t a / 2^(log + 1)) + c_i 2^64, gives
    2^64 c(x) within n units for 0 < x < 1, on limbs as wide as the
    coefficients: that settles every sign whose value exceeds n units.
    The rest are evaluated exactly (_exact_signs).
    """
    var n = len(c) - 1
    var l = log + 1
    var words = _Words(c)
    var top = words.top
    var w = (top + _GUARD + 2 * _bits(Integer(n + 1)) + 4) // 64 + 1
    var t = List[UInt64](length=w + 1, fill=0)
    var p = t.unsafe_ptr()
    var signs = List[Int](capacity=len(points))
    var unsure = List[UInt64]()
    for j in range(len(points)):
        var a = points[j]
        for k in range(w + 1):
            p.unsafe_offset(k).unsafe_store(0)
        _add_shifted(p, w, words, n, _GUARD)
        for i in range(n - 1, -1, -1):
            # Sign-extend, multiply by a modulo 2^(64 (w + 1)), and shift
            # right arithmetically: floor(t a / 2^l) fits w limbs again.
            p.unsafe_offset(w).unsafe_store(UInt64.MAX if (p.unsafe_offset(w - 1).unsafe_load() >> 63) != 0 else UInt64(0))
            var carry = UInt128(0)
            for k in range(w + 1):
                var product = UInt128(p.unsafe_offset(k).unsafe_load()) * UInt128(a) + carry
                p.unsafe_offset(k).unsafe_store(UInt64(product & UInt128(UInt64.MAX)))
                carry = product >> 64
            for k in range(w):
                p.unsafe_offset(k).unsafe_store((p.unsafe_offset(k).unsafe_load() >> UInt64(l)) | (p.unsafe_offset(k + 1).unsafe_load() << UInt64(64 - l)))
            _add_shifted(p, w, words, i, _GUARD)
        var sign = _limb_sign(p, 0, w)
        if _beyond(p, w, sign, n):
            signs.append(sign)
        else:
            signs.append(2)
            unsure.append(a)
    if len(unsure):
        var exact = _exact_signs(c, log, unsure)
        var next = 0
        for j in range(len(signs)):
            if signs[j] == 2:
                signs[j] = exact[next]
                next += 1
    return signs^


def _beyond(p: Pointer[mut=True, UInt64, _], w: Int, sign: Int, bound: Int) -> Bool:
    # Whether the w-limb two's complement value exceeds bound in magnitude.
    if sign == 0:
        return False
    var flip = UInt64.MAX if sign < 0 else UInt64(0)
    for k in range(1, w):
        if (p.unsafe_offset(k).unsafe_load() ^ flip) != 0:
            return True
    var low = p.unsafe_offset(0).unsafe_load() ^ flip
    # For a negative value the magnitude is ~low + 1.
    return low > UInt64(bound) if sign > 0 else low >= UInt64(bound)


def _exact_signs(c: List[Integer], log: Int, points: List[UInt64]) -> List[Int]:
    """Exact signs of c at a / 2^(log + 1) for each a in points, with log + 1 <= 62.

    2^((log + 1) n) c(a / 2^(log + 1)) is Horner's
    t <- a t + c_i 2^((log + 1)(n - i)), on fixed-width two's complement
    limbs wide enough for the whole sum, so every step is a multiply by
    one word and a shifted add, without allocation.
    """
    var n = len(c) - 1
    var l = log + 1
    var words = _Words(c)
    var top = words.top
    var w = (top + l * n + 2 * _bits(Integer(n + 1)) + 2) // 64 + 1
    # After step i the partial sum needs only its first used[i] limbs; work
    # modulo 2^(64 used[i]) and sign-extend as the width grows.
    var used = List[Int](length=n + 1, fill=w)
    for i in range(n + 1):
        used[i] = min(w, (top + l * (n - i) + 2 * _bits(Integer(n + 1)) + 2) // 64 + 1)
    var t = List[UInt64](length=w, fill=0)
    var p = t.unsafe_ptr()
    var signs = List[Int](capacity=len(points))
    for j in range(len(points)):
        var a = points[j]
        var active = used[n]
        for k in range(active):
            p.unsafe_offset(k).unsafe_store(0)
        _add_shifted(p, active, words, n, 0)
        for i in range(n - 1, -1, -1):
            var wider = used[i]
            if wider > active:
                var fill = UInt64.MAX if (p.unsafe_offset(active - 1).unsafe_load() >> 63) != 0 else UInt64(0)
                for k in range(active, wider):
                    p.unsafe_offset(k).unsafe_store(fill)
                active = wider
            var carry = UInt128(0)
            for k in range(active):
                var product = UInt128(p.unsafe_offset(k).unsafe_load()) * UInt128(a) + carry
                p.unsafe_offset(k).unsafe_store(UInt64(product & UInt128(UInt64.MAX)))
                carry = product >> 64
            _add_shifted(p, active, words, i, l * (n - i))
        signs.append(_limb_sign(p, 0, active))
    return signs^


def _unit_sign(c: List[Integer], x: Rational) raises -> Int:
    """The sign of c at x for 0 <= x <= 1, as _sign computes it with a tighter error bound.

    With |x| <= 1, the n floor roundings of fixed-point Horner err by less
    than n units in all, so log2(n) + 1 bits of error suffice, where _sign
    allows for |x| up to 2.
    """
    var n = len(c) - 1
    if n <= 0:
        return c[0].sign() if n == 0 else 0
    var p = x.numerator()
    var q = x.denominator()
    if not p:
        return c[0].sign()
    var shift = _power_of_two(q)
    var error = _bits(Integer(n)) + 1
    var exact_bits = n * (_bits(p) + _bits(q))
    var extra = 64 + _bits(q)
    while error + extra < exact_bits:
        var t = _fixed(c, p, q, shift, error + extra)
        if t and _bits(t) - 1 > error + 1:
            return t.sign()
        extra *= 4
    return _exact_sign(c, p, q)


def _grid(count: Int, crowded: Bool) -> Tuple[List[UInt64], Int]:
    """count points in (0, 1) as numerators over 2^(log + 1), and log.

    Plain points are the midpoints of 2^log >= count equal parts. Crowded
    points are 3 s^2 - 2 s^3 at the midpoints s of count equal parts:
    spacing 1.5 / count in the middle and about 1 / count^2 at the ends,
    where roots of polynomials with only real roots crowd. Exact integer
    arithmetic places them, the same on every machine.
    """
    if not crowded:
        var parts = 1
        var plain = 0
        while parts < count:
            parts *= 2
            plain += 1
        var midpoints = List[UInt64](capacity=parts)
        for j in range(parts):
            midpoints.append(UInt64(2 * j + 1))
        return (midpoints^, plain)
    var log = 2
    while (1 << (log + 1)) < 8 * count * count:
        log += 1
    var scale = UInt128(1) << UInt128(log + 1)
    var m = UInt128(2 * count)
    var cube = m * m * m
    var points = List[UInt64](capacity=count)
    for j in range(count):
        # s = (2j + 1) / (2 count): x = s^2 (3 - 2 s) = u^2 (3 m - 2 u) / m^3.
        var u = UInt128(2 * j + 1)
        var x = u * u * (UInt128(3) * m - UInt128(2) * u) * scale // cube
        points.append(UInt64(x) | 1)
    return (points^, log)


def _grid_brackets(poly: List[Integer], index: Integer, depth: Int, width: Int, points: Int, crowded: Bool) -> Tuple[List[Integer], Int, Bool]:
    """Sign-change brackets in the node (index / 2^depth, index / 2^depth + 2^-width) with local polynomial poly.

    Signs come at the node's ends and at `points` interior points (_grid).
    Returns the numerators of the brackets' ends over 2^D, lowest first, D,
    and whether no interior point was a root; a root there ends the search
    early.
    """
    if points > 1 << 20:
        return (List[Integer](), 0, False)
    var grid = _grid(points, crowded)
    ref at_points = grid[0]
    var log = grid[1]
    var inner = _local_signs(poly, log, at_points)
    var d = max(depth, width + log + 1)
    var base: Integer
    var unit: Integer
    try:
        base = index << (d - depth)
        unit = Integer(1) << (d - width - log - 1)
    except:
        return (List[Integer](), 0, False)
    # The node's lower end, the points, and its upper end; an end that is a
    # root (recorded where it was split) drops out.
    var total = Integer(0)
    for i in range(len(poly)):
        try:
            total += poly[i].copy()
        except:
            return (List[Integer](), 0, False)
    var ends = List[UInt64]()
    var clean = True
    var last = 0
    var at = UInt64(0)
    var parts = len(at_points)
    for j in range(-1, parts + 1):
        var sign: Int
        if j < 0:
            sign = poly[0].sign()
        elif j < parts:
            sign = inner[j]
            if sign == 0:
                clean = False
                break
        else:
            sign = total.sign()
        if sign == 0:
            continue
        var position = at_points[j] if j >= 0 and j < parts else (UInt64(0) if j < 0 else UInt64(1) << UInt64(log + 1))
        if last != 0 and sign != last:
            ends.append(at)
            ends.append(position)
        last = sign
        at = position
    var out = List[Integer](capacity=len(ends))
    try:
        for i in range(len(ends)):
            out.append(base + unit * Integer(ends[i]))
    except:
        return (List[Integer](), 0, False)
    return (out^, d, clean)


def _pair_point(b: List[Integer], db: List[Integer], low: Rational, high: Rational, start: Rational, limit: Int) raises -> Optional[Rational]:
    """A dyadic point in (low, high) where b has the opposite sign to b(low), or None.

    With two close roots r +- d and the other roots far away, the Newton step
    x - 2 b(x) / b'(x) for a double root moves an error e from their centre
    to about d^2 / e plus a term in e^2: iterated from outside, it falls
    between the roots, where b changes sign. For a complex pair the step
    stops shrinking, which ends the search. Steps use fixed-point values;
    only the final sign must be exact.
    """
    var n = len(b) - 1
    var left = _sign(b, low)
    if left == 0 or _sign(b, high) != left:
        return None
    var x = start.copy()
    var previous = -(1 << 30)
    for _ in range(64):
        var p = x.numerator()
        var q = x.denominator()
        var s = _power_of_two(q)
        if s < 0:
            return None
        var error = _error_bits(n, p, q)
        var scale = 2 * s + 128
        var value: Integer
        var slope: Integer
        while True:
            value = _fixed(b, p, q, s, scale)
            slope = _fixed(db, p, q, s, scale)
            if _bits(value) > error + 8 and _bits(slope) > error + 8:
                break
            if scale > 4 * limit:
                return None
            scale *= 2
        # The step 2 value / slope is near 2^-small.
        var small = _bits(slope) - _bits(value) - 1
        if small <= previous:
            return None
        previous = small
        var t = max(s, 2 * small + 64)
        if t > limit:
            return None
        var step = (value << (t + 1)) // slope
        x = Rational((p << (t - s)) - step, Integer(1) << t)
        if x <= low or x >= high:
            return None
        var sign = _sign(b, x)
        if sign == 0:
            return None
        if sign != left:
            return x^
    return None


def _isolate_side(g: List[Integer], k: Int, negative: Bool, mut leaves: List[DescartesLeaf], mut roots: List[QQ], mut brackets: List[RootInterval], mut splits: Int, max_splits: Int, mut budget: Budget) raises PolyError:
    """Append the leaves, exact roots and brackets of g on (0, 2^k), or on (-2^k, 0), in increasing order.

    A node whose bound did not drop below its parent's tries a Newton step
    first (Sagraloff): the estimate of a cluster's position from either end
    selects a window of 4 / 2^(m+1) of the interval, and when the window's
    bound equals the node's, the rest of the node holds no roots, because
    Descartes' bounds of disjoint subintervals add up to at most the whole
    one's. Success squares the grid size for the window and failure takes
    its square root.

    A stalled node with bound two also looks for a point between two close
    real roots (_pair_point); signs there resolve the node without the deep
    transforms that separating the pair by bisection would need.
    """
    var side_leaves = List[DescartesLeaf]()
    var side_roots = List[QQ]()
    var side_brackets = List[RootInterval]()
    try:
        var b = _scaled(g, k, negative)
        var db = List[Integer](capacity=len(b))
        var top = 0
        for i in range(len(b)):
            top = max(top, _bits(b[i]))
            if i:
                db.append(b[i] * Integer(i))
        # Distinct roots of an integer polynomial of degree n and height H are
        # more than (n H)^-n apart, about.
        var n_b = len(b) - 1
        var limit = n_b * (top + _bits(Integer(n_b)) + 2) + 64
        var stack = List[_Node]()
        # The top starts as if a search had seen no sign changes under its
        # exact bound: the first search waits until 2^width reaches it. A
        # bound below half the degree comes with too few real roots for an
        # early search to pay, so the wait is then for the degree.
        budget.charge(UInt64(len(b) * len(b) // 2 + 1))
        var top_bound = _descartes_bound(b, len(b))
        if 2 * top_bound < n_b:
            top_bound = n_b
        stack.append(_Node(b.copy(), Integer(0), 0, 0, _SEARCH, 2, 0, -1, 0, 0, top_bound, Integer(0)))
        while len(stack):
            budget.charge(1)
            var node = stack.pop()
            if node.kind == _ROOT:
                _allocation[QQ](1, budget)
                side_roots.append(QQ(_point(node.index, node.depth, k)))
                continue
            var v = 0
            # The bound is capped a little above the parent's, which keeps
            # the shift short; below the cap it is exact.
            var cap = 2 if node.parent < 2 else node.parent + 1
            if node.kind == _KNOWN:
                v = Int(node.end)
                cap = v + 1
                node.kind = _SEARCH
            elif node.kind == _SEARCH:
                var n = len(node.poly) - 1
                budget.charge(UInt64(n * n // 2 + 1))
                v = _descartes_bound(node.poly, cap)
                if v >= 2 and splits >= max_splits:
                    v = _descartes_bound(node.poly, n + 1)
            # A root's leaf needs endpoints where g is nonzero, so that its
            # bracket starts with a sign change.
            var resolved = node.kind == _EMPTY or v == 0 or (v == 1 and not _root_endpoint(node.poly))
            if resolved or splits >= max_splits:
                _allocation[DescartesLeaf](1, budget)
                var low = _point(node.index, node.depth, k)
                var high = _point(node.end.copy() if node.kind == _EMPTY else node.index + (Integer(1) << (node.depth - node.width)), node.depth, k)
                if v == 1:
                    _allocation[RootInterval](1, budget)
                    side_brackets.append(RootInterval(QQ(low.copy()), QQ(high.copy())))
                side_leaves.append(DescartesLeaf(QQ(low^), QQ(high^), v, v <= 1))
                continue
            # Signs at about 4 v points may already show v sign changes,
            # so that v brackets resolve the node at once; when most show,
            # 6 v points crowded at the ends follow. A failed search saw
            # some changes; each halving of the width about doubles them, so
            # the next search waits until they could reach the exact bound of
            # the last one, which bounds this node's.
            var due = False
            if v >= 2:
                var gap = node.width - node.seen_width
                var needed = v if v < cap else max(v, node.seen_bound)
                due = gap >= 62 or max(node.seen, 1) << gap >= needed
            if due:
                var exact = v if v < cap else _descartes_bound(node.poly, len(node.poly))
                budget.charge(UInt64(len(b) * len(b) // 2 + 1))
                var points = 4 * exact
                var grid = _grid_brackets(node.poly, node.index, node.depth, node.width, points, False)
                budget.charge(UInt64(len(b) * points))
                var found = len(grid[0]) // 2
                if grid[2] and found != exact and exact >= 8 and 4 * found >= 3 * exact:
                    # Most roots showed; the rest may crowd at the ends.
                    grid = _grid_brackets(node.poly, node.index, node.depth, node.width, 6 * exact, True)
                    budget.charge(UInt64(len(b) * 6 * exact))
                    found = len(grid[0]) // 2
                if not grid[2] or found != exact:
                    node.seen = found
                    node.seen_width = node.width
                    node.seen_bound = exact
                else:
                    v = exact
                    _allocation[RootInterval](v, budget)
                    var ends = grid[0].copy()
                    var d = grid[1]
                    for j in range(v):
                        side_brackets.append(RootInterval(QQ(_point(ends[2 * j], d, k)), QQ(_point(ends[2 * j + 1], d, k))))
                    var a = QQ(_point(node.index, node.depth, k))
                    var z = QQ(_point(node.index + (Integer(1) << (node.depth - node.width)), node.depth, k))
                    _allocation[DescartesLeaf](1, budget)
                    side_leaves.append(DescartesLeaf(a^, z^, v, True))
                    continue
            if v == 2 and node.parent == 2 and (node.paired < 0 or node.width >= 2 * node.paired + 8) and not _root_endpoint(node.poly):
                var low = Rational(node.index.copy(), Integer(1) << node.depth)
                var high = Rational(node.index + (Integer(1) << (node.depth - node.width)), Integer(1) << node.depth)
                var cut = _pair_point(b, db, low, high, low, limit)
                if not cut:
                    cut = _pair_point(b, db, low, high, high, limit)
                if cut:
                    budget.charge(UInt64(len(b)))
                    var x = cut.take()
                    var t = _power_of_two(x.denominator())
                    var middle = QQ(_point(x.numerator(), t, k))
                    var a = QQ(_point(node.index, node.depth, k))
                    var z = QQ(_point(node.index + (Integer(1) << (node.depth - node.width)), node.depth, k))
                    _allocation[RootInterval](2, budget)
                    side_brackets.append(RootInterval(a.copy(), middle.copy()))
                    side_brackets.append(RootInterval(middle^, z.copy()))
                    _allocation[DescartesLeaf](1, budget)
                    side_leaves.append(DescartesLeaf(a^, z^, 2, True))
                    continue
                node.paired = node.width
            splits += 1
            var n = len(node.poly) - 1
            if v >= 2 and v == node.parent:
                var jumped = False
                var m = node.newton
                var near = _newton_step(node.poly, v, m, 0)
                var far = _newton_step(node.poly, v, m, 1)
                if near and far:
                    # The two estimates differ by about the cluster's size:
                    # coarsen the grid until one window holds both.
                    var gap = near.value() - far.value()
                    var bits = abs(gap).magnitude_bit_length()
                    if bits > 1:
                        var coarse = max(2, m - bits + 1)
                        var mean = (near.value() + far.value()) >> (m - coarse + 1)
                        near = Optional[Integer](mean^)
                        m = coarse
                elif far:
                    near = far^
                if near:
                    var parts = Integer(1) << (m + 1)
                    var start = near.value() - Integer(2)
                    if start.sign() < 0:
                        start = Integer(0)
                    if start > parts - Integer(4):
                        start = parts - Integer(4)
                    budget.charge(UInt64(n * n + 1))
                    var window = _window(node.poly, start, m)
                    if _descartes_bound(window, v + 1) == v:
                        # The window and the root-free pieces on either side.
                        var depth = max(node.depth, node.width + m + 1)
                        var base = node.index << (depth - node.depth)
                        var unit = Integer(1) << (depth - node.width - m - 1)
                        var low = base + unit * start
                        var high = base + unit * (start + Integer(4))
                        var width = node.width + m - 1
                        _allocation[_Node](3, budget)
                        if start + Integer(4) < parts:
                            var edge = base + (Integer(1) << (depth - node.width))
                            stack.append(_Node(List[Integer](), high.copy(), depth, 0, _EMPTY, 0, 0, -1, -1, 0, 0, edge^))
                        # The window's bound is known: v.
                        stack.append(_Node(window^, low.copy(), depth, width, _KNOWN, min(2 * m, _MAX_NEWTON), v, node.paired, node.seen, node.seen_width, node.seen_bound, Integer(v)))
                        if start:
                            stack.append(_Node(List[Integer](), base.copy(), depth, 0, _EMPTY, 0, 0, -1, -1, 0, 0, low.copy()))
                        jumped = True
                    else:
                        node.newton = m
                if jumped:
                    continue
                node.newton = max(2, node.newton // 2)
            var left = _halve(node.poly)
            var right = left.copy()
            budget.charge(UInt64(n * n // 2 + 1))
            _shift_one(right)
            var depth = max(node.depth, node.width + 1)
            var low = node.index << (depth - node.depth)
            var middle = low + (Integer(1) << (depth - node.width - 1))
            var root = not right[0]
            _allocation[_Node](3, budget)
            stack.append(_Node(right^, middle.copy(), depth, node.width + 1, _SEARCH, node.newton, v, node.paired, node.seen, node.seen_width, node.seen_bound, Integer(0)))
            if root:
                stack.append(_Node(List[Integer](), middle.copy(), depth, 0, _ROOT, 0, 0, -1, -1, 0, 0, Integer(0)))
            stack.append(_Node(left^, low^, depth, node.width + 1, _SEARCH, node.newton, v, node.paired, node.seen, node.seen_width, node.seen_bound, Integer(0)))
    except e:
        raise _error(e)
    if negative:
        # Mirror the leaves of g(-x), so they increase.
        for i in range(len(side_leaves) - 1, -1, -1):
            ref leaf = side_leaves[i]
            leaves.append(DescartesLeaf(leaf.high.neg(budget), leaf.low.neg(budget), leaf.variations, leaf.resolved))
        for i in range(len(side_roots) - 1, -1, -1):
            roots.append(side_roots[i].neg(budget))
        for i in range(len(side_brackets) - 1, -1, -1):
            ref bracket = side_brackets[i]
            brackets.append(RootInterval(bracket.high.neg(budget), bracket.low.neg(budget)))
    else:
        for i in range(len(side_leaves)):
            leaves.append(side_leaves[i].copy())
        for i in range(len(side_roots)):
            roots.append(side_roots[i].copy())
        for i in range(len(side_brackets)):
            brackets.append(side_brackets[i].copy())


def _integer_list(p: DensePoly[ZZ]) -> List[Integer]:
    var values = List[Integer](capacity=len(p.coefficients))
    for c in p.coefficients:
        values.append(c.value.copy())
    return values^


def _squarefree(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> Tuple[DensePoly[ZZ], DensePoly[ZZ]]:
    """The primitive positive-leading g = f / gcd(f, f') and gcd(f, f')."""
    var p = _primitive(_integer_dense(f, budget), budget)
    if p.coefficients[p.degree()].value.sign() < 0:
        p = p.neg(budget)
    var one = DensePoly[ZZ].from_coefficients(p.variable, NoContext(), [ZZ(Integer(1))], budget)
    if p.degree() < 2:
        return (p^, one^)
    _allocation[ZZ](p.degree(), budget)
    var dp = List[ZZ](capacity=p.degree())
    try:
        for i in range(1, len(p.coefficients)):
            dp.append(ZZ(p.coefficients[i].value * Integer(i)))
    except e:
        raise _error(e)
    var derivative = DensePoly[ZZ].from_coefficients(p.variable, NoContext(), dp^, budget)
    var h = _univariate_gcd(p, derivative, budget)
    if h.degree() == 0:
        return (p^, one^)
    var g = _trial_divide(p, h, budget)
    if not g:
        raise PolyError(PolyErrorKind.DOMAIN, "gcd(f, f') does not divide f")
    return (g.take(), h^)


def _bound_exponent(g: List[Integer]) -> Int:
    # Every root of g has absolute value below 2^k. A monomial, whose only
    # root is zero, gets k = 0 rather than Fujiwara's empty-maximum sentinel.
    for i in range(len(g) - 1):
        if g[i]:
            return _fujiwara_exponent(g)
    return 0


def _descartes_isolate(f: DensePoly[QQ], options: NumericOptions, mut budget: Budget) raises PolyError -> DescartesRootResult:
    """Isolate the distinct real roots of a nonzero rational polynomial in brackets of their leaves."""
    if f.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "real roots of the zero polynomial")
    var parts = _squarefree(f, budget)
    var g = _integer_list(parts[0])
    var k = _bound_exponent(g)
    var leaves = List[DescartesLeaf]()
    var roots = List[QQ]()
    var splits = 0
    if len(g) < 2:
        leaves.append(DescartesLeaf(QQ(Rational(-1)), QQ(Rational(1)), 0, True))
        var certificate = DescartesCertificate(QQ(Rational(-1)), QQ(Rational(1)), leaves^, roots^, List[RootInterval]())
        return DescartesRootResult(NumericOutcome.SUCCESS, True, True, True, parts[0].copy(), parts[1].copy(), 0, certificate^)
    var brackets = List[RootInterval]()
    # A root at zero is the shared endpoint of both halves.
    _isolate_side(g, k, True, leaves, roots, brackets, splits, options.max_subdivisions, budget)
    if not g[0]:
        roots.append(QQ(Rational(0)))
    _isolate_side(g, k, False, leaves, roots, brackets, splits, options.max_subdivisions, budget)
    var a: Rational
    var b: Rational
    try:
        a = _point(Integer(-1), 0, k)
        b = _point(Integer(1), 0, k)
    except e:
        raise _error(e)
    var complete = True
    for leaf in leaves:
        complete = complete and leaf.resolved
    var certificate = DescartesCertificate(QQ(a^), QQ(b^), leaves^, roots^, brackets^)
    var result = DescartesRootResult(NumericOutcome.SUCCESS, complete, False, complete, parts[0].copy(), parts[1].copy(), k, certificate^)
    return _descartes_refine(result, options, budget, splits)


def _descartes_refine(previous: DescartesRootResult, options: NumericOptions, mut budget: Budget, used: Int = 0) raises PolyError -> DescartesRootResult:
    """Tighten every bracket of a previous result until it meets the options' widths.

    Each step of quadratic interval refinement, with its bisection fallback,
    counts as a subdivision; with `used` already spent, the rest of the
    options' allowance bounds them, and brackets it leaves wide make the
    result inaccurate.
    """
    options._validate()
    var result = previous.copy()
    var steps = used
    var accurate = True
    if len(result.certificate.brackets):
        var g = _integer_list(result.squarefree)
        var dg = List[Integer](capacity=len(g))
        try:
            for i in range(1, len(g)):
                dg.append(g[i] * Integer(i))
        except e:
            raise _error(e)
        var oracle = _RootOracle(g^, dg^)
        # Absolute widths alone give a fixed target, which also bounds the
        # refinement grid as _tighten does.
        var target = Optional[Rational]()
        if options.relative_bits == 0 and options.absolute_bits > 0:
            try:
                target = Rational(Integer(1), Integer(1) << options.absolute_bits)
            except e:
                raise _error(e)
        for i in range(len(result.certificate.brackets)):
            budget.charge(1)
            ref bracket = result.certificate.brackets[i]
            if _bracket_accurate(bracket, options, budget):
                continue
            var enclosure = _enclosure(oracle, bracket.low, bracket.high, budget)
            if _check_sign(oracle, bracket.high.value, budget) != -enclosure.left:
                raise PolyError(PolyErrorKind.DOMAIN, "a root bracket lost its sign change")
            while not enclosure.exact and not _accurate(NumericBox(QQ(enclosure.low.copy()), QQ(enclosure.high.copy()), QQ(Rational(0)), QQ(Rational(0))), options, budget):
                if steps >= options.max_subdivisions:
                    break
                steps += 1
                budget.charge(1)
                if target:
                    enclosure.grid = max(2, min(enclosure.grid, _ratio_bits(enclosure, target.value()) + 1))
                if _qir_step(oracle, enclosure, budget):
                    enclosure.grid = min(enclosure.grid * 2, 1 << 20)
                else:
                    enclosure.grid = max(2, enclosure.grid // 2)
                    if not enclosure.exact:
                        var middle: Rational
                        try:
                            middle = (enclosure.low + enclosure.high) / Rational(2)
                        except e:
                            raise _error(e)
                        _ = _side(oracle, enclosure, QQ(middle^), budget)
            _charge_size((_bits(enclosure.low.numerator()) + _bits(enclosure.low.denominator())) // 8 + 1, budget)
            _charge_size((_bits(enclosure.high.numerator()) + _bits(enclosure.high.denominator())) // 8 + 1, budget)
            bracket.low = QQ(enclosure.low.copy())
            bracket.high = QQ(enclosure.high.copy())
            accurate = accurate and _bracket_accurate(bracket, options, budget)
    result.accurate = accurate
    result.outcome = _outcome(result.complete, accurate)
    return result^


def _interval_bound(g: List[Integer], low: Rational, high: Rational, cap: Int) raises -> Int:
    """Descartes' bound for the roots of g in the open interval (low, high)."""
    # With low = p / q and high = r / q, q^n g((p + (r - p) x) / q) is
    # G(p + w x) for G(y) = q^n g(y / q) and w = r - p.
    var q = low.denominator()
    var dh = high.denominator()
    if q != dh:
        q = q * dh // _gcd(q, dh)
    var p = low.numerator() * (q // low.denominator())
    var w = high.numerator() * (q // dh) - p
    var n = len(g) - 1
    var c = List[Integer](capacity=n + 1)
    var power = Integer(1)
    for _ in range(n + 1):
        c.append(Integer(0))
    for i in range(n, -1, -1):
        c[i] = g[i] * power
        power = power * q
    _shift(c, p)
    power = Integer(1)
    for i in range(n + 1):
        c[i] = c[i] * power
        power = power * w
    return _descartes_bound(c, cap)


def _leaf_bound(c: List[Integer], leaf: DescartesLeaf, mut budget: Budget) raises PolyError -> Int:
    try:
        return _interval_bound(c, leaf.low.value, leaf.high.value, max(leaf.variations + 1, 2))
    except e:
        raise _error(e)


def _gcd(a: Integer, b: Integer) raises -> Integer:
    var x = abs(a)
    var y = abs(b)
    while y:
        var t = x % y
        x = y^
        y = t^
    return x^


def check_descartes(g: DensePoly[ZZ], certificate: DescartesCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Verify a Descartes certificate for an integer polynomial without searching.

    Checks that the leaves partition the certificate's interval at nonroots
    and at the listed exact roots, recomputes every leaf's Descartes bound,
    and checks the brackets of each resolved leaf: disjoint, inside it, and
    each with a sign change or a root. Unresolved leaves are accepted; they
    leave a result incomplete, not wrong.

    Args:
        g: The polynomial the certificate is for; squarefree in practice, though the checks do not need it.
        certificate: The proposed certificate.
        budget: The shared work allowance.

    Returns:
        Whether every claim of the certificate holds for g.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    var c = _integer_list(g)
    var oracle = _RootOracle(c.copy(), List[Integer]())
    ref leaves = certificate.leaves
    if not len(leaves) or leaves[0].low != certificate.low or leaves[len(leaves) - 1].high != certificate.high:
        return False
    var next_root = 0
    var next = 0
    for i in range(len(leaves)):
        ref leaf = leaves[i]
        var root = Optional[QQ]()
        if i:
            if leaves[i - 1].high != leaf.low:
                return False
            if next_root < len(certificate.exact_roots) and certificate.exact_roots[next_root] == leaf.low:
                root = certificate.exact_roots[next_root].copy()
                next_root += 1
        var count = leaf.variations if leaf.resolved else 0
        if next + count > len(certificate.brackets):
            return False
        if not _check_leaf(c, oracle, leaf, root, certificate.brackets, next, i == 0, budget):
            return False
        next += count
    if _check_sign(oracle, certificate.high.value, budget) == 0:
        return False
    return next_root == len(certificate.exact_roots) and next == len(certificate.brackets)


def _check_leaf(c: List[Integer], oracle: _RootOracle, leaf: DescartesLeaf, root: Optional[QQ], brackets: List[RootInterval], first: Int, outer: Bool, mut budget: Budget) raises PolyError -> Bool:
    """Check one leaf: its lower end against the claimed exact root, its recomputed bound, and its brackets.

    The lower end of the outer leaf is never a root; a shared lower end is
    a root exactly when `root` claims it.
    """
    budget.charge(UInt64(len(c) * len(c) // 2 + 1))
    if leaf.variations < 0 or (leaf.variations <= 1 and not leaf.resolved) or _compare(leaf.low, leaf.high, budget) >= 0:
        return False
    var s = _check_sign(oracle, leaf.low.value, budget)
    if root:
        if outer or s != 0 or root.value() != leaf.low:
            return False
    elif s == 0:
        return False
    if len(c) > 1 and _leaf_bound(c, leaf, budget) != leaf.variations:
        return False
    if not leaf.resolved:
        return True
    # Disjoint brackets inside the leaf, each with an odd root count.
    var floor = leaf.low.copy()
    for j in range(leaf.variations):
        ref bracket = brackets[first + j]
        if _compare(floor, bracket.low, budget) > 0 or _compare(bracket.high, leaf.high, budget) > 0:
            return False
        var order = _compare(bracket.low, bracket.high, budget)
        if order > 0:
            return False
        if order == 0:
            if bracket.low == leaf.low or bracket.low == leaf.high or (j and bracket.low == floor) or _check_sign(oracle, bracket.low.value, budget) != 0:
                return False
        else:
            var left = _check_sign(oracle, bracket.low.value, budget)
            if left == 0 or _check_sign(oracle, bracket.high.value, budget) != -left:
                return False
        floor = bracket.high.copy()
    return True


def _check_header(f: DensePoly[QQ], squarefree: DensePoly[ZZ], cofactor: DensePoly[ZZ], bound_exponent: Int, mut budget: Budget) raises PolyError -> Bool:
    # g = f_prim / gcd(f, f') and gcd(f, f'), recomputed, and the bound.
    if f.is_zero():
        return False
    var parts = _squarefree(f, budget)
    if parts[0] != squarefree or parts[1] != cofactor:
        return False
    return bound_exponent >= _bound_exponent(_integer_list(squarefree))


def _outer(g: DensePoly[ZZ], bound_exponent: Int) raises PolyError -> Tuple[QQ, QQ]:
    # The certificate's interval: (-2^k, 2^k), or (-1, 1) for a constant.
    try:
        if g.degree() < 1:
            return (QQ(Rational(-1)), QQ(Rational(1)))
        return (QQ(_point(Integer(-1), 0, bound_exponent)), QQ(_point(Integer(1), 0, bound_exponent)))
    except e:
        raise _error(e)


def _check_descartes_roots(f: DensePoly[QQ], result: DescartesRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    """Verify a Descartes result for f without repeating its search."""
    options._validate()
    if not _check_header(f, result.squarefree, result.cofactor, result.bound_exponent, budget):
        return False
    var outer = _outer(result.squarefree, result.bound_exponent)
    if result.certificate.low != outer[0] or result.certificate.high != outer[1]:
        return False
    if not check_descartes(result.squarefree, result.certificate, budget):
        return False
    var complete = True
    for leaf in result.certificate.leaves:
        complete = complete and leaf.resolved
    var accurate = True
    for bracket in result.certificate.brackets:
        accurate = accurate and _bracket_accurate(bracket, options, budget)
    return result.complete == complete and result.isolated == complete and result.accurate == accurate and result.outcome == _outcome(complete, accurate)


def _bracket_accurate(bracket: RootInterval, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    return _accurate(NumericBox(bracket.low.copy(), bracket.high.copy(), QQ(Rational(0)), QQ(Rational(0))), options, budget)


def _outcome(complete: Bool, accurate: Bool) -> NumericOutcome:
    return NumericOutcome.SUCCESS if complete and accurate else NumericOutcome.SUBDIVISION_LIMIT


def descartes_intervals(certificate: DescartesCertificate, mut budget: Budget) raises PolyError -> List[RootInterval]:
    """List the certified roots of a Descartes certificate in increasing order.

    Args:
        certificate: A certificate, checked or produced by this module.
        budget: The shared work allowance.

    Returns:
        One interval per root: the point (r, r) for an exact root and the bracket otherwise; roots of unresolved leaves are omitted.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    ref exact = certificate.exact_roots
    _allocation[RootInterval](len(exact) + len(certificate.brackets), budget)
    var out = List[RootInterval]()
    var next_root = 0
    var next_bracket = 0
    for leaf in certificate.leaves:
        budget.charge(1)
        while next_root < len(exact) and _compare(exact[next_root], leaf.low, budget) <= 0:
            out.append(RootInterval(exact[next_root].copy(), exact[next_root].copy()))
            next_root += 1
        if leaf.resolved:
            for _ in range(leaf.variations):
                out.append(certificate.brackets[next_bracket].copy())
                next_bracket += 1
    while next_root < len(exact):
        out.append(RootInterval(exact[next_root].copy(), exact[next_root].copy()))
        next_root += 1
    return out^


def _descartes_roots[C: CoefficientRing](poly: DensePoly[C], options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> DescartesRootResult:
    """Isolate and refine the distinct real roots of an exact ZZ/QQ polynomial."""
    options._validate()
    var exact = _rational_polynomial(_exact_dense_boxes(poly, budget), budget)
    if not exact:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "Descartes isolation needs exact coefficients")
    var result = _descartes_isolate(exact.value(), options, budget)
    _require_success(result.outcome, complete_only)
    return result^


def _exact_input[C: CoefficientRing](poly: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    var exact = _rational_polynomial(_exact_dense_boxes(poly, budget), budget)
    if not exact:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "Descartes isolation needs exact coefficients")
    return exact.take()


def _descartes_resume(f: DensePoly[QQ], previous: DescartesRootResult, previous_options: NumericOptions, options: NumericOptions, mut budget: Budget, complete_only: Bool, trusted: Bool) raises PolyError -> DescartesRootResult:
    """Validate a previous result unless trusted, then tighten its brackets; an incomplete one restarts with the new options."""
    options._validate()
    if not trusted and not _check_descartes_roots(f, previous, previous_options, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "previous real-root evidence is invalid")
    var result = _descartes_refine(previous, options, budget) if previous.complete else _descartes_isolate(f, options, budget)
    _require_success(result.outcome, complete_only)
    return result^


@fieldwise_init
struct DescartesRootItem(Copyable, Movable):
    """One leaf of a streamed Descartes certificate with the evidence that belongs to it."""
    var leaf: DescartesLeaf
    """The next leaf, starting where the previous one ended."""
    var root: Optional[QQ]
    """The exact root at the leaf's lower end, shared with the previous leaf, if g vanishes there."""
    var brackets: List[RootInterval]
    """The leaf's brackets: as many as its bound when resolved, none otherwise."""


trait DescartesRootSink:
    """Consume a Descartes certificate leaf by leaf without retaining its prefix."""

    def begin(mut self, squarefree: DensePoly[ZZ], cofactor: DensePoly[ZZ], bound_exponent: Int, mut budget: Budget) raises PolyError:
        """Receive the squarefree part and the bound before any leaves.

        Args:
            squarefree: The primitive positive-leading g = f / gcd(f, f').
            cofactor: The cofactor gcd(f, f').
            bound_exponent: The k with every root in (-2^k, 2^k).
            budget: The shared transport allowance.

        Raises:
            PolyError: On transport or resource failure.
        """
        ...

    def write(mut self, var item: DescartesRootItem, mut budget: Budget) raises PolyError:
        """Consume the next leaf with its exact root and brackets.

        Args:
            item: The next owned leaf in increasing order.
            budget: The shared transport allowance.

        Raises:
            PolyError: On transport or resource failure.
        """
        ...


@fieldwise_init
struct DescartesStreamChecker(DescartesRootSink, Movable):
    """Verify a streamed Descartes certificate leaf by leaf, retaining only the current edge."""
    var _exact: DensePoly[QQ]
    var _options: NumericOptions
    var _poly: List[Integer]
    var _oracle: _RootOracle
    var _edge: QQ
    var _end: QQ
    var _items: Int
    var _complete: Bool
    var _accurate: Bool
    var _started: Bool
    var _failed: Bool

    @staticmethod
    def from_dense[C: CoefficientRing](poly: DensePoly[C], options: NumericOptions, mut budget: Budget) raises PolyError -> Self:
        """Prepare independent checking for an exact polynomial.

        Parameters:
            C: ZZ or QQ.

        Args:
            poly: The original polynomial.
            options: The original requested guarantees.
            budget: The shared preparation allowance.

        Returns:
            A checker that accepts one header and then the leaves in order.

        Raises:
            PolyError: DOMAIN for invalid input or options, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
        """
        options._validate()
        var exact = _exact_input(poly, budget)
        if exact.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "real roots of the zero polynomial")
        return Self(exact^, options, List[Integer](), _RootOracle(List[Integer](), List[Integer]()), QQ(Rational(0)), QQ(Rational(0)), 0, True, True, False, False)

    def begin(mut self, squarefree: DensePoly[ZZ], cofactor: DensePoly[ZZ], bound_exponent: Int, mut budget: Budget) raises PolyError:
        """Check the squarefree part and the bound once, against the original polynomial.

        Args:
            squarefree: The claimed g = f / gcd(f, f').
            cofactor: The claimed gcd(f, f').
            bound_exponent: The claimed k with every root in (-2^k, 2^k).
            budget: The shared checking allowance.

        Raises:
            PolyError: DOMAIN for invalid or repeated headers, LIMIT or CAPACITY on work limits; failure invalidates this checker.
        """
        if self._failed or self._started:
            self._failed = True
            raise PolyError(PolyErrorKind.DOMAIN, "invalid root stream header")
        self._failed = True
        if not _check_header(self._exact, squarefree, cofactor, bound_exponent, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid root stream squarefree part or bound")
        var outer = _outer(squarefree, bound_exponent)
        self._poly = _integer_list(squarefree)
        self._oracle = _RootOracle(self._poly.copy(), List[Integer]())
        self._edge = outer[0].copy()
        self._end = outer[1].copy()
        self._started = True
        self._failed = False

    def write(mut self, var item: DescartesRootItem, mut budget: Budget) raises PolyError:
        """Check a leaf directly as a streaming destination.

        Args:
            item: The next owned leaf.
            budget: The shared checking allowance.

        Raises:
            PolyError: DOMAIN for rejected evidence, LIMIT or CAPACITY on work limits.
        """
        if not self.accept(item^, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid root stream witness")

    def accept(mut self, var item: DescartesRootItem, mut budget: Budget) raises PolyError -> Bool:
        """Check and discard a leaf while retaining coverage state.

        Args:
            item: The next owned leaf.
            budget: The shared checking allowance.

        Returns:
            Whether the leaf is valid; rejection permanently invalidates the checker.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits; failure invalidates the checker.
        """
        if self._failed:
            return False
        self._failed = True
        budget.charge(1)
        ref leaf = item.leaf
        if not self._started or leaf.low != self._edge or _compare(leaf.high, self._end, budget) > 0:
            return False
        var count = leaf.variations if leaf.resolved else 0
        if len(item.brackets) != count:
            return False
        if not _check_leaf(self._poly, self._oracle, leaf, item.root, item.brackets, 0, self._items == 0, budget):
            return False
        for bracket in item.brackets:
            self._accurate = self._accurate and _bracket_accurate(bracket, self._options, budget)
        self._complete = self._complete and leaf.resolved
        self._edge = leaf.high.copy()
        self._items = _sum(self._items, 1)
        self._failed = False
        return True

    def finish(self, receipt: RootStreamReceipt, mut budget: Budget) raises PolyError -> Bool:
        """Check full coverage and all receipt claims after consumption.

        Args:
            receipt: The producer's completed receipt.
            budget: The shared checking allowance.

        Returns:
            True only for a valid stream covering the whole interval with consistent outcome flags.

        Raises:
            PolyError: LIMIT on exhaustion.
        """
        budget.charge(1)
        if self._failed or not self._started or not self._items or self._edge != self._end:
            return False
        return receipt.items == self._items and receipt.complete == self._complete and receipt.isolated == self._complete and receipt.accurate == self._accurate and receipt.outcome == _outcome(self._complete, self._accurate)


def _descartes_roots_to[S: DescartesRootSink](f: DensePoly[QQ], options: NumericOptions, mut sink: S, mut budget: Budget, complete_only: Bool) raises PolyError -> RootStreamReceipt:
    """Isolate and refine, then stream the certificate leaf by leaf.

    The negative half is isolated outward from zero and mirrored, so the
    producer holds the certificate before streaming it; the sink need not.
    """
    var result = _descartes_isolate(f, options, budget)
    sink.begin(result.squarefree, result.cofactor, result.bound_exponent, budget)
    ref certificate = result.certificate
    var next_root = 0
    var next = 0
    for i in range(len(certificate.leaves)):
        ref leaf = certificate.leaves[i]
        var root = Optional[QQ]()
        if i and next_root < len(certificate.exact_roots) and certificate.exact_roots[next_root] == leaf.low:
            root = certificate.exact_roots[next_root].copy()
            next_root += 1
        var count = leaf.variations if leaf.resolved else 0
        _allocation[RootInterval](count, budget)
        var brackets = List[RootInterval](capacity=count)
        for j in range(count):
            brackets.append(certificate.brackets[next + j].copy())
        next += count
        sink.write(DescartesRootItem(leaf.copy(), root^, brackets^), budget)
    var receipt = RootStreamReceipt(result.outcome, result.complete, result.accurate, result.isolated, len(certificate.leaves))
    _require_success(result.outcome, complete_only)
    return receipt^


def _root_isolation(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> RootIsolation:
    """Every rational root exactly and disjoint open intervals narrower than 2^-40 for the others.

    Rational roots have denominators dividing the leading coefficient D of
    the primitive squarefree part. Brackets narrower than 1/(2 D^2) make any
    such root a continued-fraction convergent of the midpoint (Legendre),
    which _rational_root tries with exact evaluation.
    """
    if f.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "real roots of the zero polynomial")
    var exact = List[QQ]()
    var intervals = List[RootInterval]()
    if f.degree() == 0:
        return RootIsolation(exact^, intervals^)
    var parts = _squarefree(f, budget)
    var bound = parts[0].coefficients[parts[0].degree()].copy()
    var bits = max(41, 2 * bound.value.magnitude_bit_length() + 2)
    var options = NumericOptions(absolute_bits=bits, max_subdivisions=1 << 60)
    var result = _descartes_isolate(f, options, budget)
    if not result.complete or not result.accurate:
        raise PolyError(PolyErrorKind.LIMIT, "real root isolation did not complete")
    var g = _rational_from(parts[0], budget)
    for root in descartes_intervals(result.certificate, budget):
        if root.low == root.high:
            _allocation[QQ](1, budget)
            exact.append(root.low.copy())
            continue
        var rational = _rational_root(g, root.low, root.high, bound, budget)
        if rational:
            _allocation[QQ](1, budget)
            exact.append(rational.take())
        else:
            _allocation[RootInterval](1, budget)
            intervals.append(root.copy())
    return RootIsolation(exact^, intervals^)


def _rational_from(p: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    _allocation[QQ](len(p.coefficients), budget)
    var values = List[QQ](capacity=len(p.coefficients))
    for c in p.coefficients:
        values.append(QQ(Rational(c.value)))
    return DensePoly[QQ].from_coefficients(p.variable, NoContext(), values^, budget)


def _refine_interval(f: DensePoly[QQ], interval: RootInterval, width: QQ, mut budget: Budget) raises PolyError -> RootInterval:
    """Validate a one-root interval of f by Descartes' rule and tighten it below width by QIR.

    The primitive squarefree part g has the roots of f, each simple. Parts of
    the interval are split until each has Descartes bound 0 or 1, the bound
    being exact at 0 and 1 and settling for a squarefree g; split points that
    are roots count as roots. Exactly one root must remain. A bound-1 part has
    a sign change of g at its nonroot ends, which QIR narrows below 5 / 32 of width;
    an exact root m gets the part of the interval within width / 4 of m.
    """
    if _compare(interval.low, interval.high, budget) >= 0 or _sign_at(f, interval.low, budget) == 0 or _sign_at(f, interval.high, budget) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "root count requires ordered nonroot endpoints")
    var c = _integer_list(_squarefree(f, budget)[0])
    var dc: List[Integer]
    try:
        dc = _derivative(c)
    except e:
        raise _error(e)
    var oracle = _RootOracle(c.copy(), dc^)
    var found = _split_roots(c, oracle, interval, 2, budget)
    if len(found) != 1:
        raise PolyError(PolyErrorKind.DOMAIN, "refinement requires exactly one root")
    if _compare(interval.high.sub(interval.low, budget), width, budget) < 0:
        return interval.copy()
    var root = found[0].copy()
    if root.low != root.high:
        var enclosure = _enclosure(oracle, root.low, root.high, budget)
        _tighten_guarded(oracle, enclosure, width, budget)
        if not enclosure.exact:
            # Below width / 8 and rounded outward to a multiple of
            # 2**-k <= width / 64 where that is shorter, the result is under
            # 5 / 32 of the width: tighter than bisection's quarter, and
            # QIR's quadratic convergence makes the extra precision cheap.
            var k = _bits(width.value.denominator()) - _bits(width.value.numerator()) + 7
            var low = _dyadic_outward(oracle, enclosure.low, k, False, interval, budget)
            var high = _dyadic_outward(oracle, enclosure.high, k, True, interval, budget)
            return RootInterval(QQ(low^), QQ(high^))
        root = RootInterval(QQ(enclosure.low.copy()), QQ(enclosure.low.copy()))
    var low: Rational
    var high: Rational
    try:
        var quarter = width.value / Rational(4)
        low = root.low.value - quarter
        high = root.low.value + quarter
    except e:
        raise _error(e)
    var lower = QQ(low^)
    var upper = QQ(high^)
    if _compare(lower, interval.low, budget) < 0:
        lower = interval.low.copy()
    if _compare(upper, interval.high, budget) > 0:
        upper = interval.high.copy()
    return RootInterval(lower^, upper^)


def _tighten_guarded(oracle: _RootOracle, mut enclosure: _Enclosure, width: QQ, mut budget: Budget) raises PolyError:
    """Shrink the enclosure below width / 8 by QIR, bisecting while it fails.

    A failed QIR step, typical next to a close root outside the enclosure,
    costs two sign evaluations; after one, bisection runs for a number of
    steps that doubles with each consecutive failure before QIR is tried again.
    """
    var target: Rational
    try:
        target = width.value / Rational(8)
    except e:
        raise _error(e)
    var wait = 0
    var backoff = 1
    while not enclosure.exact and _wider(enclosure, target):
        budget.charge(1)
        if wait == 0:
            # Never subdivide below the target: 2**grid parts of the enclosure.
            enclosure.grid = max(2, min(enclosure.grid, _ratio_bits(enclosure, target) + 1))
            if _qir_step(oracle, enclosure, budget):
                enclosure.grid = min(enclosure.grid * 2, 1 << 20)
                backoff = 1
                continue
            enclosure.grid = max(2, enclosure.grid // 2)
            wait = backoff
            backoff = min(backoff * 2, 64)
            if enclosure.exact or not _wider(enclosure, target):
                break
        wait -= 1
        var middle: Rational
        try:
            middle = (enclosure.low + enclosure.high) / Rational(2)
        except e:
            raise _error(e)
        _ = _side(oracle, enclosure, QQ(middle^), budget)


def _dyadic_outward(oracle: _RootOracle, x: Rational, k: Int, up: Bool, interval: RootInterval, mut budget: Budget) raises PolyError -> Rational:
    """x rounded down, or up, to a multiple of 2**-k when that is shorter, inside interval and not a root."""
    budget.charge(1)
    var rounded: Rational
    try:
        var scaled = x.numerator() << k
        var q = x.denominator()
        var n = scaled // q
        if n * q > scaled:
            n = n - Integer(1)
        if up and n * q < scaled:
            n = n + Integer(1)
        rounded = Rational(n, Integer(1) << k)
        if _bits(rounded.numerator()) + _bits(rounded.denominator()) >= _bits(x.numerator()) + _bits(x.denominator()):
            return x.copy()
        if _order(rounded, interval.low.value) <= 0 or _order(interval.high.value, rounded) <= 0:
            return x.copy()
    except e:
        raise _error(e)
    if _check_sign(oracle, rounded, budget) == 0:
        return x.copy()
    return rounded^


def _split_roots(c: List[Integer], oracle: _RootOracle, interval: RootInterval, limit: Int, mut budget: Budget) raises PolyError -> List[RootInterval]:
    """Up to `limit` roots of the squarefree c in an open interval with nonroot ends.

    Parts are split at midpoints until each has Descartes bound 0 or 1, which
    settles for a squarefree polynomial; the bound is exact at 0 and 1, and
    split points that are roots count as roots. Returns bound-1 parts and
    exact points, stopping once `limit` are found.
    """
    var found = List[RootInterval]()
    var pending = List[RootInterval]()
    pending.append(interval.copy())
    while len(pending) and len(found) < limit:
        var part = pending.pop()
        budget.charge(UInt64(len(c) * len(c) // 2 + 1))
        var bound: Int
        try:
            bound = _interval_bound(c, part.low.value, part.high.value, 2)
        except e:
            raise _error(e)
        if bound == 1:
            found.append(part^)
        elif bound > 1:
            var middle: Rational
            try:
                middle = (part.low.value + part.high.value) / Rational(2)
            except e:
                raise _error(e)
            var m = QQ(middle^)
            if _check_sign(oracle, m.value, budget) == 0:
                found.append(RootInterval(m.copy(), m.copy()))
            pending.append(RootInterval(m.copy(), part.high.copy()))
            pending.append(RootInterval(part.low.copy(), m^))
    return found^


def _rational_bound(f: DensePoly[QQ], low: QQ, high: QQ, cap: Int, mut budget: Budget) raises PolyError -> Int:
    """Descartes' bound, capped, for the roots of f in (low, high), counted with multiplicity."""
    var c = _integer_list(_primitive(_integer_dense(f, budget), budget))
    budget.charge(UInt64(len(c) * len(c) // 2 + 1))
    try:
        return _interval_bound(c, low.value, high.value, cap)
    except e:
        raise _error(e)


def _one_root(f: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> Bool:
    """Whether f has a single root in (low, high): opposite nonzero end signs and Descartes bound 1.

    The bound counts multiplicity, so it also excludes a repeated root inside.
    """
    if _compare(low, high, budget) >= 0:
        return False
    var left = _sign_at(f, low, budget)
    var right = _sign_at(f, high, budget)
    return left != 0 and right == -left and _rational_bound(f, low, high, 2, budget) == 1


def _unique_root(f: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> Bool:
    """Whether f has exactly one distinct real root in (low, high), whose ends must be nonroots.

    Descartes' bound of f decides 0 and 1 directly; otherwise the squarefree
    part is subdivided until the count settles.
    """
    if _compare(low, high, budget) >= 0 or _sign_at(f, low, budget) == 0 or _sign_at(f, high, budget) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "root count requires ordered nonroot endpoints")
    var bound = _rational_bound(f, low, high, 2, budget)
    if bound <= 1:
        return bound == 1
    var c = _integer_list(_squarefree(f, budget)[0])
    var oracle = _RootOracle(c.copy(), List[Integer]())
    return len(_split_roots(c, oracle, RootInterval(low.copy(), high.copy()), 2, budget)) == 1


def _narrow_to_one(f: DensePoly[QQ], mut low: QQ, mut high: QQ, mut budget: Budget) raises PolyError:
    """Bisect a squarefree f's single sign-changing root until its Descartes bound is 1.

    The ends keep opposite signs; a midpoint root is enclosed by the two
    quarter points, which then have opposite signs.
    """
    while _rational_bound(f, low, high, 2, budget) != 1:
        budget.charge(1)
        var middle = _midpoint(low, high, budget)
        var at_middle = _sign_at(f, middle, budget)
        if at_middle == 0:
            low = _midpoint(low, middle, budget)
            high = _midpoint(middle, high, budget)
        elif at_middle == _sign_at(f, low, budget):
            low = middle^
        else:
            high = middle^


@fieldwise_init
struct AlgebraicDescartesCertificate(Copyable, Movable, Writable):
    """A complete Descartes isolation of every distinct real root of a rational polynomial.

    The squarefree part and cofactor are recomputed by the checker; every leaf
    of the certificate is resolved, so its brackets and exact roots are all
    the real roots, each once.
    """
    var squarefree: DensePoly[ZZ]
    """The primitive positive-leading squarefree part g = f / gcd(f, f')."""
    var cofactor: DensePoly[ZZ]
    """The primitive gcd(f, f'), one for a squarefree input."""
    var bound_exponent: Int
    """The k with every root of g in (-2^k, 2^k)."""
    var certificate: DescartesCertificate
    """The resolved leaves partitioning (-2^k, 2^k), with brackets and exact roots."""

    def write_to[W: Writer](self, mut writer: W):
        """Write the squarefree part, the bound exponent, the leaf count and every root.

        Parameters:
            W: The writer type.

        Args:
            writer: The destination.
        """
        writer.write("AlgebraicDescartesCertificate(squarefree=", self.squarefree, ", bound_exponent=", self.bound_exponent, ", leaves=", len(self.certificate.leaves), ", exact_roots=[")
        for i in range(len(self.certificate.exact_roots)):
            writer.write(", " if i else "", self.certificate.exact_roots[i])
        writer.write("], brackets=[")
        for i in range(len(self.certificate.brackets)):
            writer.write(", " if i else "", self.certificate.brackets[i])
        writer.write("])")


def _algebraic_descartes(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> AlgebraicDescartesCertificate:
    # Isolation to completion; brackets need no refinement beyond width 1/2.
    var result = _descartes_isolate(f, NumericOptions(absolute_bits=1, max_subdivisions=1 << 60), budget)
    if not result.complete:
        raise PolyError(PolyErrorKind.LIMIT, "real root isolation did not complete")
    return AlgebraicDescartesCertificate(result.squarefree.copy(), result.cofactor.copy(), result.bound_exponent, result.certificate.copy())


def _count_in(isolation: AlgebraicDescartesCertificate, a: QQ, b: QQ, mut budget: Budget) raises PolyError -> Int:
    """The roots of a complete isolation strictly between nonroot a < b.

    Exact roots and point brackets compare directly. An open bracket holds one
    root with a sign change of g across it, so where it contains a (or b) the
    root lies above a exactly when g has the same sign at a as at the
    bracket's lower end (below b exactly when g at b matches its upper end).
    """
    var oracle = _RootOracle(_integer_list(isolation.squarefree), List[Integer]())
    var count = 0
    for r in isolation.certificate.exact_roots:
        budget.charge(1)
        if _compare(a, r, budget) < 0 and _compare(r, b, budget) < 0:
            count += 1
    for bracket in isolation.certificate.brackets:
        budget.charge(1)
        if _compare(bracket.high, a, budget) <= 0 or _compare(b, bracket.low, budget) <= 0:
            continue
        if bracket.low == bracket.high:
            count += 1
            continue
        var above = _compare(a, bracket.low, budget) <= 0 or _check_sign(oracle, a.value, budget) == _check_sign(oracle, bracket.low.value, budget)
        var below = _compare(bracket.high, b, budget) <= 0 or _check_sign(oracle, b.value, budget) == _check_sign(oracle, bracket.high.value, budget)
        if above and below:
            count += 1
    return count


@fieldwise_init
struct DescartesCountCertificate(Copyable, Movable, Writable):
    """A complete Descartes isolation proving a distinct-root count on an open interval.

    The isolation brackets every real root once; the count is the exact roots
    strictly between the endpoints and the brackets inside, a bracket holding
    an endpoint being decided by the squarefree part's sign there.
    """
    var isolation: AlgebraicDescartesCertificate
    """The complete isolation of every distinct real root."""
    var a: QQ
    """The lower endpoint, where the input does not vanish."""
    var b: QQ
    """The upper endpoint, where the input does not vanish."""
    var count: Int
    """The claimed number of distinct roots strictly between the endpoints."""

    def write_to[W: Writer](self, mut writer: W):
        """Write the endpoints, the count and the isolation.

        Parameters:
            W: The writer type.

        Args:
            writer: The destination.
        """
        writer.write("DescartesCountCertificate(a=", self.a, ", b=", self.b, ", count=", self.count, ", isolation=", self.isolation, ")")


def descartes_count(f: DensePoly[QQ], a: QQ, b: QQ, mut budget: Budget) raises PolyError -> DescartesCountCertificate:
    """Count the distinct real roots of f strictly between nonroot endpoints by Descartes' rule.

    Args:
        f: The nonzero polynomial.
        a: The lower endpoint, not a root of f.
        b: The upper endpoint, above a and not a root of f.
        budget: The shared work allowance.

    Returns:
        The count with a complete Descartes isolation that check_descartes_count verifies.

    Raises:
        PolyError: DOMAIN for a zero polynomial or unordered or root endpoints, LIMIT or CAPACITY on work limits.
    """
    if f.is_zero() or _compare(a, b, budget) >= 0 or _sign_at(f, a, budget) == 0 or _sign_at(f, b, budget) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "root count requires ordered nonroot endpoints")
    var isolation = _algebraic_descartes(f, budget)
    var count = _count_in(isolation, a, b, budget)
    return DescartesCountCertificate(isolation^, a.copy(), b.copy(), count)


def check_descartes_count(f: DensePoly[QQ], certificate: DescartesCountCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check a Descartes root count without searching.

    Verifies the isolation with check_algebraic_descartes, requires ordered
    endpoints where f does not vanish, and recounts.

    Args:
        f: The polynomial whose distinct roots are counted.
        certificate: The isolation, endpoints and count.
        budget: The shared work allowance.

    Returns:
        Whether the evidence establishes the count on the open interval.

    Raises:
        PolyError: DOMAIN for malformed polynomial data, LIMIT or CAPACITY on work limits.
    """
    ref a = certificate.a
    ref b = certificate.b
    if f.is_zero() or certificate.count < 0 or _compare(a, b, budget) >= 0 or _sign_at(f, a, budget) == 0 or _sign_at(f, b, budget) == 0:
        return False
    if not check_algebraic_descartes(f, certificate.isolation, budget):
        return False
    return _count_in(certificate.isolation, a, b, budget) == certificate.count


def check_algebraic_descartes(f: DensePoly[QQ], certificate: AlgebraicDescartesCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check a complete Descartes isolation of the real roots of a rational polynomial.

    Recomputes the squarefree part and cofactor, checks the interval against
    the root bound, verifies the certificate without searching, and requires
    every leaf resolved.

    Args:
        f: The nonzero polynomial.
        certificate: The proposed complete isolation.
        budget: The shared work allowance.

    Returns:
        Whether the certificate isolates every distinct real root of f.

    Raises:
        PolyError: DOMAIN for malformed polynomial data, LIMIT or CAPACITY on work limits.
    """
    if f.degree() < 0 or not _check_header(f, certificate.squarefree, certificate.cofactor, certificate.bound_exponent, budget):
        return False
    var outer = _outer(certificate.squarefree, certificate.bound_exponent)
    if certificate.certificate.low != outer[0] or certificate.certificate.high != outer[1]:
        return False
    if not check_descartes(certificate.squarefree, certificate.certificate, budget):
        return False
    for leaf in certificate.certificate.leaves:
        if not leaf.resolved:
            return False
    return True


def _witness_intervals(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[RootInterval]:
    """Increasing intervals, each holding one real root of p with Descartes bound 1 and a sign change.

    For p without rational roots: a bracket inside a leaf of bound 1 has bound
    at most 1 by subadditivity, and a sign change makes it exactly 1; brackets
    of wider leaves are narrowed on their dyadic grid until the bound is 1.
    """
    var result = _descartes_isolate(p, NumericOptions(absolute_bits=1, max_subdivisions=1 << 60), budget)
    if not result.complete:
        raise PolyError(PolyErrorKind.LIMIT, "real root isolation did not complete")
    ref certificate = result.certificate
    var g = _rational_from(result.squarefree, budget)
    var intervals = List[RootInterval]()
    var next = 0
    for leaf in certificate.leaves:
        budget.charge(1)
        for _ in range(leaf.variations):
            ref bracket = certificate.brackets[next]
            next += 1
            if bracket.low == bracket.high:
                raise PolyError(PolyErrorKind.DOMAIN, "witness intervals need a polynomial without rational roots")
            var low = bracket.low.copy()
            var high = bracket.high.copy()
            if leaf.variations > 1:
                _narrow_to_one(g, low, high, budget)
            _allocation[RootInterval](1, budget)
            intervals.append(RootInterval(low^, high^))
    if len(certificate.exact_roots):
        raise PolyError(PolyErrorKind.DOMAIN, "witness intervals need a polynomial without rational roots")
    return intervals^


def _strip_factor(f: DensePoly[QQ], p: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    """F with every factor p divided out exactly: f = p^k q with p not dividing q."""
    var q = f.copy()
    if p.degree() < 1:
        return q^
    while q.degree() >= p.degree():
        budget.charge(1)
        var division = divmod(q, p, budget)
        if not division.remainder.is_zero():
            break
        q = division.quotient.copy()
    return q^
