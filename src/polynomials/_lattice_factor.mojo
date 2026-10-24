"""Van Hoeij trace-knapsack recombination with exact lattice reduction.

Fresh implementation from Mark van Hoeij (2002), sections 2.1 and 2.2,
DOI 10.1006/jnth.2001.2763. Scaled traces handle nonmonic inputs.
"""
from apn_mojo import div_rem_euclid, Integer, inverse_mod
from .coefficients import ZZ, QQ, ModInt, NoContext, _work, _pow
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _product, _sum, _charge_size
from ._coefficient_ops import _exact, _quotient, _copy as _zcopy
from ._dense_ops import _copy, _constant
from ._hensel import _HenselState, _hensel_start, _hensel_extend, _truncate, _primitive, _trial_divide, _divmod_monic
from ._dense_ops import _derivative
from .primes import _symmetric
from ._lattice import _lll, _lll_knapsack, _q, _exact_norms, _modular_norms
from ._root_ops import _compare


def _greater(a: ZZ, b: ZZ, mut budget: Budget) raises PolyError -> Bool:
    _work(a.size(), b.size(), False, budget)
    return a.value > b.value


def _traces(f: DensePoly[ZZ], count: Int, modulus: ZZ, mut budget: Budget) raises PolyError -> List[ZZ]:
    # Newton sums of the roots of the monic modular factor.
    _allocation[ZZ](_sum(count, 1), budget)
    var sums: List[ZZ] = [ZZ.zero(NoContext())]
    for k in range(1, _sum(count, 1)):
        budget.charge(1)
        var total = ZZ.zero(NoContext())
        for j in range(1, min(k, _sum(f.degree(), 1))):
            budget.charge(1)
            total = total.add(f.coefficients[f.degree() - j].mul(sums[k - j], budget), budget)
            total = ZZ(_symmetric(total.value, modulus.value, budget))
        if k <= f.degree():
            total = total.add(ZZ(Integer(k)).mul(f.coefficients[f.degree() - k], budget), budget)
        sums.append(ZZ(_symmetric(total.neg(budget).value, modulus.value, budget)))
    return sums^


comptime _FILTER_PRIME = UInt64(2305843009213693951)
"""2^61 - 1."""


def _modular_partition(basis: List[List[ZZ]], count: Int, mut budget: Budget) raises PolyError -> List[List[QQ]]:
    """The groups of the projection's reduced row echelon form modulo
    2^61 - 1 when it is a 0/1 partition, else no groups.

    Every true factor's indicator is an integer combination of the rows, so
    modulo the prime it lies in the span of the groups; a 0/1 vector there is
    a union of groups. Hence when every group's product reconstructs a true
    factor, the groups are exactly the irreducible factors.
    """
    var p = _FILTER_PRIME
    var rows = len(basis)
    var modulus = Integer(Int(p))
    var m = List[List[UInt64]](capacity=rows)
    for i in range(rows):
        budget.charge(1)
        var row = List[UInt64](capacity=count)
        for j in range(count):
            budget.charge(1)
            try:
                row.append(UInt64(Int(div_rem_euclid(basis[i][j].value, modulus)[1])))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
        m.append(row^)
    var rank = 0
    var pp = UInt128(p)
    for column in range(count):
        budget.charge(1)
        var pivot = rank
        while pivot < rows and m[pivot][column] == 0:
            pivot += 1
        if pivot == rows:
            continue
        if pivot != rank:
            var t = m[pivot].copy()
            m[pivot] = m[rank].copy()
            m[rank] = t^
        # inverse of the pivot by Fermat
        var inverse = UInt64(1)
        var base = m[rank][column]
        var e = p - 2
        while e:
            if e & 1:
                inverse = UInt64((UInt128(inverse) * UInt128(base)) % pp)
            base = UInt64((UInt128(base) * UInt128(base)) % pp)
            e >>= 1
        _charge_size(count, budget)
        for j in range(column, count):
            m[rank][j] = UInt64((UInt128(m[rank][j]) * UInt128(inverse)) % pp)
        for i in range(rows):
            if i == rank or m[i][column] == 0:
                continue
            var factor = m[i][column]
            _charge_size(count, budget)
            for j in range(column, count):
                var y = UInt64((UInt128(factor) * UInt128(m[rank][j])) % pp)
                var z = m[i][j]
                m[i][j] = z - y if z >= y else z + (p - y)
        rank += 1
        if rank == rows:
            break
    var none = List[List[QQ]]()
    for j in range(count):
        budget.charge(1)
        var occurrences = 0
        for i in range(rank):
            if m[i][j] == 0:
                continue
            if m[i][j] != 1:
                return none^
            occurrences += 1
        if occurrences != 1:
            return none^
    _allocation[List[QQ]](rank, budget)
    _allocation[QQ](_product(rank, count), budget)
    var groups = List[List[QQ]](capacity=rank)
    for i in range(rank):
        budget.charge(1)
        var row = List[QQ](capacity=count)
        for j in range(count):
            row.append(QQ.one(NoContext()) if m[i][j] == 1 else QQ.zero(NoContext()))
        groups.append(row^)
    return groups^


def _shifted_upper(value: Integer, shift: Int) raises PolyError -> Integer:
    """An upper bound on value * 2^shift for value >= 0."""
    try:
        if shift >= 0:
            return value << shift
        return (value >> (-shift)) + Integer(1)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _cld_bound(f: DensePoly[ZZ], j: Int, mut budget: Budget) raises PolyError -> ZZ:
    """A bound on |coefficient j of f g'/g| over every factor g of f.

    f g'/g sums f/(x - a) over the roots a of g, and coefficient j of
    f/(x - a) equals sum_(i>j) f_i a^(i-j-1) = -sum_(i<=j) f_i a^(i-j-1). The
    first sum grows with |a| and the second shrinks, so for any t = 2^e the
    larger of the two at t bounds every root; deg(f) times it bounds the sum.
    The exponent e is chosen near the crossing from bit lengths.
    """
    var n = f.degree()
    var bits = List[Int](capacity=n + 1)
    for i in range(n + 1):
        bits.append(f.coefficients[i].value.magnitude_bit_length())
    var hi = 0
    for b in bits:
        hi = max(hi, b + 64)
    var lo = -hi
    # Smallest e whose estimated first sum reaches the second.
    while lo < hi:
        budget.charge(1)
        var mid = (lo + hi) // 2 if (lo + hi) >= 0 else -((-(lo + hi) + 1) // 2)
        var la = -1 << 30
        var lb = -1 << 30
        for i in range(n + 1):
            if bits[i] == 0:
                continue
            var size = bits[i] + mid * (i - j - 1)
            if i > j:
                la = max(la, size)
            else:
                lb = max(lb, size)
        if la >= lb:
            hi = mid
        else:
            lo = mid + 1
    var best = ZZ.zero(NoContext())
    var found = False
    for e in [lo - 1, lo]:
        _work(n, 1, False, budget)
        var first = Integer(0)
        var second = Integer(0)
        for i in range(n + 1):
            var magnitude = f.coefficients[i].value.copy()
            if magnitude.sign() == 0:
                continue
            if magnitude.sign() < 0:
                magnitude = -magnitude
            var term = _shifted_upper(magnitude, e * (i - j - 1))
            try:
                if i > j:
                    first = first + term
                else:
                    second = second + term
            except err:
                raise PolyError(PolyErrorKind.CAPACITY, String(err))
        var larger = ZZ(first.copy() if first > second else second.copy())
        if not found or _greater(best, larger, budget):
            best = larger^
            found = True
    return best.mul(ZZ(Integer(n)), budget)


def _may_divide(f: DensePoly[ZZ], g: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Bool:
    """False only when g cannot divide f: an exact divisor's constant term
    divides f's, which rejects most false candidates before a division."""
    if f.coefficients[0].is_zero():
        return True
    if g.coefficients[0].is_zero():
        return False
    return Bool(_quotient(f.coefficients[0], g.coefficients[0], budget))


comptime _CLD_SERIES = 64
"""Coefficients from either end of g' (f/g) up to which truncated series
replace the division of f by g."""


def _head(p: DensePoly[ZZ], count: Int, mut budget: Budget) raises PolyError -> List[ZZ]:
    """The first count coefficients of p, fewer when p is shorter."""
    var length = min(count, len(p.coefficients))
    _allocation[ZZ](length, budget)
    var out = List[ZZ](capacity=length)
    for i in range(length):
        out.append(p.coefficients[i].copy())
    return out^


def _series_coefficients(a: List[ZZ], inverse: List[ZZ], count: Int, modulus: ZZ, mut budget: Budget) raises PolyError -> List[ZZ]:
    """The first count coefficients of a * inverse modulo the modulus."""
    _allocation[ZZ](count, budget)
    var out = List[ZZ](capacity=count)
    for t in range(count):
        budget.charge(1)
        var total = ZZ.zero(NoContext())
        for u in range(min(t + 1, len(a))):
            if t - u < len(inverse) and not a[u].is_zero():
                total = total.add(a[u].mul(inverse[t - u], budget), budget)
        out.append(ZZ(_symmetric(total.value, modulus.value, budget)))
    return out^


def _series_inverse_of(a: List[ZZ], head: ZZ, count: Int, modulus: ZZ, mut budget: Budget) raises PolyError -> List[ZZ]:
    """The first count coefficients of 1/a modulo the modulus, with head = 1/a_0."""
    _allocation[ZZ](count, budget)
    var out = List[ZZ](capacity=count)
    out.append(_zcopy(head, budget))
    for k in range(1, count):
        budget.charge(1)
        var total = ZZ.zero(NoContext())
        for t in range(1, min(k, len(a) - 1) + 1):
            if not a[t].is_zero():
                total = total.add(a[t].mul(out[k - t], budget), budget)
        out.append(ZZ(_symmetric(total.mul(head, budget).neg(budget).value, modulus.value, budget)))
    return out^


def _cld_values(f: DensePoly[ZZ], g: DensePoly[ZZ], indices: List[Int], modulus: ZZ, mut budget: Budget) raises PolyError -> Optional[List[ZZ]]:
    """Coefficients `indices` of g' (f/g), symmetric modulo the modulus, for
    a monic g that divides f modulo it, as lifted factors do, or no value.

    With f = q g exactly modulo the modulus, q = f / g as power series from
    the bottom when g(0) is a unit, and rev(q) = rev(f) / rev(g) from the
    top, so the coefficients near either end need only that many terms.
    """
    var n = f.degree()
    var d = g.degree()
    var bottom = 0
    var top = 0
    for j in indices:
        if 2 * j <= n - 1:
            bottom = max(bottom, j + 1)
        else:
            top = max(top, n - j)
    if bottom > _CLD_SERIES or top > _CLD_SERIES:
        return None
    var derivative = _derivative(g, budget)
    var low = List[ZZ]()
    if bottom > 0:
        var head: Integer
        try:
            head = inverse_mod(g.coefficients[0].value, modulus.value)
        except:
            return None
        var inverse = _series_inverse_of(_head(g, bottom, budget), ZZ(head^), bottom, modulus, budget)
        var quotient = _series_coefficients(_head(f, bottom, budget), inverse, bottom, modulus, budget)
        low = _series_coefficients(_head(derivative, bottom, budget), quotient, bottom, modulus, budget)
    var high = List[ZZ]()
    if top > 0:
        _allocation[ZZ](_sum(top, d) + 1, budget)
        var rf = List[ZZ](capacity=top)
        for t in range(min(top, n + 1)):
            rf.append(f.coefficients[n - t].copy())
        var rg = List[ZZ](capacity=d + 1)
        for t in range(d + 1):
            rg.append(g.coefficients[d - t].copy())
        var rd = List[ZZ](capacity=d)
        for t in range(d):
            rd.append(derivative.coefficients[d - 1 - t].copy() if d - 1 - t < len(derivative.coefficients) else ZZ.zero(NoContext()))
        var inverse = _series_inverse_of(rg, ZZ.one(NoContext()), top, modulus, budget)
        var quotient = _series_coefficients(rf, inverse, top, modulus, budget)
        high = _series_coefficients(rd, quotient, top, modulus, budget)
    _allocation[ZZ](len(indices), budget)
    var values = List[ZZ](capacity=len(indices))
    for j in indices:
        budget.charge(1)
        values.append(low[j].copy() if 2 * j <= n - 1 else high[n - 1 - j].copy())
    return values^


def _reconstruct(f: DensePoly[ZZ], partition: List[List[QQ]], leading: ZZ, state: _HenselState, mut budget: Budget) raises PolyError -> Optional[List[DensePoly[ZZ]]]:
    """The factors lc(f) prod_(j in group) g_j modulo the lifted modulus,
    when every one divides f exactly with its group's degree."""
    var remaining = _copy(f, budget)
    var factors = List[DensePoly[ZZ]]()
    for row in partition:
        budget.charge(1)
        var candidate = _constant(f, leading, budget)
        var degree = 0
        for j in range(len(row)):
            budget.charge(1)
            if not row[j].is_zero():
                degree = _sum(degree, state.lifted[j].degree())
                candidate = _truncate(candidate.mul(state.lifted[j], budget), state.modulus, budget)
        candidate = _primitive(candidate, budget)
        if candidate.degree() != degree or degree < 1 or not _may_divide(remaining, candidate, budget):
            return None
        var quotient = _trial_divide(remaining, candidate, budget)
        if not quotient:
            return None
        remaining = quotient.take()
        _allocation[DensePoly[ZZ]](1, budget)
        factors.append(candidate^)
    if remaining.degree() != 0:
        return None
    return factors^


def _root_bound(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> ZZ:
    """A bound on |lc(f) z| over the roots z of f."""
    var root_bound = ZZ.zero(NoContext())
    for coefficient_index in range(len(f.coefficients)):
        ref coefficient = f.coefficients[coefficient_index]
        budget.charge(1)
        var magnitude = coefficient.neg(budget) if coefficient.value.sign() < 0 else _zcopy(coefficient, budget)
        if _greater(magnitude, root_bound, budget):
            root_bound = magnitude^
    ref lc = f.coefficients[f.degree()]
    var leading = lc.neg(budget) if lc.value.sign() < 0 else _zcopy(lc, budget)
    root_bound = root_bound.add(leading, budget)
    # |lc z| <= |lc| + max |f_i| (Cauchy) bounds every root z; Fujiwara's
    # |z| <= 2 max_i |f_(n-i) / f_n|^(1/i) is usually far smaller.
    var fujiwara = _fujiwara_power(f, budget)
    var scaled = _zcopy(root_bound, budget)
    if fujiwara < 4096:
        scaled = leading.mul(_pow(ZZ(Integer(2)), UInt64(fujiwara), budget), budget)
    if _greater(root_bound, scaled, budget):
        root_bound = scaled^
    return root_bound^


def _padding(count: Int, p: ZZ, mut budget: Budget) raises PolyError -> ZZ:
    """The least power of p with at least count + 20 bits."""
    var padding = _zcopy(p, budget)
    while padding.value.magnitude_bit_length() < _sum(count, 20):
        budget.charge(1)
        padding = padding.mul(p, budget)
    return padding^


def _trace_precision(f: DensePoly[ZZ], count: Int, p: ZZ, mut budget: Budget) raises PolyError -> ZZ:
    """The modulus the first two power sums need, as _van_hoeij_lifted
    computes it."""
    var root_bound = _root_bound(f, budget)
    var power = ZZ.one(NoContext())
    var cut = ZZ.one(NoContext())
    for _ in range(min(f.degree(), 2)):
        budget.charge(1)
        power = power.mul(root_bound, budget)
        var bound = power.mul(ZZ(Integer(_product(2, f.degree()))), budget)
        cut = ZZ.one(NoContext())
        while not _greater(cut, bound, budget):
            budget.charge(1)
            cut = cut.mul(p, budget)
    return cut.mul(_padding(count, p, budget), budget)


def _cld_precision(f: DensePoly[ZZ], count: Int, p: ZZ, mut budget: Budget) raises PolyError -> ZZ:
    """The modulus the first CLD columns, at both ends, need."""
    var largest_cut = ZZ.one(NoContext())
    for index in [0, f.degree() - 1]:
        budget.charge(1)
        var bound = _cld_bound(f, index, budget).mul(ZZ(Integer(2)), budget)
        var cut = ZZ.one(NoContext())
        while not _greater(cut, bound, budget):
            budget.charge(1)
            cut = cut.mul(p, budget)
        if _greater(cut, largest_cut, budget):
            largest_cut = cut^
    return largest_cut.mul(_padding(count, p, budget), budget)


def _uses_traces(f: DensePoly[ZZ], count: Int, p: ZZ, mut budget: Budget) raises PolyError -> Bool:
    """Whether the lattice takes power sums rather than CLD columns: only
    when their first precision is the smaller, which a large leading
    coefficient or large roots rule out."""
    if f.coefficients[f.degree()].value.magnitude_bit_length() >= 1024:
        return False
    return not _greater(_trace_precision(f, count, p, budget), _cld_precision(f, count, p, budget), budget)


def _lattice_precision(f: DensePoly[ZZ], count: Int, p: ZZ, mut budget: Budget) raises PolyError -> Optional[ZZ]:
    """The modulus the first lattice round lifts to, for the columns
    _uses_traces picks."""
    if _uses_traces(f, count, p, budget):
        return _trace_precision(f, count, p, budget)
    return _cld_precision(f, count, p, budget)


def _fujiwara_power(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Int:
    """e with every root of f at most 2^e in absolute value: Fujiwara's
    bound 2 max_i |f_(n-i) / f_n|^(1/i), each ratio below
    2^(bits(f_(n-i)) - bits(f_n) + 1), its i-th root rounded up."""
    var n = f.degree()
    var lead_bits = f.coefficients[n].value.magnitude_bit_length()
    var best = -(1 << 30)
    for i in range(1, n + 1):
        budget.charge(1)
        ref c = f.coefficients[n - i]
        if c.value.sign() == 0:
            continue
        var numerator = c.value.magnitude_bit_length() - lead_bits + 1
        # ceil(numerator / i), also for negative numerators
        var e = -((-numerator) // i)
        best = max(best, e)
    return max(best, 0) + 1


def _extend_to(f: DensePoly[ZZ], images: List[DensePoly[ModInt]], target: ZZ, square: Bool, mut state: _HenselState, mut budget: Budget) raises PolyError:
    """Lift to a multiple of target, a power of p. A modulus equal to it
    serves the columns, and _hensel_extend lifts strictly past its bound,
    so it gets target - 1. With square, a lift goes at least to the square
    of the modulus, which costs a trace lift about as much as any target
    below it, and later batches use the headroom."""
    if not _greater(target, state.modulus, budget):
        return
    var goal = target.sub(ZZ.one(NoContext()), budget)
    var squared = state.modulus.mul(state.modulus, budget)
    if square and _greater(squared, goal, budget):
        goal = squared.sub(ZZ.one(NoContext()), budget)
    _hensel_extend(f, images, goal, state, budget)


def _van_hoeij(f: DensePoly[ZZ], images: List[DensePoly[ModInt]], reconstruction_bound: ZZ, mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    var state = _hensel_start(images, budget)
    return _van_hoeij_lifted(f, images, reconstruction_bound, state, budget)


def _van_hoeij_lifted(f: DensePoly[ZZ], images: List[DensePoly[ModInt]], reconstruction_bound: ZZ, mut state: _HenselState, mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    var count = len(images)
    if count < 2 or f.degree() < 2:
        raise PolyError(PolyErrorKind.DOMAIN, "lattice recombination requires at least two modular factors")
    var p = ZZ(Integer(Int(images[0].ring.value)))
    ref lc = f.coefficients[f.degree()]
    var leading = lc.neg(budget) if lc.value.sign() < 0 else _zcopy(lc, budget)
    var root_bound = _root_bound(f, budget)
    # Power sums of the roots need precision about (|lc|+max|f_i|)^k for
    # the k-th one. With a large leading coefficient or large roots, that
    # drives lifting far beyond the reconstruction bound, which lc(f) does
    # not enter, and coefficients of the logarithmic derivative (CLD),
    # bounded near the size of f, are cheaper.
    var monic = _uses_traces(f, count, p, budget)
    # Gradual feeding: the reduced basis rows persist, and each round
    # appends one data column at its own precision cut_j * padding, so every
    # new entry is about the padding in size. Rows whose exact Gram-Schmidt
    # norm exceeds the bound for true indicators are dropped, which never
    # loses one, and a partition can end the search before every column is in.
    _allocation[List[ZZ]](count, budget)
    var rows = List[List[ZZ]](length=count, fill=List[ZZ](length=count, fill=ZZ.zero(NoContext())))
    for i in range(count):
        rows[i][i] = ZZ.one(NoContext())
    var done = 0
    var traces_count = min(f.degree(), 2)
    # CLD indices from both ends inward: 0, n-1, 1, n-2, ...
    var order = List[Int]()
    var low_index = 0
    var high_index = f.degree() - 1
    while low_index <= high_index:
        order.append(low_index)
        if high_index != low_index:
            order.append(high_index)
        low_index += 1
        high_index -= 1
    # Each column must add clearly more than about one bit per modular factor,
    # or LLL cannot push the false indicators above the bound.
    var padding = _zcopy(p, budget)
    while padding.value.magnitude_bit_length() < _sum(count, 20):
        budget.charge(1)
        padding = padding.mul(p, budget)
    # Columns wait in pending, computed a batch of traces at a time.
    var pending = List[List[ZZ]]()
    var pending_padding = List[ZZ]()
    var fed = 0
    while True:
        budget.charge(1)
        if len(pending) == 0:
            if done >= traces_count:
                if traces_count < f.degree():
                    traces_count = min(f.degree(), _product(2, traces_count))
                else:
                    # Every trace is in; feed them again with more precision.
                    padding = padding.mul(p, budget)
                    done = 0
            budget.charge(1)
            var fresh = traces_count - done
            _allocation[List[ZZ]](count, budget)
            var columns = List[List[ZZ]](capacity=count)
            if monic:
                # For a rational factor, lc(f)^k * Tr_k is an integer with
                # absolute value at most degree(f)*(|lc(f)|+max |f_i|)^k.
                _allocation[ZZ](traces_count, budget)
                var cuts = List[ZZ]()
                var power = ZZ.one(NoContext())
                var largest_cut = ZZ.one(NoContext())
                for _ in range(traces_count):
                    budget.charge(1)
                    power = power.mul(root_bound, budget)
                    var bound = power.mul(ZZ(Integer(_product(2, f.degree()))), budget)
                    var cut = ZZ.one(NoContext())
                    while not _greater(cut, bound, budget):
                        budget.charge(1)
                        cut = cut.mul(p, budget)
                    largest_cut = _zcopy(cut, budget)
                    cuts.append(cut^)
                var target = largest_cut.mul(padding, budget)
                _extend_to(f, images, target, True, state, budget)
                ref lifted = state.lifted
                ref modulus = state.modulus
                # columns[i][j] = (lc^(j+1) Tr_(j+1)(g_i) mod cut_j*padding - low) / cut_j.
                for i in range(count):
                    budget.charge(1)
                    var sums = _traces(lifted[i], traces_count, modulus, budget)
                    var scale = ZZ.one(NoContext())
                    var row = List[ZZ](capacity=fresh)
                    for j in range(traces_count):
                        budget.charge(1)
                        scale = ZZ(_symmetric(scale.mul(leading, budget).value, modulus.value, budget))
                        if j < done:
                            continue
                        var precision = cuts[j].mul(padding, budget)
                        var trace = ZZ(_symmetric(sums[j + 1].mul(scale, budget).value, precision.value, budget))
                        var low = ZZ(_symmetric(trace.value, cuts[j].value, budget))
                        row.append(_exact(trace.sub(low, budget), cuts[j], budget))
                    columns.append(row^)
            else:
                # Data columns are coefficients of the logarithmic derivative
                # f g'/g = g' (f/g), additive over modular factors and bounded by
                # about the size of f's coefficients for every index (CLD).
                _allocation[ZZ](traces_count, budget)
                var cuts = List[ZZ]()
                var largest_cut = ZZ.one(NoContext())
                for k in range(traces_count):
                    budget.charge(1)
                    var bound = _cld_bound(f, order[k], budget).mul(ZZ(Integer(2)), budget)
                    var cut = ZZ.one(NoContext())
                    while not _greater(cut, bound, budget):
                        budget.charge(1)
                        cut = cut.mul(p, budget)
                    if _greater(cut, largest_cut, budget):
                        largest_cut = _zcopy(cut, budget)
                    cuts.append(cut^)
                var target = largest_cut.mul(padding, budget)
                _extend_to(f, images, target, False, state, budget)
                ref lifted = state.lifted
                ref modulus = state.modulus
                # columns[i][k] = (coefficient order[k] of g_i' f/g_i mod cut_k*padding - low) / cut_k.
                var needed = List[Int](capacity=fresh)
                for k in range(done, traces_count):
                    needed.append(order[k])
                for i in range(count):
                    budget.charge(1)
                    var values = _cld_values(f, lifted[i], needed, modulus, budget)
                    if not values:
                        var quotient = _divmod_monic(f, lifted[i], modulus, budget)[0].copy()
                        var cld = _truncate(_derivative(lifted[i], budget).mul(quotient, budget), modulus, budget)
                        _allocation[ZZ](fresh, budget)
                        var all = List[ZZ](capacity=fresh)
                        for index in needed:
                            all.append(cld.coefficients[index].copy() if index < len(cld.coefficients) else ZZ.zero(NoContext()))
                        values = all^
                    var row = List[ZZ](capacity=fresh)
                    for k in range(done, traces_count):
                        budget.charge(1)
                        var value = values.value()[k - done].copy()
                        var precision = cuts[k].mul(padding, budget)
                        var data = ZZ(_symmetric(value.value, precision.value, budget))
                        var low = ZZ(_symmetric(data.value, cuts[k].value, budget))
                        row.append(_exact(data.sub(low, budget), cuts[k], budget))
                    columns.append(row^)
            for t in range(fresh):
                budget.charge(1)
                var column = List[ZZ](capacity=count)
                for i in range(count):
                    column.append(columns[i][t].copy())
                pending.append(column^)
                pending_padding.append(_zcopy(padding, budget))
            done = traces_count
        var column = pending.pop(0)
        var column_padding = pending_padding.pop(0)
        for a in range(len(rows)):
            budget.charge(1)
            var total = ZZ.zero(NoContext())
            for i in range(count):
                if not rows[a][i].is_zero():
                    total = total.add(rows[a][i].mul(column[i], budget), budget)
            rows[a].append(ZZ(_symmetric(total.value, column_padding.value, budget)))
        var modulus_row = List[ZZ](length=len(rows[0]), fill=ZZ.zero(NoContext()))
        modulus_row[len(modulus_row) - 1] = column_padding^
        rows.append(modulus_row^)
        fed += 1
        var dimension = len(rows)
        var reduced = _lll_knapsack(rows, budget, False)
        # A true factor's indicator u has |u|^2 <= n, and its entry in each
        # column is (T - sum of lows)/cut with |T| < cut/2 and every low at
        # most cut/2, so at most (n + 1)/2. Dropping only a suffix of larger
        # Gram-Schmidt norms retains every such short vector.
        var entry = (count + 2) // 2
        var weight = _sum(count, _product(fed, _product(entry, entry)))
        var bound_squared = _q(ZZ(Integer(weight)), budget)
        var retained = dimension
        # Removal drops only a suffix, so it can start only if the last row
        # may exceed the bound; keeping rows is always safe, and the exact
        # norms decide every row that is dropped.
        ref estimates = reduced.estimates
        if len(estimates) < dimension or estimates[dimension - 1] > Float64(weight):
            var modular = _modular_norms(reduced.basis, budget)
            var norms = modular.take() if modular else _exact_norms(reduced.basis, budget)
            while retained and _compare(norms[retained - 1], bound_squared, budget) > 0:
                budget.charge(1)
                retained -= 1
        if retained == 0:
            raise PolyError(PolyErrorKind.DOMAIN, "lattice reduction lost every factor indicator")
        _allocation[List[ZZ]](retained, budget)
        var kept = List[List[ZZ]](capacity=retained)
        for a in range(retained):
            budget.charge(1)
            kept.append(reduced.basis[a].copy())
        rows = kept^
        # The indicators span at most count dimensions, so more rows hold
        # vectors beyond them; a partition is tried once at most count rows
        # remain, or once every column is in.
        var exhausted = traces_count >= f.degree() and done >= traces_count and len(pending) == 0
        var partition = List[List[QQ]]()
        if len(rows) <= count or exhausted:
            partition = _modular_partition(rows, count, budget)
        if len(partition) == 1:
            # Every true factor is a union of groups, so with one group f
            # is irreducible, at any precision. As _reconstruct would give
            # it: primitive, with a positive leading coefficient.
            var whole = _primitive(_copy(f, budget), budget)
            if whole.coefficients[whole.degree()].value.sign() < 0:
                whole = whole.scale(ZZ(Integer(-1)), budget)
            var irreducible = List[DensePoly[ZZ]]()
            irreducible.append(whole^)
            return irreducible^
        if partition:
            # Try the lattice's precision first. A divisor of f with the
            # modular image of a group and its full degree is the unique
            # such factor, so a lower precision cannot change the answer.
            var factors = _reconstruct(f, partition, leading, state, budget)
            # Lift for the full bound only once the lattice has merged some
            # factors or holds every column; until then, feed it further.
            var settled = len(partition) < count or exhausted
            if not factors and settled and _greater(reconstruction_bound, state.modulus, budget):
                _hensel_extend(f, images, reconstruction_bound, state, budget)
                factors = _reconstruct(f, partition, leading, state, budget)
            if factors:
                return factors.take()
