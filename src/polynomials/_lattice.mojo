"""Exact rational Gram-Schmidt and integer LLL for van Hoeij recombination.

Fresh implementation of the lattice step in van Hoeij (2002), section 2.1.
All Gram-Schmidt calculations are exact and fraction-free; no floating-point
decisions occur.
"""
from std.bit import count_leading_zeros
from apn_mojo import Rational, Integer, div_rem_floor
from .coefficients import ZZ, QQ, NoContext, _work
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _product, _sum, _charge_size
from ._root_ops import _compare
from ._coefficient_ops import _copy
from ._nmod import _Field, _field, _mulm, _subm, _addm, _to, _from, _invm
from .primes import prime


@fieldwise_init
struct _Orthogonal(Movable):
    var mu: List[List[QQ]]
    var norms: List[QQ]


@fieldwise_init
struct _LLL(Movable):
    var basis: List[List[ZZ]]
    var transform: List[List[ZZ]]
    var mu: List[List[QQ]]
    var norms: List[QQ]
    var estimates: List[Float64]
    """Floating squared Gram-Schmidt norms, when the reduction kept them."""


def _q(a: ZZ, mut budget: Budget) raises PolyError -> QQ:
    _charge_size(a.size(), budget)
    return QQ(Rational(a.value.copy()))


def _copy_matrix(a: List[List[ZZ]], mut budget: Budget) raises PolyError -> List[List[ZZ]]:
    _allocation[List[ZZ]](len(a), budget)
    var result = List[List[ZZ]]()
    for row in a:
        budget.charge(1)
        _allocation[ZZ](len(row), budget)
        var copied = List[ZZ]()
        for value in row:
            _charge_size(value.size(), budget)
            copied.append(value.copy())
        result.append(copied^)
    return result^


def _orthogonal(a: List[List[ZZ]], mut budget: Budget) raises PolyError -> _Orthogonal:
    var rows = len(a)
    var columns = len(a[0])
    _allocation[List[QQ]](rows, budget)
    _allocation[QQ](_product(rows, rows), budget)
    _allocation[QQ](_product(rows, columns), budget)
    _allocation[QQ](rows, budget)
    var mu = List[List[QQ]](length=rows, fill=List[QQ](length=rows, fill=QQ.zero(NoContext())))
    var orthogonal = List[List[QQ]]()
    var norms = List[QQ]()
    _allocation[QQ](rows, budget)
    var inverse_norms = List[QQ](capacity=rows)
    for i in range(rows):
        budget.charge(1)
        var vector = List[QQ]()
        for value in a[i]:
            budget.charge(1)
            vector.append(_q(value, budget))
        for j in range(i):
            budget.charge(1)
            var dot = QQ.zero(NoContext())
            for k in range(columns):
                budget.charge(1)
                dot = dot.add(_q(a[i][k], budget).mul(orthogonal[j][k], budget), budget)
            mu[i][j] = dot.mul(inverse_norms[j], budget)
            for k in range(columns):
                budget.charge(1)
                vector[k] = vector[k].sub(mu[i][j].mul(orthogonal[j][k], budget), budget)
        var norm = QQ.zero(NoContext())
        for value in vector:
            budget.charge(1)
            norm = norm.add(value.mul(value, budget), budget)
        if norm.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "LLL requires independent nonzero rows")
        if i + 1 < rows:
            inverse_norms.append(norm.inverse(budget))
        norms.append(norm^)
        orthogonal.append(vector^)
    return _Orthogonal(mu^, norms^)


def _nearest(q: QQ, mut budget: Budget) raises PolyError -> ZZ:
    _work(q.size(), 1, True, budget)
    try:
        return ZZ(q.value.round())
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _swap_lists[T: Copyable & Deinitable](mut a: List[List[T]], i: Int, j: Int):
    """Exchange two rows by moving their buffers, without copying."""
    var temporary = List[T]()
    swap(temporary, a[i])
    swap(temporary, a[j])
    swap(temporary, a[i])


def _swap_rows(mut a: List[List[ZZ]], i: Int, j: Int):
    var temporary = List[ZZ]()
    swap(temporary, a[i])
    swap(temporary, a[j])
    swap(temporary, a[i])


def _imul(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    _work(a.magnitude_bit_length() // 64 + 1, b.magnitude_bit_length() // 64 + 1, True, budget)
    try:
        return a * b
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _iadd(a: Integer, b: Integer, subtract: Bool, mut budget: Budget) raises PolyError -> Integer:
    _work(a.magnitude_bit_length() // 64 + 1, b.magnitude_bit_length() // 64 + 1, False, budget)
    try:
        return a - b if subtract else a + b
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _ratio(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> QQ:
    _work(a.magnitude_bit_length() // 64 + 1, b.magnitude_bit_length() // 64 + 1, True, budget)
    try:
        return QQ(Rational(a.copy(), b.copy()))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _iexact(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    """The exact quotient a / b of the fraction-free recurrences."""
    _work(a.magnitude_bit_length() // 64 + 1, b.magnitude_bit_length() // 64 + 1, True, budget)
    try:
        return div_rem_floor(a, b)[0].copy()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _rounded(a: Integer, d: Integer, mut budget: Budget) raises PolyError -> Integer:
    """floor(a / d + 1/2) for d > 0."""
    _work(a.magnitude_bit_length() // 64 + 1, d.magnitude_bit_length() // 64 + 1, True, budget)
    try:
        return div_rem_floor(a + a + d, d + d)[0].copy()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _row_dot(a: List[ZZ], b: List[ZZ], mut budget: Budget) raises PolyError -> Integer:
    var total = Integer(0)
    for i in range(len(a)):
        budget.charge(1)
        if a[i].value.sign() == 0 or b[i].value.sign() == 0:
            continue
        total = _iadd(total, _imul(a[i].value, b[i].value, budget), False, budget)
    return total^


def _reduce[carry: Bool](k: Int, l: Int, mut basis: List[List[ZZ]], mut transform: List[List[ZZ]], mut lam: List[List[Integer]], d: List[Integer], mut budget: Budget) raises PolyError:
    budget.charge(1)
    var twice = _iadd(lam[k][l], lam[k][l], False, budget)
    if twice.sign() < 0:
        twice = -twice
    if not (twice > d[l + 1]):
        return
    var q = _rounded(lam[k][l], d[l + 1], budget)
    var zq = ZZ(q.copy())
    for column in range(len(basis[k])):
        budget.charge(1)
        if basis[l][column].value.sign() != 0:
            basis[k][column] = basis[k][column].sub(zq.mul(basis[l][column], budget), budget)
    comptime if carry:
        for column in range(len(transform[k])):
            budget.charge(1)
            if transform[l][column].value.sign() != 0:
                transform[k][column] = transform[k][column].sub(zq.mul(transform[l][column], budget), budget)
    lam[k][l] = _iadd(lam[k][l], _imul(q, d[l + 1], budget), True, budget)
    for i in range(l):
        budget.charge(1)
        if lam[l][i].sign() != 0:
            lam[k][i] = _iadd(lam[k][i], _imul(q, lam[l][i], budget), True, budget)


def _lll[carry: Bool = True](input: List[List[ZZ]], mut budget: Budget) raises PolyError -> _LLL:
    """Integral LLL with delta 3/4 (Cohen, Algorithm 2.6.7, de Weger).

    Gram-Schmidt data stay exact as the integers d_i (leading Gram
    determinants) and lambda_kj = d_(j+1) mu_kj, so no rational gcds occur;
    every decision equals the rational algorithm's.
    """
    if len(input) == 0 or len(input) > len(input[0]):
        raise PolyError(PolyErrorKind.DOMAIN, "LLL requires a nonempty full-row-rank matrix")
    var rows = len(input)
    var columns = len(input[0])
    for row in input:
        budget.charge(1)
        if len(row) != columns:
            raise PolyError(PolyErrorKind.DOMAIN, "LLL matrix has unequal row lengths")
    var basis = _copy_matrix(input, budget)
    var transform = List[List[ZZ]]()
    comptime if carry:
        _allocation[List[ZZ]](rows, budget)
        _allocation[ZZ](_product(rows, rows), budget)
        transform = List[List[ZZ]](length=rows, fill=List[ZZ](length=rows, fill=ZZ.zero(NoContext())))
        for i in range(rows):
            budget.charge(1)
            transform[i][i] = ZZ.one(NoContext())
    _allocation[Integer](_product(rows, rows), budget)
    var lam = List[List[Integer]](length=rows, fill=List[Integer](length=rows, fill=Integer(0)))
    # d[i] is the Gram determinant of the first i rows; d[0] = 1.
    var d = List[Integer](length=rows + 1, fill=Integer(0))
    d[0] = Integer(1)

    var kmax = 0
    d[1] = _row_dot(basis[0], basis[0], budget)
    if d[1].sign() == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "LLL requires independent nonzero rows")
    var k = 1
    while k < rows:
        budget.charge(1)
        if k > kmax:
            kmax = k
            for j in range(k + 1):
                budget.charge(1)
                var u = _row_dot(basis[k], basis[j], budget)
                for i in range(j):
                    u = _iexact(_iadd(_imul(d[i + 1], u, budget), _imul(lam[k][i], lam[j][i], budget), True, budget), d[i], budget)
                if j < k:
                    lam[k][j] = u^
                else:
                    if u.sign() == 0:
                        raise PolyError(PolyErrorKind.DOMAIN, "LLL requires independent nonzero rows")
                    d[k + 1] = u^
        _reduce[carry](k, k - 1, basis, transform, lam, d, budget)
        var left = _imul(_imul(Integer(4), d[k + 1], budget), d[k - 1], budget)
        var right = _iadd(_imul(_imul(Integer(3), d[k], budget), d[k], budget), _imul(_imul(Integer(4), lam[k][k - 1], budget), lam[k][k - 1], budget), True, budget)
        if left < right:
            _swap_rows(basis, k, k - 1)
            comptime if carry:
                _swap_rows(transform, k, k - 1)
            for j in range(k - 1):
                budget.charge(1)
                var t = lam[k][j].copy()
                lam[k][j] = lam[k - 1][j].copy()
                lam[k - 1][j] = t^
            var l = lam[k][k - 1].copy()
            var b = _iexact(_iadd(_imul(d[k - 1], d[k + 1], budget), _imul(l, l, budget), False, budget), d[k], budget)
            for i in range(k + 1, kmax + 1):
                budget.charge(1)
                var t = lam[i][k].copy()
                lam[i][k] = _iexact(_iadd(_imul(d[k + 1], lam[i][k - 1], budget), _imul(l, t, budget), True, budget), d[k], budget)
                lam[i][k - 1] = _iexact(_iadd(_imul(b, t, budget), _imul(l, lam[i][k], budget), False, budget), d[k + 1], budget)
            d[k] = b^
            k = max(k - 1, 1)
        else:
            var l = k - 2
            while l >= 0:
                _reduce[carry](k, l, basis, transform, lam, d, budget)
                l -= 1
            k += 1
    _allocation[QQ](_product(rows, rows), budget)
    var mu = List[List[QQ]](length=rows, fill=List[QQ](length=rows, fill=QQ.zero(NoContext())))
    var norms = List[QQ](capacity=rows)
    for i in range(rows):
        budget.charge(1)
        for j in range(i):
            mu[i][j] = _ratio(lam[i][j], d[j + 1], budget)
        norms.append(_ratio(d[i + 1], d[i], budget))
    return _LLL(basis^, transform^, mu^, norms^, List[Float64]())


def _to_float(a: Integer) raises PolyError -> Float64:
    var bits = a.magnitude_bit_length()
    try:
        if bits <= 62:
            return Float64(Int(a))
        var shift = bits - 62
        return Float64(Int(a >> shift)) * (Float64(2.0) ** Float64(shift))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _from_float(x: Float64) raises PolyError -> Integer:
    """The integer nearest x (half away from zero)."""
    var r = x + 0.5 if x >= 0 else x - 0.5
    var magnitude = r if r >= 0 else -r
    if magnitude < 4.0e18:
        return Integer(Int(r))
    var shift = 0
    while magnitude >= 4.0e18:
        magnitude *= 0.5
        shift += 1
    try:
        var value = Integer(Int(magnitude)) << shift
        return -value if r < 0 else value^
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _float_vector(row: List[ZZ], mut budget: Budget) raises PolyError -> List[Float64]:
    _allocation[Float64](len(row), budget)
    var result = List[Float64](capacity=len(row))
    for value in row:
        budget.charge(1)
        result.append(_to_float(value.value))
    return result^


def _float_dot(a: List[Float64], b: List[Float64]) -> Float64:
    # Four interleaved partial sums in a fixed order, then the tail.
    var pa = a.unsafe_ptr()
    var pb = b.unsafe_ptr()
    var n = len(a)
    var partial = SIMD[DType.float64, 4](0)
    var c = 0
    while c + 4 <= n:
        partial += pa.unsafe_offset(c).unsafe_load[width=4]() * pb.unsafe_offset(c).unsafe_load[width=4]()
        c += 4
    var total = (partial[0] + partial[1]) + (partial[2] + partial[3])
    while c < n:
        total += pa.unsafe_offset(c).unsafe_load() * pb.unsafe_offset(c).unsafe_load()
        c += 1
    return total


def _gram_row(k: Int, fb: List[List[Float64]], mut gram: List[List[Float64]], mut budget: Budget) raises PolyError:
    _charge_size(_product(len(fb), len(fb[k])), budget)
    for i in range(len(fb)):
        var value = _float_dot(fb[k], fb[i])
        gram[k][i] = value
        gram[i][k] = value


def _float_row(k: Int, gram: List[List[Float64]], mut r: List[List[Float64]], mut mu: List[List[Float64]], mut budget: Budget) raises PolyError:
    # Floating Gram-Schmidt for row k from the floating Gram matrix; the
    # inner sums run in the same order over raw row pointers.
    _charge_size(_product(k + 1, k + 1), budget)
    var rk = r[k].unsafe_ptr()
    var mk = mu[k].unsafe_ptr()
    var gk = gram[k].unsafe_ptr()
    for j in range(k + 1):
        var value = gk.unsafe_offset(j).unsafe_load()
        var mj = mu[j].unsafe_ptr()
        for i in range(j):
            value -= mj.unsafe_offset(i).unsafe_load() * rk.unsafe_offset(i).unsafe_load()
        rk.unsafe_offset(j).unsafe_store(value)
        if j < k:
            var rjj = r[j][j]
            mk.unsafe_offset(j).unsafe_store(value / rjj if rjj != 0 else 0.0)


def _diagonal(r: List[List[Float64]]) -> List[Float64]:
    var result = List[Float64](capacity=len(r))
    for i in range(len(r)):
        result.append(r[i][i])
    return result^


def _exact_norms(basis: List[List[ZZ]], mut budget: Budget) raises PolyError -> List[QQ]:
    """Exact squared Gram-Schmidt norms by the fraction-free recurrence."""
    var rows = len(basis)
    _allocation[Integer](_product(rows, rows), budget)
    var lam = List[List[Integer]](length=rows, fill=List[Integer](length=rows, fill=Integer(0)))
    var d = List[Integer](length=rows + 1, fill=Integer(0))
    d[0] = Integer(1)
    var norms = List[QQ](capacity=rows)
    for k in range(rows):
        budget.charge(1)
        for j in range(k + 1):
            budget.charge(1)
            var u = _row_dot(basis[k], basis[j], budget)
            for i in range(j):
                u = _iexact(_iadd(_imul(d[i + 1], u, budget), _imul(lam[k][i], lam[j][i], budget), True, budget), d[i], budget)
            if j < k:
                lam[k][j] = u^
            else:
                if u.sign() == 0:
                    raise PolyError(PolyErrorKind.DOMAIN, "LLL requires independent nonzero rows")
                d[k + 1] = u^
        norms.append(_ratio(d[k + 1], d[k], budget))
    return norms^


def _gram_entry(a: List[ZZ], b: List[ZZ], small_a: List[Int128], small_b: List[Int128], a_bits: Int, b_bits: Int, mut budget: Budget) raises PolyError -> Integer:
    """The dot product of two rows, in Int128 when no partial sum can pass
    2^126; small_a and small_b are the rows in Int128 when they fit."""
    var width_bits = 64 - Int(count_leading_zeros(UInt64(len(a) + 1)))
    if a_bits + b_bits + width_bits > 126:
        return _row_dot(a, b, budget)
    _charge_size(len(a), budget)
    var total = Int128(0)
    var pa = small_a.unsafe_ptr()
    var pb = small_b.unsafe_ptr()
    for i in range(len(a)):
        total += pa.unsafe_offset(i).unsafe_load() * pb.unsafe_offset(i).unsafe_load()
    var m = UInt128(total) if total >= 0 else UInt128(-total)
    if m < (UInt128(1) << 62):
        return Integer(Int(Int64(total)))
    var value: Integer
    try:
        value = (Integer(Int(UInt64(m >> 64) >> 1)) << 65) + (Integer(Int(UInt64(m >> 64) & 1)) << 64) + (Integer(Int(UInt64(m & UInt128(UInt64.MAX)) >> 1)) << 1) + Integer(Int(UInt64(m) & 1))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return -value if total < 0 else value^


def _word_residue(value: Integer, f: _Field) raises PolyError -> UInt64:
    """value modulo f.p, in Montgomery form."""
    var r: UInt64
    if value.magnitude_bit_length() <= 63:
        try:
            var word = Int(value)
            r = UInt64(word) % f.p if word >= 0 else (f.p - UInt64(-word) % f.p) % f.p
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    else:
        try:
            r = UInt64(Int(div_rem_floor(value, Integer(Int(f.p)))[1]))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return _to(r, f)


def _modular_norms(basis: List[List[ZZ]], mut budget: Budget) raises PolyError -> Optional[List[QQ]]:
    """_exact_norms through the leading principal minors d_k of the Gram
    matrix: d_(k+1)/d_k is the k-th squared norm. Each 0 < d_k is at most
    the product of the Gram diagonal (Hadamard), so primes above 2^61 whose
    product passes it fix every d_k by its residues. No value when a minor
    vanishes modulo one of them."""
    var rows = len(basis)
    var bits = List[Int](capacity=rows)
    var small = List[List[Int128]](capacity=rows)
    for row in basis:
        _charge_size(len(row), budget)
        var largest = 0
        for v in row:
            largest = max(largest, v.value.magnitude_bit_length())
        bits.append(largest)
        var values = List[Int128](capacity=len(row) if largest < 126 else 0)
        if largest < 126:
            for v in row:
                values.append(_int128(v.value))
        small.append(values^)
    _allocation[Integer](_product(rows, rows), budget)
    var gram = List[List[Integer]](length=rows, fill=List[Integer]())
    var hadamard = 1
    for i in range(rows):
        budget.charge(1)
        var row = List[Integer](capacity=rows)
        for j in range(rows):
            row.append(_gram_entry(basis[i], basis[j], small[i], small[j], bits[i], bits[j], budget) if j >= i else gram[j][i].copy())
        hadamard = _sum(hadamard, row[i].magnitude_bit_length())
        gram[i] = row^
    var count = hadamard // 61 + 1
    _allocation[UInt64](_product(_sum(rows, 1), count), budget)
    var primes = List[UInt64](capacity=count)
    var fields = List[_Field](capacity=count)
    for t in range(count):
        primes.append(prime(t, budget))
        fields.append(_field(primes[t]))
    # residues[k][t] = d_(k+1) modulo primes[t], in Montgomery form.
    var residues = List[List[UInt64]](length=rows, fill=List[UInt64](length=count, fill=0))
    _allocation[UInt64](_product(rows, rows), budget)
    for t in range(count):
        ref f = fields[t]
        var g = List[List[UInt64]](capacity=rows)
        for i in range(rows):
            _charge_size(rows, budget)
            var row = List[UInt64](capacity=rows)
            for j in range(rows):
                row.append(_word_residue(gram[i][j], f) if j >= i else 0)
            g.append(row^)
        # Symmetric elimination on the upper triangle: the pivots' running
        # products are the leading minors.
        var minor = f.one
        for c in range(rows):
            var pivot = g[c][c]
            if pivot == 0:
                return None
            minor = _mulm(minor, pivot, f)
            residues[c][t] = minor
            var inverse = _invm(pivot, f, budget)
            for i in range(c + 1, rows):
                var factor = _mulm(g[c][i], inverse, f)
                if factor == 0:
                    continue
                _charge_size(rows - i, budget)
                var pc = g[c].unsafe_ptr()
                var pi = g[i].unsafe_ptr()
                for j in range(i, rows):
                    pi.unsafe_offset(j).unsafe_store(_subm(pi.unsafe_offset(j).unsafe_load(), _mulm(factor, pc.unsafe_offset(j).unsafe_load(), f), f))
    # Garner: inverses[t] = (p_0 ... p_(t-1))^-1 modulo p_t.
    var inverses = List[UInt64](capacity=count)
    for t in range(count):
        ref f = fields[t]
        var product = f.one
        for j in range(t):
            product = _mulm(product, _to(primes[j], f), f)
        inverses.append(_invm(product, f, budget) if t else f.one)
    var minors = List[Integer](capacity=rows + 1)
    minors.append(Integer(1))
    for k in range(rows):
        _charge_size(_product(count, count), budget)
        var digits = List[UInt64](capacity=count)
        for t in range(count):
            ref f = fields[t]
            # The value of the digits so far, modulo p_t, by Horner.
            var partial = UInt64(0)
            var j = t - 1
            while j >= 0:
                partial = _addm(_mulm(partial, _to(primes[j], f), f), _to(digits[j], f), f)
                j -= 1
            digits.append(_from(_mulm(_subm(residues[k][t], partial, f), inverses[t], f), f))
        var value = Integer(Int(digits[count - 1]))
        var j = count - 2
        while j >= 0:
            value = _iadd(_imul(value, Integer(Int(primes[j])), budget), Integer(Int(digits[j])), False, budget)
            j -= 1
        if value.sign() == 0:
            raise PolyError(PolyErrorKind.DOMAIN, "LLL requires independent nonzero rows")
        minors.append(value^)
    var norms = List[QQ](capacity=rows)
    for k in range(rows):
        norms.append(_ratio(minors[k + 1], minors[k], budget))
    return norms^


comptime _DELTA = 0.9
"""The Lovasz parameter of the floating-steered reduction."""


comptime _SMALL_BITS = 100
"""Entries below 2^100 run the knapsack reduction on Int128; a reduction
step whose result could pass 2^126 abandons it for Integer arithmetic."""


def _bits128(v: Int128) -> Int:
    var m = UInt128(v) if v >= 0 else UInt128(-v)
    var high = UInt64(m >> 64)
    if high != 0:
        return 128 - Int(count_leading_zeros(high))
    var low = UInt64(m)
    return 64 - Int(count_leading_zeros(low)) if low != 0 else 0


def _float128(v: Int128) -> Float64:
    # As _to_float: values beyond 62 bits keep their top 62 by a floor shift.
    var bits = _bits128(v)
    if bits <= 62:
        return Float64(Int64(v))
    var shift = bits - 62
    return Float64(Int64(v >> Int128(shift))) * (Float64(2.0) ** Float64(shift))


def _int128(a: Integer) raises PolyError -> Int128:
    """a, of magnitude below 2^126, by 32-bit chunks of its magnitude."""
    if a.magnitude_bit_length() <= 62:
        try:
            return Int128(Int(a))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var m = a.copy()
    if a.sign() < 0:
        m = -m
    var result = Int128(0)
    try:
        for k in range(4):
            var chunk = Int((m >> (32 * k)) & Integer(0xFFFFFFFF))
            result |= Int128(chunk) << Int128(32 * k)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return -result if a.sign() < 0 else result


def _row_bits(row: List[Int128]) -> Int:
    # The bit length of the largest magnitude: OR of magnitudes suffices.
    var acc = UInt128(0)
    for v in row:
        acc |= UInt128(v) if v >= 0 else UInt128(-v)
    var high = UInt64(acc >> 64)
    if high != 0:
        return 128 - Int(count_leading_zeros(high))
    var low = UInt64(acc)
    return 64 - Int(count_leading_zeros(low)) if low != 0 else 0


def _lll_knapsack_small(input: List[List[ZZ]], mut budget: Budget, exact_norms: Bool = True) raises PolyError -> Optional[_LLL]:
    """_lll_knapsack on Int128 entries: the same floating decisions and
    exact row operations, so the same reduced basis, or no value when an
    entry or a step could overflow."""
    var rows = len(input)
    var basis = List[List[Int128]](capacity=rows)
    for row in input:
        budget.charge(1)
        var values = List[Int128](capacity=len(row))
        for v in row:
            if v.value.magnitude_bit_length() >= _SMALL_BITS:
                return None
            values.append(_int128(v.value))
        basis.append(values^)
    var widths = List[Int](capacity=rows)
    for row in basis:
        widths.append(_row_bits(row))
    _allocation[Float64](_product(2 * rows, rows), budget)
    var r = List[List[Float64]](length=rows, fill=List[Float64](length=rows, fill=0.0))
    var mu = List[List[Float64]](length=rows, fill=List[Float64](length=rows, fill=0.0))
    var fb = List[List[Float64]](capacity=rows)
    for row in basis:
        var values = List[Float64](capacity=len(row))
        for v in row:
            values.append(_float128(v))
        fb.append(values^)
    var gram = List[List[Float64]](length=rows, fill=List[Float64](length=rows, fill=0.0))
    for i in range(rows):
        _gram_row(i, fb, gram, budget)
    _float_row(0, gram, r, mu, budget)
    var k = 1
    # Whether row k's floating Gram-Schmidt data is current: a swap leaves
    # it so for the row that becomes k, and computing it again would
    # repeat the same arithmetic on the same Gram row.
    var current = False
    var iterations = 0
    var cap = 10000 + 100 * rows * rows
    while k < rows:
        budget.charge(1)
        iterations += 1
        if iterations > cap:
            return None
        var rounds = 0
        while True:
            budget.charge(1)
            if not current:
                _float_row(k, gram, r, mu, budget)
            current = False
            var reduced = False
            var j = k - 1
            while j >= 0:
                if mu[k][j] > 0.51 or mu[k][j] < -0.51:
                    # As _from_float, _int128 and _to_float, without an
                    # Integer while the rounded value fits a machine word.
                    var xs: Int128
                    var fx: Float64
                    var rounded = mu[k][j] + 0.5 if mu[k][j] >= 0 else mu[k][j] - 0.5
                    if (rounded if rounded >= 0 else -rounded) < 4.0e18:
                        var word = Int(rounded)
                        xs = Int128(word)
                        fx = Float64(word)
                    else:
                        var x = _from_float(mu[k][j])
                        if x.magnitude_bit_length() > 125:
                            return None
                        xs = _int128(x)
                        fx = _to_float(x)
                    # widths holds upper bounds on the row bit lengths; the
                    # exact ones decide whether the step could overflow.
                    var xbits = _bits128(xs)
                    if xbits + widths[j] > 125 or widths[k] > 125:
                        widths[j] = _row_bits(basis[j])
                        widths[k] = _row_bits(basis[k])
                        if xbits + widths[j] > 125 or widths[k] > 125:
                            return None
                    _charge_size(len(basis[k]), budget)
                    var target = basis[k].unsafe_ptr()
                    var source = basis[j].unsafe_ptr()
                    for column in range(len(basis[k])):
                        var b = source.unsafe_offset(column).unsafe_load()
                        if b != 0:
                            target.unsafe_offset(column).unsafe_store(target.unsafe_offset(column).unsafe_load() - xs * b)
                    widths[k] = max(widths[k], xbits + widths[j]) + 1
                    var mk = mu[k].unsafe_ptr()
                    var mj = mu[j].unsafe_ptr()
                    for i in range(j):
                        mk.unsafe_offset(i).unsafe_store(mk.unsafe_offset(i).unsafe_load() - fx * mj.unsafe_offset(i).unsafe_load())
                    mu[k][j] -= fx
                    reduced = True
                j -= 1
            if reduced:
                for column in range(len(basis[k])):
                    fb[k][column] = _float128(basis[k][column])
                _gram_row(k, fb, gram, budget)
            rounds += 1
            if not reduced or rounds > 50:
                break
        var lhs = _DELTA * r[k - 1][k - 1]
        var rhs = r[k][k] + mu[k][k - 1] * mu[k][k - 1] * r[k - 1][k - 1]
        if lhs > rhs:
            _swap_lists(basis, k, k - 1)
            var held_width = widths[k]
            widths[k] = widths[k - 1]
            widths[k - 1] = held_width
            _swap_lists(fb, k, k - 1)
            _swap_lists(gram, k, k - 1)
            for i in range(rows):
                var gi = gram[i].unsafe_ptr()
                var t = gi.unsafe_offset(k).unsafe_load()
                gi.unsafe_offset(k).unsafe_store(gi.unsafe_offset(k - 1).unsafe_load())
                gi.unsafe_offset(k - 1).unsafe_store(t)
            _float_row(k - 1, gram, r, mu, budget)
            current = k > 1
            k = max(k - 1, 1)
        else:
            k += 1
    var exact = List[List[ZZ]](capacity=rows)
    for row in basis:
        budget.charge(1)
        var values = List[ZZ](capacity=len(row))
        for v in row:
            var m = UInt128(v) if v >= 0 else UInt128(-v)
            var value: Integer
            if m < (UInt128(1) << 62):
                value = Integer(Int(UInt64(m)))
                values.append(ZZ(value^) if v >= 0 else ZZ(value^).neg(budget))
                continue
            try:
                value = (Integer(Int(UInt64(m >> 64) >> 1)) << 65) + (Integer(Int(UInt64(m >> 64) & 1)) << 64) + (Integer(Int(UInt64(m & UInt128(UInt64.MAX)) >> 1)) << 1) + Integer(Int(UInt64(m) & 1))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            values.append(ZZ(value^) if v >= 0 else ZZ(value^).neg(budget))
        exact.append(values^)
    var norms = _exact_norms(exact, budget) if exact_norms else List[QQ]()
    return _LLL(exact^, List[List[ZZ]](), List[List[QQ]](), norms^, _diagonal(r))


def _lll_knapsack(input: List[List[ZZ]], mut budget: Budget, exact_norms: Bool = True) raises PolyError -> _LLL:
    """LLL steered by floating-point Gram-Schmidt data, with exact norms.

    Row operations are exact integer ones, so the result is always a basis
    of the input lattice; floating point only chooses them. The returned
    squared norms are recomputed exactly, so every decision a caller bases
    on them is rigorous. A deterministic iteration cap falls back to the
    exact integral algorithm when floating point stalls.
    """
    if len(input) == 0 or len(input) > len(input[0]):
        raise PolyError(PolyErrorKind.DOMAIN, "LLL requires a nonempty full-row-rank matrix")
    var small = _lll_knapsack_small(input, budget, exact_norms)
    if small:
        return small.take()
    return _lll_knapsack_integer(input, budget, exact_norms)


def _lll_knapsack_integer(input: List[List[ZZ]], mut budget: Budget, exact_norms: Bool = True) raises PolyError -> _LLL:
    """_lll_knapsack on Integer entries, for any size."""
    var rows = len(input)
    var basis = _copy_matrix(input, budget)
    _allocation[Float64](_product(2 * rows, rows), budget)
    var r = List[List[Float64]](length=rows, fill=List[Float64](length=rows, fill=0.0))
    var mu = List[List[Float64]](length=rows, fill=List[Float64](length=rows, fill=0.0))
    var fb = List[List[Float64]](capacity=rows)
    for row in basis:
        fb.append(_float_vector(row, budget))
    var gram = List[List[Float64]](length=rows, fill=List[Float64](length=rows, fill=0.0))
    for i in range(rows):
        _gram_row(i, fb, gram, budget)
    _float_row(0, gram, r, mu, budget)
    var k = 1
    # Whether row k's floating Gram-Schmidt data is current: a swap leaves
    # it so for the row that becomes k, and computing it again would
    # repeat the same arithmetic on the same Gram row.
    var current = False
    var iterations = 0
    var cap = 10000 + 100 * rows * rows
    while k < rows:
        budget.charge(1)
        iterations += 1
        if iterations > cap:
            return _lll[False](input, budget)
        var rounds = 0
        while True:
            budget.charge(1)
            if not current:
                _float_row(k, gram, r, mu, budget)
            current = False
            var reduced = False
            var j = k - 1
            while j >= 0:
                if mu[k][j] > 0.51 or mu[k][j] < -0.51:
                    var x = _from_float(mu[k][j])
                    var zx = ZZ(x.copy())
                    for column in range(len(basis[k])):
                        budget.charge(1)
                        if basis[j][column].value.sign() != 0:
                            basis[k][column] = basis[k][column].sub(zx.mul(basis[j][column], budget), budget)
                    var fx = _to_float(x)
                    for i in range(j):
                        mu[k][i] -= fx * mu[j][i]
                    mu[k][j] -= fx
                    reduced = True
                j -= 1
            if reduced:
                fb[k] = _float_vector(basis[k], budget)
                _gram_row(k, fb, gram, budget)
            rounds += 1
            if not reduced or rounds > 50:
                break
        var lhs = _DELTA * r[k - 1][k - 1]
        var rhs = r[k][k] + mu[k][k - 1] * mu[k][k - 1] * r[k - 1][k - 1]
        if lhs > rhs:
            _swap_rows(basis, k, k - 1)
            _swap_lists(fb, k, k - 1)
            _swap_lists(gram, k, k - 1)
            for i in range(rows):
                var t = gram[i][k]
                gram[i][k] = gram[i][k - 1]
                gram[i][k - 1] = t
            _float_row(k - 1, gram, r, mu, budget)
            current = k > 1
            k = max(k - 1, 1)
        else:
            k += 1
    var norms = _exact_norms(basis, budget) if exact_norms else List[QQ]()
    return _LLL(basis^, List[List[ZZ]](), List[List[QQ]](), norms^, _diagonal(r))
