"""Exact endpoint signs for the private final-CAD verifier.

Eight homogeneous Horner evaluations share a SIMD instruction stream. A
verifier-computed absolute-value bound proves that every intermediate fits
signed 64-bit arithmetic before the vector kernel runs. Larger batches use
arbitrary-precision integers. No floating-point decision enters acceptance.
"""
from apn_mojo import Integer, lcm
from .coefficients import QQ, _work
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum


struct _SignStats(Copyable, Movable):
    var word_points: Int
    var big_points: Int
    var vector_batches: Int

    def __init__(out self):
        self.word_points = 0
        self.big_points = 0
        self.vector_batches = 0


def _integer_coefficients(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[Integer]:
    _allocation[Integer](len(p.coefficients), budget)
    var common = Integer(1)
    var result = List[Integer](capacity=len(p.coefficients))
    for c in p.coefficients:
        _work(common.magnitude_bit_length() // 8 + 1, c.size(), True, budget)
        try:
            common = lcm(common, c.value.denominator())
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    for c in p.coefficients:
        _work(common.magnitude_bit_length() // 8 + 1, c.size(), True, budget)
        try:
            result.append(c.value.numerator() * (common // c.value.denominator()))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return result^


def _native(x: Integer) raises PolyError -> Optional[Int64]:
    # Deliberately exclude Int64.MIN so abs(x) is representable too.
    if x.magnitude_bit_length() > 63:
        return None
    try:
        return Int64(Int(x))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


@always_inline
def _abs_word(x: Int64) -> UInt64:
    return UInt64(-x if x < 0 else x)


def _fits(c: List[Int64], a: UInt64, b: UInt64, mut budget: Budget) raises PolyError -> Bool:
    """Bound power, each product, and each sum in homogeneous Horner."""
    if len(c) == 0:
        return True
    var limit = UInt64(Int64.MAX)
    var value = _abs_word(c[len(c) - 1])
    var power = UInt64(1)
    for i in range(len(c) - 2, -1, -1):
        budget.charge(1)
        if b > limit // power:
            return False
        power *= b
        if a and value > limit // a:
            return False
        value *= a
        var coefficient = _abs_word(c[i])
        if coefficient > (limit - value) // power:
            return False
        value += coefficient * power
    return True


def _big_sign(c: List[Integer], x: QQ, mut budget: Budget) raises PolyError -> Int:
    if len(c) == 0:
        return 0
    var value = c[len(c) - 1].copy()
    var power = Integer(1)
    var a = x.value.numerator()
    var b = x.value.denominator()
    for i in range(len(c) - 2, -1, -1):
        _work(value.magnitude_bit_length() // 8 + 1, x.size(), True, budget)
        _work(power.magnitude_bit_length() // 8 + 1, x.size(), True, budget)
        try:
            power = power * b
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        _work(c[i].magnitude_bit_length() // 8 + 1, power.magnitude_bit_length() // 8 + 1, True, budget)
        try:
            value = value * a + c[i] * power
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return value.sign()


def _batches[lanes: Int](c: List[Integer], small: List[Int64], points: List[QQ], mut stats: _SignStats, mut budget: Budget) raises PolyError -> List[Int]:
    _allocation[Int](len(points), budget)
    var result = List[Int](capacity=len(points))
    for start in range(0, len(points), lanes):
        var count = min(lanes, len(points) - start)
        var numerators = SIMD[DType.int64, lanes](0)
        var denominators = SIMD[DType.int64, lanes](1)
        var largest_a = UInt64(0)
        var largest_b = UInt64(1)
        var fits = len(small) == len(c)
        for j in range(count):
            budget.charge(1)
            var a = _native(points[start + j].value.numerator())
            var b = _native(points[start + j].value.denominator())
            if not a or not b:
                fits = False
            else:
                numerators[j] = a.value()
                denominators[j] = b.value()
                largest_a = max(largest_a, _abs_word(a.value()))
                largest_b = max(largest_b, UInt64(b.value()))
        if fits and _fits(small, largest_a, largest_b, budget):
            var value = SIMD[DType.int64, lanes](small[len(small) - 1] if len(small) else 0)
            var power = SIMD[DType.int64, lanes](1)
            for i in range(len(small) - 2, -1, -1):
                budget.charge(UInt64(count))
                power *= denominators
                value = value * numerators + SIMD[DType.int64, lanes](small[i]) * power
            for j in range(count):
                result.append(1 if value[j] > 0 else (-1 if value[j] < 0 else 0))
            stats.word_points = _sum(stats.word_points, count)
        else:
            for j in range(count):
                result.append(_big_sign(c, points[start + j], budget))
            stats.big_points = _sum(stats.big_points, count)
    return result^


@fieldwise_init
struct _PreparedSigns(Copyable, Movable):
    """Verifier-derived denominator clearing and native coefficients."""
    var integers: List[Integer]
    var native: List[Int64]


@fieldwise_init
struct _EndpointBatch(Copyable, Movable):
    var numerators: SIMD[DType.int64, 8]
    var denominators: SIMD[DType.int64, 8]
    var largest_a: UInt64
    var largest_b: UInt64
    var native: Bool
    var count: Int


def _prepare_signs(p: DensePoly[QQ], mut budget: Budget) raises PolyError -> _PreparedSigns:
    p.check_invariants(budget)
    var c = _integer_coefficients(p, budget)
    _allocation[Int64](len(c), budget)
    var small = List[Int64](capacity=len(c))
    for coefficient in c:
        budget.charge(1)
        var value = _native(coefficient)
        if not value:
            small.clear()
            break
        small.append(value.value())
    return _PreparedSigns(c^, small^)


def _pack_endpoints(points: List[QQ], mut budget: Budget) raises PolyError -> List[_EndpointBatch]:
    _allocation[_EndpointBatch]((len(points) + 7) // 8, budget)
    var result = List[_EndpointBatch]()
    for start in range(0, len(points), 8):
        var count = min(8, len(points) - start)
        var a = SIMD[DType.int64, 8](0)
        var b = SIMD[DType.int64, 8](1)
        var largest_a = UInt64(0)
        var largest_b = UInt64(1)
        var native = True
        for j in range(count):
            budget.charge(1)
            var numerator = _native(points[start + j].value.numerator())
            var denominator = _native(points[start + j].value.denominator())
            if not numerator or not denominator:
                native = False
            else:
                a[j] = numerator.value()
                b[j] = denominator.value()
                largest_a = max(largest_a, _abs_word(numerator.value()))
                largest_b = max(largest_b, UInt64(denominator.value()))
        result.append(_EndpointBatch(a, b, largest_a, largest_b, native, count))
    return result^


def _prepared_signs_at(prepared: _PreparedSigns, points: List[QQ], batches: List[_EndpointBatch], mut stats: _SignStats, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> List[Int]:
    """Reuse verifier-owned coefficient and endpoint packing across a chain."""
    ref c = prepared.integers
    ref small = prepared.native
    if len(c) <= 1:
        _allocation[Int](len(points), budget)
        budget.charge(UInt64(len(points)))
        if len(small) == len(c):
            stats.word_points = _sum(stats.word_points, len(points))
        else:
            stats.big_points = _sum(stats.big_points, len(points))
        return List[Int](length=len(points), fill=c[0].sign() if len(c) else 0)
    if not vectorized or len(points) < 8:
        return _batches[1](c, small, points, stats, budget)
    if len(batches) != (len(points) + 7) // 8:
        raise PolyError(PolyErrorKind.DOMAIN, "incomplete endpoint packing")
    _allocation[Int](len(points), budget)
    var result = List[Int](capacity=len(points))
    var start = 0
    for batch in batches:
        if batch.native and len(small) == len(c) and _fits(small, batch.largest_a, batch.largest_b, budget):
            var value = SIMD[DType.int64, 8](small[len(small) - 1] if len(small) else 0)
            var power = SIMD[DType.int64, 8](1)
            for i in range(len(small) - 2, -1, -1):
                budget.charge(UInt64(batch.count))
                power *= batch.denominators
                value = value * batch.numerators + SIMD[DType.int64, 8](small[i]) * power
            for j in range(batch.count):
                result.append(1 if value[j] > 0 else (-1 if value[j] < 0 else 0))
            stats.word_points = _sum(stats.word_points, batch.count)
            stats.vector_batches = _sum(stats.vector_batches, 1)
        else:
            for j in range(batch.count):
                result.append(_big_sign(c, points[start + j], budget))
            stats.big_points = _sum(stats.big_points, batch.count)
        start += batch.count
    return result^


def _signs_at(p: DensePoly[QQ], points: List[QQ], mut stats: _SignStats, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> List[Int]:
    """Exact signs; scalar remains the word/arbitrary-integer reference."""
    var prepared = _prepare_signs(p, budget)
    var batches = List[_EndpointBatch]()
    if vectorized and len(points) >= 8:
        batches = _pack_endpoints(points, budget)
    return _prepared_signs_at(prepared, points, batches, stats, budget, vectorized=vectorized)
