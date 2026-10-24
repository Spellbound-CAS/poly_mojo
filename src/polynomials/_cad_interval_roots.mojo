"""Bounded interval-first producer for relative field roots.

Adapted from this repository's exp-roots part-2 producer, commit
465ec59ba0cbd79822d7fafc1727bb0b66533a46. No certificate or search replay is
retained. First try checked Int64 Descartes transforms, then wider exact dyadic
enclosures. An undecided attempt returns None for the exact Sturm fallback;
the final stack verifier independently checks all roots and coverage.
"""
from apn_mojo import Integer, Rational
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product, _charge_size
from ._algebraic import _refine
from ._algebraic_filter import _dense_interval, _constant
from ._dyadic import _bounds, _left
from ._root_ops import _compare, _midpoint
from ._cad_sample_field import _SampleField, _kevaluate, _kexact, _ksign, _kevaluate_sign
from ._cad_final import _FinalRoot
from ._cad_final_signs import _native
from ._cad_native import _shift_native


struct _IntervalRootStats(Copyable, Movable):
    var native_nodes: Int
    var wide_nodes: Int
    var attempts: Int
    var exact_fallbacks: Int

    def __init__(out self):
        self.native_nodes = 0
        self.wide_nodes = 0
        self.attempts = 0
        self.exact_fallbacks = 0


def _iadd(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    _charge_size(_sum(a.magnitude_bit_length(), b.magnitude_bit_length()) // 64 + 1, budget)
    try:
        return a + b
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _ishift(value: Integer, shift: Int, mut budget: Budget) raises PolyError -> Integer:
    _charge_size(_sum(value.magnitude_bit_length(), shift if value else 0) // 64 + 1, budget)
    return _left(value, shift)


def _isign(low: Integer, high: Integer) -> Int:
    if low.sign() > 0:
        return 1
    if high.sign() < 0:
        return -1
    return 0 if not low and not high else 2


def _shift_one(mut low: List[Integer], mut high: List[Integer], mut budget: Budget) raises PolyError:
    var n = len(low) - 1
    for i in range(n):
        for j in range(n - 1, i - 1, -1):
            low[j] = _iadd(low[j], low[j + 1], budget)
            high[j] = _iadd(high[j], high[j + 1], budget)


def _descartes[paired: Bool = False](low: List[Integer], high: List[Integer], native: Bool, mut stats: _IntervalRootStats, mut budget: Budget) raises PolyError -> Optional[Int]:
    """Variations after reversal and shift by one; -1 means undecided.

    Paired native SIMD is an explicit experiment. Complete-CAD measurements
    did not establish an improvement, so the default keeps the scalar path.
    """
    var n = len(low) - 1
    var previous = 0
    var changes = 0
    if native:
        var largest = 0
        for i in range(len(low)):
            budget.charge(1)
            largest = max(largest, max(low[i].magnitude_bit_length(), high[i].magnitude_bit_length()))
        # Every Taylor-shift intermediate is a sum of binomial multiples,
        # bounded in absolute value by 2^n times the largest input magnitude.
        if _sum(largest, n) > 63:
            return None
        _allocation[Int64](_product(2, len(low)), budget)
        comptime if paired:
            # Pair the independent lower/upper transforms in native lanes. The
            # guard above bounds every intermediate; retain the same abstract
            # allocation/work charge as the two separate Int64 vectors.
            var coefficients = List[SIMD[DType.int64, 2]](capacity=len(low))
            for i in range(n, -1, -1):
                var lo = _native(low[i])
                var hi = _native(high[i])
                if not lo or not hi:
                    return None
                coefficients.append(SIMD[DType.int64, 2](lo.value(), hi.value()))
            budget.charge(UInt64(_product(n, n + 1)))
            _shift_native[True](coefficients)
            stats.native_nodes = _sum(stats.native_nodes, 1)
            for i in range(n + 1):
                var a = coefficients[i][0]
                var b = coefficients[i][1]
                var sign = 1 if a > 0 else (-1 if b < 0 else (0 if a == 0 and b == 0 else 2))
                if sign == 2:
                    return -1
                if sign:
                    changes += Int(previous != 0 and previous != sign)
                    previous = sign
        else:
            var a = List[Int64]()
            var b = List[Int64]()
            for i in range(n, -1, -1):
                var lo = _native(low[i])
                var hi = _native(high[i])
                if not lo or not hi:
                    return None
                a.append(lo.value())
                b.append(hi.value())
            for i in range(n):
                for j in range(n - 1, i - 1, -1):
                    budget.charge(2)
                    a[j] += a[j + 1]
                    b[j] += b[j + 1]
            stats.native_nodes = _sum(stats.native_nodes, 1)
            for i in range(n + 1):
                var sign = 1 if a[i] > 0 else (-1 if b[i] < 0 else (0 if a[i] == 0 and b[i] == 0 else 2))
                if sign == 2:
                    return -1
                if sign:
                    changes += Int(previous != 0 and previous != sign)
                    previous = sign
    else:
        _allocation[Integer](_product(2, len(low)), budget)
        var a = List[Integer]()
        var b = List[Integer]()
        for i in range(n, -1, -1):
            a.append(low[i].copy())
            b.append(high[i].copy())
        _shift_one(a, b, budget)
        stats.wide_nodes = _sum(stats.wide_nodes, 1)
        for i in range(n + 1):
            var sign = _isign(a[i], b[i])
            if sign == 2:
                return -1
            if sign:
                changes += Int(previous != 0 and previous != sign)
                previous = sign
    return changes


def _end(index: Integer, depth: Int, bound: Int, negative: Bool, mut budget: Budget) raises PolyError -> QQ:
    var exponent = bound - depth
    var numerator = _ishift(index, max(exponent, 0), budget)
    var denominator = _ishift(Integer(1), max(-exponent, 0), budget)
    try:
        var value = Rational(numerator^, denominator^)
        return QQ(-value if negative else value)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


@fieldwise_init
struct _IntervalNode(Copyable, Movable):
    var low: List[Integer]
    var high: List[Integer]
    var index: Integer
    var depth: Int


def _isolate_side(low: List[Integer], high: List[Integer], bound: Int, negative: Bool, native: Bool, precision: Int, mut brackets: List[Tuple[QQ, QQ]], mut split: QQ, mut stats: _IntervalRootStats, mut budget: Budget) raises PolyError -> Int:
    # Status 0 = complete, 1 = retry, 2 = inspect an undecided split exactly.
    var n = len(low) - 1
    _allocation[Integer](_product(2, len(low)), budget)
    var a = List[Integer]()
    var b = List[Integer]()
    for i in range(n + 1):
        var lo = _ishift(low[i], _product(bound, i), budget)
        var hi = _ishift(high[i], _product(bound, i), budget)
        if negative and i % 2:
            a.append(-hi)
            b.append(-lo)
        else:
            a.append(lo^)
            b.append(hi^)
    _allocation[_IntervalNode](1, budget)
    var pending: List[_IntervalNode] = [_IntervalNode(a^, b^, Integer(0), 0)]
    var nodes = 0
    while len(pending):
        budget.charge(1)
        nodes += 1
        if nodes > 4096:
            return 1
        var node = pending.pop()
        var variation = _descartes(node.low, node.high, native, stats, budget)
        if not variation:
            return 1
        var count = variation.value()
        if count == 0:
            continue
        if count == 1:
            var lo = _end(node.index, node.depth, bound, negative, budget)
            var hi = _end(_iadd(node.index, Integer(1), budget), node.depth, bound, negative, budget)
            _allocation[Tuple[QQ, QQ]](1, budget)
            brackets.append((hi^, lo^) if negative else (lo^, hi^))
            continue
        if node.depth >= precision or (count < 0 and node.depth >= 48):
            return 1
        _allocation[Integer](_product(4, len(low)), budget)
        var left_low = List[Integer]()
        var left_high = List[Integer]()
        for i in range(n + 1):
            left_low.append(_ishift(node.low[i], n - i, budget))
            left_high.append(_ishift(node.high[i], n - i, budget))
        var right_low = left_low.copy()
        var right_high = left_high.copy()
        _shift_one(right_low, right_high, budget)
        var index = _ishift(node.index, 1, budget)
        var right_index = _iadd(index, Integer(1), budget)
        if _isign(right_low[0], right_high[0]) == 0 or _isign(right_low[0], right_high[0]) == 2:
            split = _end(right_index, node.depth + 1, bound, negative, budget)
            return 2
        _allocation[_IntervalNode](2, budget)
        pending.append(_IntervalNode(right_low^, right_high^, right_index^, node.depth + 1))
        pending.append(_IntervalNode(left_low^, left_high^, index^, node.depth + 1))
    return 0


def _isolate_intervals(f: List[DensePoly[QQ]], field: _SampleField, precision: Int, native: Bool, mut brackets: List[Tuple[QQ, QQ]], mut split: QQ, mut stats: _IntervalRootStats, mut budget: Budget) raises PolyError -> Int:
    var n = len(f) - 1
    _allocation[Tuple[Integer, Integer, Int]](len(f), budget)
    var mantissas = List[Tuple[Integer, Integer, Int]]()
    var exponent = Int.MAX
    for coefficient in f:
        var low = QQ.zero(NoContext())
        var high = QQ.zero(NoContext())
        if field.embedding:
            var box = _dense_interval(coefficient, field.embedding.value(), budget)
            low = box[0].copy()
            high = box[1].copy()
        elif not coefficient.is_zero():
            low = coefficient.coefficients[0].copy()
            high = low.copy()
        var dyadic = _bounds(low, high, precision, budget)
        if dyadic.low or dyadic.high:
            exponent = min(exponent, dyadic.exponent)
        mantissas.append((dyadic.low.copy(), dyadic.high.copy(), dyadic.exponent))
    _allocation[Integer](_product(2, len(f)), budget)
    var low = List[Integer]()
    var high = List[Integer]()
    for m in mantissas:
        low.append(_ishift(m[0], max(0, m[2] - exponent), budget))
        high.append(_ishift(m[1], max(0, m[2] - exponent), budget))
    var leading = _isign(low[n], high[n])
    var constant = _isign(low[0], high[0])
    if leading == 0 or leading == 2 or constant == 0 or constant == 2:
        return 1
    var top = 0
    for i in range(n):
        top = max(top, max(low[i].magnitude_bit_length(), high[i].magnitude_bit_length()))
    var least = min(low[n].magnitude_bit_length(), high[n].magnitude_bit_length())
    var bound = max(top - least + 2, 1)
    var status = _isolate_side(low, high, bound, True, native, precision, brackets, split, stats, budget)
    if status:
        return status
    return _isolate_side(low, high, bound, False, native, precision, brackets, split, stats, budget)


def _deflate_at(f: List[DensePoly[QQ]], root: QQ, field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    var divisor: List[DensePoly[QQ]] = [_constant(root.neg(budget), budget), _constant(QQ.one(NoContext()), budget)]
    return _kexact(f, divisor, field, budget)


def _one_basis(f: List[DensePoly[QQ]], owner: Int, mut field: _SampleField, mut stats: _IntervalRootStats, mut budget: Budget) raises PolyError -> Optional[List[_FinalRoot]]:
    # A monic linear root is already an element of the sample field. Enclose
    # that element directly instead of transforming the polynomial or doing
    # a fixed number of generator refinements for every linear basis member.
    if len(f) == 2:
        var low: QQ
        var high: QQ
        if f[0].degree() <= 0:
            low = f[0].coefficients[0].neg(budget) if not f[0].is_zero() else QQ.zero(NoContext())
            high = low.copy()
        else:
            if not field.embedding:
                raise PolyError(PolyErrorKind.DOMAIN, "algebraic linear root without an embedding")
            var enclosure = _dense_interval(f[0], field.embedding.value(), budget)
            low = enclosure[1].neg(budget)
            high = enclosure[0].neg(budget)
        _allocation[_FinalRoot](1, budget)
        var result: List[_FinalRoot] = [_FinalRoot(owner, low^, high^)]
        return result^
    _allocation[DensePoly[QQ]](len(f), budget)
    var remaining = f.copy()
    var exact = List[_FinalRoot]()
    var zero = QQ.zero(NoContext())
    if remaining[0].is_zero():
        _allocation[_FinalRoot](1, budget)
        exact.append(_FinalRoot(owner, zero.copy(), zero.copy()))
        remaining = _deflate_at(remaining, zero, field, budget)
    var precisions: List[Int] = [24, 64, 128, 256]
    var refinements: List[Int] = [8, 16, 40, 96]
    for attempt in range(len(precisions)):
        stats.attempts = _sum(stats.attempts, 1)
        if field.embedding:
            for _ in range(refinements[attempt]):
                _refine(field.embedding.value(), budget)
        while True:
            if len(remaining) <= 1:
                return exact^
            var brackets = List[Tuple[QQ, QQ]]()
            var split = zero.copy()
            var status = _isolate_intervals(remaining, field, precisions[attempt], attempt == 0, brackets, split, stats, budget)
            if status == 2 and _kevaluate_sign(remaining, split, field, budget) == 0:
                _allocation[_FinalRoot](1, budget)
                exact.append(_FinalRoot(owner, split.copy(), split.copy()))
                remaining = _deflate_at(remaining, split, field, budget)
                continue
            if status == 0:
                var removed_count = len(exact)
                for bracket in brackets:
                    var low = bracket[0].copy()
                    var high = bracket[1].copy()
                    # The Descartes brackets isolate the deflated polynomial.
                    # Exclude every removed root, including endpoints, before
                    # assigning these brackets back to the original basis.
                    for removed_index in range(removed_count):
                        ref removed = exact[removed_index]
                        while _compare(low, removed.low, budget) <= 0 and _compare(removed.low, high, budget) <= 0 and low != high:
                            var middle = _midpoint(low, high, budget)
                            var sign = _kevaluate_sign(remaining, middle, field, budget)
                            if sign == 0:
                                low = middle.copy()
                                high = middle^
                            elif sign == _kevaluate_sign(remaining, low, field, budget):
                                low = middle^
                            else:
                                high = middle^
                    _allocation[_FinalRoot](1, budget)
                    exact.append(_FinalRoot(owner, low^, high^))
                return exact^
            break
    return None


def _refine_relative(mut root: _FinalRoot, basis: List[List[DensePoly[QQ]]], mut field: _SampleField, mut budget: Budget) raises PolyError:
    if root.low == root.high:
        return
    var middle = _midpoint(root.low, root.high, budget)
    var sign = _kevaluate_sign(basis[root.basis], middle, field, budget)
    if sign == 0:
        root.low = middle.copy()
        root.high = middle^
    elif sign == _kevaluate_sign(basis[root.basis], root.low, field, budget):
        root.low = middle^
    else:
        root.high = middle^


def _interval_field_roots(basis: List[List[DensePoly[QQ]]], mut field: _SampleField, mut stats: _IntervalRootStats, mut budget: Budget, *, cap: UInt64 = UInt64(2000000)) raises PolyError -> Optional[List[_FinalRoot]]:
    var available = budget.remaining()
    var allowance = min(cap, available)
    var local = Budget(allowance)
    try:
        var roots = List[_FinalRoot]()
        for k in range(len(basis)):
            var found = _one_basis(basis[k], k, field, stats, local)
            if not found:
                budget.charge(allowance - local.remaining())
                return None
            for root in found.value():
                _allocation[_FinalRoot](1, local)
                roots.append(root.copy())
        for i in range(1, len(roots)):
            var j = i
            while j:
                while _compare(roots[j - 1].high, roots[j].low, local) >= 0 and _compare(roots[j].high, roots[j - 1].low, local) >= 0:
                    if roots[j - 1].low == roots[j - 1].high and roots[j].low == roots[j].high:
                        raise PolyError(PolyErrorKind.DOMAIN, "duplicate roots in interval basis")
                    _refine_relative(roots[j - 1], basis, field, local)
                    _refine_relative(roots[j], basis, field, local)
                if _compare(roots[j - 1].high, roots[j].low, local) < 0:
                    break
                var temporary = roots[j - 1].copy()
                roots[j - 1] = roots[j].copy()
                roots[j] = temporary^
                j -= 1
        budget.charge(allowance - local.remaining())
        return roots^
    except e:
        budget.charge(allowance - local.remaining())
        if e.kind == PolyErrorKind.LIMIT and cap < available:
            return None
        raise e^
