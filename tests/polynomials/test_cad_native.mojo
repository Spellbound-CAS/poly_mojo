"""Paired native transforms agree with independent binomial arithmetic."""
from apn_mojo import Integer
from polynomials import Budget, PolyErrorKind
from polynomials._cad_native import _shift_native
from polynomials._cad_interval_roots import _descartes, _IntervalRootStats
from support.check import Checker


def _binomial(n: Int, k: Int) raises -> Integer:
    var value = Integer(1)
    for i in range(k):
        value = value * Integer(n - i) // Integer(i + 1)
    return value^


def main() raises:
    var c = Checker("native CAD transforms")
    for degree in [0, 1, 2, 3, 4, 8, 16, 31, 62, 63]:
        for pattern in range(4):
            var pairs = List[SIMD[DType.int64, 2]]()
            var low = List[Integer]()
            var high = List[Integer]()
            for i in range(degree + 1):
                var lo = Int64(0)
                var hi = Int64(0)
                if degree < 62:
                    var scale = Int64(1) << Int64(max(0, 59 - degree))
                    lo = Int64((i * 7 + pattern * 3) % 5 - 2) * scale
                    hi = lo + Int64(pattern % 2) * scale
                elif degree == 62:
                    lo = -1
                    hi = 1
                pairs.append(SIMD[DType.int64, 2](lo, hi))
                low.append(Integer(lo))
                high.append(Integer(hi))
            var scalar = pairs.copy()
            var vector = pairs.copy()
            _shift_native[False](scalar)
            _shift_native[True](vector)
            for k in range(degree + 1):
                var expected_low = Integer(0)
                var expected_high = Integer(0)
                for j in range(k, degree + 1):
                    var choose = _binomial(j, k)
                    expected_low += low[j] * choose
                    expected_high += high[j] * choose
                c.check(Integer(scalar[k][0]) == expected_low and Integer(vector[k][0]) == expected_low, "lower transform equals exact binomial sum")
                c.check(Integer(scalar[k][1]) == expected_high and Integer(vector[k][1]) == expected_high, "upper transform equals exact binomial sum")

            var native_stats = _IntervalRootStats()
            var default_stats = _IntervalRootStats()
            var wide_stats = _IntervalRootStats()
            var native_work = Budget()
            var default_work = Budget()
            var wide_work = Budget()
            var native = _descartes[True](low, high, True, native_stats, native_work)
            var original = _descartes(low, high, True, default_stats, default_work)
            var wide = _descartes(low, high, False, wide_stats, wide_work)
            c.check(Bool(native) and Bool(wide) and native.value() == wide.value(), "native variation including uncertain signs agrees with arbitrary precision")
            c.check(Bool(original) and original.value() == native.value() and default_work.remaining() == native_work.remaining(), "default scalar and opt-in paired paths agree and charge the same completed work")
            c.check(native_stats.native_nodes == 1 and wide_stats.wide_nodes == 1, "both transform paths were exercised")
            c.check(UInt64.MAX - native_work.remaining() == UInt64(3 * (degree + 1) + degree * (degree + 1)), "successful native transform preserves coefficient work charge")

    for value in [Integer(Int64.MIN), Integer(1) << 63, -(Integer(1) << 90)]:
        var stats = _IntervalRootStats()
        var work = Budget()
        var low: List[Integer] = [value.copy(), Integer(1)]
        var high = low.copy()
        c.check(not _descartes[True](low, high, True, stats, work), "unrepresentable or unsafe paired transform requests wider arithmetic")
        c.check(not _descartes(low, high, True, stats, work), "unrepresentable or unsafe default transform requests wider arithmetic")
        c.check(stats.native_nodes == 0, "overflow never enters the native kernel")
    # At degree one, 62 input bits fit the 63-bit bound; 63 do not.
    for bits in [61, 62]:
        var stats = _IntervalRootStats()
        var work = Budget()
        var low: List[Integer] = [Integer(1) << bits, Integer(1)]
        var high = low.copy()
        c.check(Bool(_descartes[True](low, high, True, stats, work)) == (bits == 61), "exact native overflow-bound boundary")
    var work = Budget(3 * 3 + 2 * 3 - 1)
    var stats = _IntervalRootStats()
    var low: List[Integer] = [Integer(-1), Integer(0), Integer(1)]
    var high = low.copy()
    var limited = False
    try:
        _ = _descartes[True](low, high, True, stats, work)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    c.check(limited and stats.native_nodes == 0, "insufficient work budget stops before the native transform")
    c.finish()
