"""Compare separate scalar, packed scalar and paired SIMD Taylor shifts.

Arguments: MODE DEGREE JOBS REPETITIONS. Packing, allocation and checksums
are outside the kernel timer; complete CAD benchmarks measure those costs.
Input coefficients are bounded so every Taylor-shift intermediate fits Int64.
"""
from polynomials._cad_native import _shift_native
from std.sys import argv
from std.time import perf_counter_ns


def _separate(mut low: List[Int64], mut high: List[Int64]):
    var n = len(low) - 1
    for i in range(n):
        for j in range(n - 1, i - 1, -1):
            low[j] += low[j + 1]
            high[j] += high[j + 1]


def main() raises:
    var args = argv()
    if len(args) != 5:
        raise Error("expected MODE DEGREE JOBS REPETITIONS")
    var mode = String(args[1])
    var degree = Int(args[2])
    var jobs = Int(args[3])
    var repetitions = Int(args[4])
    if degree < 0 or degree > 32 or jobs <= 0 or repetitions <= 0:
        raise Error("invalid native transform dimensions")
    if mode != "separate" and mode != "packed-scalar" and mode != "simd":
        raise Error("invalid native transform mode")
    for repetition in range(repetitions):
        var lows = List[List[Int64]]()
        var highs = List[List[Int64]]()
        var pairs = List[List[SIMD[DType.int64, 2]]]()
        for job in range(jobs):
            var low = List[Int64](capacity=degree + 1)
            var high = List[Int64](capacity=degree + 1)
            var paired = List[SIMD[DType.int64, 2]](capacity=degree + 1)
            for i in range(degree + 1):
                var lo = Int64((job * 7 + i * 3) % 5 - 2)
                var hi = lo + 1
                if mode == "separate":
                    low.append(lo)
                    high.append(hi)
                else:
                    paired.append(SIMD[DType.int64, 2](lo, hi))
            if mode == "separate":
                lows.append(low^)
                highs.append(high^)
            else:
                pairs.append(paired^)
        var start = perf_counter_ns()
        if mode == "separate":
            for job in range(jobs):
                _separate(lows[job], highs[job])
        elif mode == "packed-scalar":
            for job in range(jobs):
                _shift_native[False](pairs[job])
        else:
            for job in range(jobs):
                _shift_native[True](pairs[job])
        var elapsed = perf_counter_ns() - start
        # Keep every transformed coefficient observable, outside the timer.
        var checksum = UInt64(0)
        for job in range(jobs):
            for i in range(degree + 1):
                if mode == "separate":
                    checksum ^= UInt64(lows[job][i]) ^ (UInt64(highs[job][i]) << 1)
                else:
                    checksum ^= UInt64(pairs[job][i][0]) ^ (UInt64(pairs[job][i][1]) << 1)
        print(String('{"mode":"', mode, '","degree":', degree,
            ',"jobs":', jobs, ',"repetition":', repetition,
            ',"elapsed_ns":', elapsed, ',"checksum":', checksum, '}'))
