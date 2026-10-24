"""Cofactor kernel controls: KIND LENGTH PRIME REPEATS MODE.

Kinds: subtract, inverse. Modes: 0 is the original scalar arithmetic,
1 is the new scalar reduction, 2 enables SIMD. Allocation and
reciprocal preparation are outside the timer; mutable outputs prevent
loop-invariant elimination. Complete-operation benchmarks include setup.
"""
from std.sys import argv
from std.time import perf_counter_ns
from polynomials._modular_simd import _Prime32, _subtract32, _inverse32


def _old_inverse(a: UInt64, p: UInt64) -> UInt64:
    var result = UInt64(1)
    var base = a
    var e = p - 2
    while e:
        if e & 1:
            result = UInt64((UInt128(result) * UInt128(base)) % UInt128(p))
        base = UInt64((UInt128(base) * UInt128(base)) % UInt128(p))
        e >>= 1
    return result


def run[mode: Int](kind: String, n: Int, p: UInt64, repetitions: Int) raises:
    var data = List[UInt64]()
    var source = List[UInt64]()
    var seed = UInt64(817891)
    for _ in range(n):
        seed = seed * 6364136223846793005 + 1442695040888963407
        data.append(1 + seed % (p - 1))
        seed = seed * 6364136223846793005 + 1442695040888963407
        source.append(seed % p)
    var small = _Prime32.of(p)
    var factor = p - 1
    var elapsed: Int
    if kind == "subtract":
        var start = perf_counter_ns()
        for _ in range(repetitions):
            comptime if mode == 0:
                for i in range(n):
                    var t = UInt64((UInt128(factor) * UInt128(source[i])) % UInt128(p))
                    data[i] = (data[i] + p - t) % p
            else:
                _subtract32[mode == 2](data.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](), source.unsafe_ptr(), n, factor, small)
        elapsed = perf_counter_ns() - start
    elif kind == "inverse":
        var start = perf_counter_ns()
        for _ in range(repetitions):
            comptime if mode == 0:
                for i in range(n):
                    data[i] = _old_inverse(data[i], p)
            else:
                _inverse32[mode == 2](data, small)
        elapsed = perf_counter_ns() - start
    else:
        raise Error("unknown kernel")
    var checksum = UInt64(0)
    for x in data:
        checksum += x
    print(String('{"ns":', elapsed, ',"checksum":', checksum, '}'))


def main() raises:
    var args = argv()
    if len(args) != 6:
        raise Error("expected KIND LENGTH PRIME REPEATS MODE")
    var kind = String(args[1])
    var n = Int(args[2])
    var p = UInt64(Int(args[3]))
    var repetitions = Int(args[4])
    var mode = Int(args[5])
    if n < 1 or repetitions < 1 or p < 2 or p >= (UInt64(1) << 32):
        raise Error("invalid dimensions or modulus")
    if mode == 0:
        run[0](kind, n, p, repetitions)
    elif mode == 1:
        run[1](kind, n, p, repetitions)
    elif mode == 2:
        run[2](kind, n, p, repetitions)
    else:
        raise Error("mode must be 0, 1, or 2")
