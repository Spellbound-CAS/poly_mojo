"""Scalar/SIMD controls for native FGLM updates, NTT pairs, and modular dots.

Arguments: kind (column/ntt/dot), length, modulus, repetitions, vectorized (0/1).
Allocation, table construction, and final checksum traversal are untimed.
Dot inputs change between calls to prevent loop-invariant call elimination.
"""
from std.sys import argv
from std.time import perf_counter_ns
from polynomials._word_simd import _column_axpy
from polynomials._ntt import _Tables, _forward, _inverse, _PRIMES, _ROOTS
from polynomials._nmod import _field, _dot


def run[vectorized: Bool](kind: String, n: Int, p: UInt64, repeats: Int) raises:
    var checksum = UInt64(0)
    var elapsed = Int(0)
    if kind == "column":
        var square = p * p
        var data = List[UInt64](length=n, fill=square - 1)
        var column = List[UInt32](capacity=n)
        for i in range(n):
            column.append(UInt32(UInt64(i * 13519 + 7919) % p))
        var factor = p - 1
        var start = perf_counter_ns()
        for _ in range(repeats):
            _column_axpy[vectorized](data.unsafe_ptr(), column.unsafe_ptr(), n, factor, square)
        elapsed = Int(perf_counter_ns() - start)
        for value in data:
            checksum += value
    elif kind == "ntt":
        var primes = materialize[_PRIMES]()
        var roots = materialize[_ROOTS]()
        var q = primes[0]
        var tables = _Tables(q, roots[0], n)
        var data = List[UInt64](capacity=n)
        for i in range(n):
            data.append(UInt64(i) * 18446744073709551557 % (2 * q))
        var start = perf_counter_ns()
        for _ in range(repeats):
            _forward[vectorized](data, tables, q)
            _inverse[vectorized](data, tables, q)
        elapsed = Int(perf_counter_ns() - start)
        for value in data:
            checksum += value
    elif kind == "dot":
        var field = _field(p)
        var a = List[UInt64](capacity=32 * n)
        var b = List[UInt64](capacity=32 * n)
        var seed = UInt64(918273645)
        for _ in range(32 * n):
            seed = seed * 6364136223846793005 + 1442695040888963407
            a.append(seed % p)
            seed = seed * 6364136223846793005 + 1442695040888963407
            b.append(seed % p)
        var start = perf_counter_ns()
        for i in range(repeats):
            var at = (i % 32) * n
            a[at] = a[at] + 1 if a[at] + 1 < p else UInt64(0)
            checksum += _dot[vectorized](a.unsafe_ptr().unsafe_offset(at), b.unsafe_ptr().unsafe_offset(at), n, field)
        elapsed = Int(perf_counter_ns() - start)
    else:
        raise Error("unknown native SIMD benchmark kind")
    print(String('{"ns": ', elapsed, ', "checksum": ', checksum, '}'))


def main() raises:
    var args = argv()
    if len(args) != 6:
        raise Error("expected KIND LENGTH MODULUS REPETITIONS VECTORIZED")
    var kind = String(args[1])
    var n = Int(String(args[2]))
    var p = UInt64(Int(String(args[3])))
    var repeats = Int(String(args[4]))
    if n <= 0 or repeats <= 0:
        raise Error("length and repetitions must be positive")
    if kind == "column" and (p < 2 or UInt128(p) * UInt128(p) >= (UInt128(1) << 63)):
        raise Error("column modulus must satisfy 2 <= p and p*p < 2**63")
    if kind == "ntt" and (n & (n - 1)) != 0:
        raise Error("NTT length must be a power of two")
    if String(args[5]) != "0" and String(args[5]) != "1":
        raise Error("vectorized must be 0 or 1")
    if String(args[5]) == "1":
        run[True](kind, n, p, repeats)
    else:
        run[False](kind, n, p, repeats)
