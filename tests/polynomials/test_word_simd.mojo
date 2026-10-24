"""Native SIMD products, carries, transforms, and FGLM updates against exact oracles."""
from apn_mojo import Integer
from polynomials._word_simd import _mul_wide, _mul_high, _sum_products, _column_axpy
from polynomials._ntt import _Tables, _forward, _inverse, _convolve, _PRIMES, _ROOTS, _PRIME_COUNT
from polynomials._nmod import _field, _to, _from, _dot
from support.check import Checker
from support.rng import SplitMix64


def products[W: Int](mut c: Checker):
    var rng = SplitMix64(19271)
    for pattern in range(260):
        var a = SIMD[DType.uint64, W](0)
        var b = SIMD[DType.uint64, W](0)
        for lane in range(W):
            a[lane] = rng.next()
            b[lane] = rng.next()
            if pattern == 0:
                a[lane] = UInt64.MAX
                b[lane] = UInt64.MAX
            elif pattern == 1:
                a[lane] = UInt64(1) << UInt64((lane * 13) % 64)
                b[lane] = UInt64.MAX
            elif pattern == 2:
                a[lane] = 0xFFFFFFFF
                b[lane] = 0x100000001
            elif pattern == 3:
                a[lane] = 0
        var pair = _mul_wide[W](a, b)
        var high = _mul_high[W](a, b)
        for lane in range(W):
            var expected = UInt128(a[lane]) * UInt128(b[lane])
            c.check(pair[0][lane] == UInt64(expected) and pair[1][lane] == UInt64(expected >> 64), "full-width unsigned product including cross-term carries")
            c.check(high[lane] == UInt64(expected >> 64), "high-only product agrees with UInt128")


def sums(mut c: Checker) raises:
    var rng = SplitMix64(9981)
    var bound = (UInt64(1) << 62) - 1
    for n in [0, 1, 7, 8, 9, 31, 32, 33, 127, 128, 129, 255, 256, 257, 1023, 8193]:
        for pattern in range(3):
            var a = List[UInt64](length=n + 2, fill=0)
            var b = List[UInt64](length=n + 2, fill=0)
            var expected = Integer(0)
            for i in range(n):
                a[i + 1] = bound if pattern == 0 else rng.next() & bound
                b[i + 1] = bound if pattern == 0 else rng.next() & bound
                if pattern == 2 and i % 2 == 0:
                    a[i + 1] = 0
                expected += Integer(a[i + 1]) * Integer(b[i + 1])
            var result = _sum_products(a.unsafe_ptr().unsafe_offset(1), b.unsafe_ptr().unsafe_offset(1), n)
            var actual = Integer(UInt64(result[0])) + (Integer(UInt64(result[0] >> 64)) << 64) + (Integer(result[1]) << 128)
            c.check(actual == expected, "unaligned dot-product sum and 128-bit overflow carries agree with arbitrary precision")


def modular_dots(mut c: Checker) raises:
    var rng = SplitMix64(8711)
    for p in [UInt64(3), UInt64(65521), UInt64(268435399), UInt64(536870909), UInt64(2147483647), UInt64(4294967291), UInt64(2305843009213693951), UInt64(4611686018427387847)]:
        var field = _field(p)
        for n in [0, 1, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257, 1023]:
            for pattern in range(3):
                var a = List[UInt64](length=n + 2, fill=0)
                var b = List[UInt64](length=n + 2, fill=0)
                var expected = UInt128(0)
                for i in range(n):
                    var x = p - 1 if pattern == 0 else rng.next() % p
                    var y = p - 1 if pattern == 0 else rng.next() % p
                    if pattern == 2 and i % 2 == 0:
                        x = 0
                    expected = (expected + UInt128(x) * UInt128(y)) % UInt128(p)
                    a[i + 1] = _to(x, field)
                    b[i + 1] = _to(y, field)
                var pa = a.unsafe_ptr().unsafe_offset(1)
                var pb = b.unsafe_ptr().unsafe_offset(1)
                var scalar = _dot[False](pa, pb, n, field)
                var vector = _dot[True](pa, pb, n, field)
                c.check(vector == scalar and _from(vector, field) == UInt64(expected), "dot products agree with scalar control and modular oracle across prime and length boundaries")
                var square = _dot[True](pa, pa, n, field)
                c.check(square == _dot[False](pa, pa, n, field), "dot inputs may alias")


def columns(mut c: Checker):
    var rng = SplitMix64(9121)
    for p in [UInt64(2), UInt64(101), UInt64(65521), UInt64(2147483647), UInt64(2147483659), UInt64(3037000499)]:
        var square = p * p
        for n in [0, 1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 127, 128, 129, 257]:
            for pattern in range(3):
                var a = List[UInt64](length=n + 2, fill=UInt64.MAX)
                var column = List[UInt32](length=n + 2, fill=UInt32.MAX)
                var factor = p - 1 if pattern == 0 else rng.next() % p
                var expected = List[UInt64](capacity=n)
                for i in range(n):
                    a[i + 1] = square - 1 if pattern == 0 else rng.next() % square
                    column[i + 1] = UInt32(p - 1 if pattern == 0 else rng.next() % p)
                    expected.append(UInt64((UInt128(a[i + 1]) + UInt128(factor) * UInt128(column[i + 1])) % UInt128(square)))
                var scalar = a.copy()
                _column_axpy[True](a.unsafe_ptr().unsafe_offset(1), column.unsafe_ptr().unsafe_offset(1), n, factor, square)
                _column_axpy[False](scalar.unsafe_ptr().unsafe_offset(1), column.unsafe_ptr().unsafe_offset(1), n, factor, square)
                c.check(a == scalar, "FGLM SIMD and scalar updates agree")
                c.check(a[0] == UInt64.MAX and a[n + 1] == UInt64.MAX, "FGLM update preserves both unaligned-buffer guards")
                for i in range(n):
                    c.check(a[i + 1] == expected[i], "FGLM update agrees with widened arithmetic at the prime boundary")


def transforms(mut c: Checker):
    var primes = materialize[_PRIMES]()
    var roots = materialize[_ROOTS]()
    var rng = SplitMix64(6531)
    for index in [0, _PRIME_COUNT // 2, _PRIME_COUNT - 1]:
        var q = primes[index]
        for n in [1, 2, 4, 8, 16, 32, 64, 128, 256, 1024]:
            var tables = _Tables(q, roots[index], n)
            for pattern in range(3):
                var input = List[UInt64](length=n, fill=0)
                for i in range(n):
                    input[i] = 2 * q - 1 if pattern == 0 else rng.next() % (2 * q)
                    if pattern == 2 and i % 2 == 0:
                        input[i] = 0
                var scalar = input.copy()
                var vector = input.copy()
                _forward[False](scalar, tables, q)
                _forward[True](vector, tables, q)
                c.check(vector == scalar, "NTT forward preserves exact lazy residues, including narrow stages")
                _inverse[False](scalar, tables, q)
                _inverse[True](vector, tables, q)
                c.check(vector == scalar, "NTT inverse preserves exact lazy residues")
                for i in range(n):
                    c.check(vector[i] < 2 * q and vector[i] % q == UInt64((UInt128(n) * UInt128(input[i])) % UInt128(q)), "NTT round trip equals n times the input and respects its range")
        for n in [8, 16, 32, 64]:
            var a = List[UInt64](capacity=n)
            var b = List[UInt64](capacity=n)
            for _ in range(n):
                a.append(rng.next() % (2 * q))
                b.append(rng.next() % (2 * q))
            var expected = List[UInt64](length=n, fill=0)
            for i in range(n):
                for j in range(n):
                    var k = (i + j) % n
                    expected[k] = UInt64((UInt128(expected[k]) + UInt128(a[i]) * UInt128(b[j])) % UInt128(q))
            var actual = _convolve(a.copy(), b.copy(), n, False, index, n)
            c.check(actual == expected, "cyclic convolution equals independent schoolbook modular arithmetic")


def main() raises:
    var c = Checker("native-word SIMD")
    products[1](c)
    products[2](c)
    products[4](c)
    products[8](c)
    sums(c)
    modular_dots(c)
    columns(c)
    transforms(c)
    c.finish()
