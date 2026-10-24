"""Resultant-cofactor SIMD against widened arithmetic and exact identities."""
from polynomials import Budget, MonomialOrder, resultant_value
from polynomials._modular_simd import _Prime32, _inverse32, _subtract32
from polynomials._resultant_cofactors import _interpolate, _interpolate32, _bivariate_resultant
from polynomials._unit_bezout import _modular_bezout
from support.check import Checker
from support.rng import SplitMix64
from support.corpus import parse_polynomial


def multiply(a: UInt64, b: UInt64, p: UInt64) -> UInt64:
    return UInt64((UInt128(a) * UInt128(b)) % UInt128(p))


def power(a: UInt64, exponent: UInt64, p: UInt64) -> UInt64:
    var r = UInt64(1)
    var b = a
    var e = exponent
    while e:
        if e & 1:
            r = multiply(r, b, p)
        b = multiply(b, b, p)
        e >>= 1
    return r


def products[W: Int](mut c: Checker):
    var rng = SplitMix64(198721)
    for p in [UInt64(2), UInt64(3), UInt64(65521), UInt64(1073741823), UInt64(1073741827), UInt64(2147483647), UInt64(2147483659), UInt64(4294967291), UInt64(4294967295)]:
        var field = _Prime32.of(p)
        for round in range(128):
            var a = SIMD[DType.uint64, W](0)
            var b = a
            for lane in range(W):
                a[lane] = p - 1 if round == 0 else rng.next() % p
                b[lane] = p - 1 if round == 0 else rng.next() % p
                if round == 1:
                    a[lane] = 0
            var got = field.mul[W](a, b)
            for lane in range(W):
                c.check(got[lane] == multiply(a[lane], b[lane], p), "32-bit reciprocal reduction covers full unsigned prime boundary")


def small_arrays(mut c: Checker):
    var rng = SplitMix64(8713)
    for p in [UInt64(2), UInt64(65521), UInt64(2147483659), UInt64(4294967291)]:
        var field = _Prime32.of(p)
        for n in [0, 1, 7, 8, 9, 15, 16, 17, 65, 257]:
            var values = List[UInt64]()
            for _ in range(n):
                values.append(1 + rng.next() % (p - 1))
            var scalar = values.copy()
            var vector = values.copy()
            _inverse32[False](scalar, field)
            _inverse32[True](vector, field)
            c.check(vector == scalar, "batched inversions agree with scalar reciprocal arithmetic")
            for i in range(n):
                c.check(multiply(values[i], vector[i], p) == 1, "each batched inverse has exact multiplicative identity")
            for alias in range(2):
                var a = List[UInt64](length=n + 2, fill=UInt64.MAX)
                var b = a.copy()
                for i in range(n):
                    a[i + 1] = rng.next() % p
                    b[i + 1] = a[i + 1] if alias else rng.next() % p
                var old = a.copy()
                var source = b.unsafe_ptr().unsafe_offset(1) if alias == 0 else a.unsafe_ptr().unsafe_offset(1)
                _subtract32(a.unsafe_ptr().unsafe_offset(1).unsafe_origin_cast[MutUntrackedOrigin](), source, n, p - 1, field)
                c.check(a[0] == UInt64.MAX and a[n + 1] == UInt64.MAX, "small subtract guards")
                for i in range(n):
                    c.check(a[i + 1] == (old[i + 1] + p - multiply(b[i + 1], p - 1, p)) % p, "small subtract agrees with UInt128 including alias and tail")


def interpolation(mut c: Checker):
    var rng = SplitMix64(98351)
    for p in [UInt64(65521), UInt64(2147483659), UInt64(4294967291), UInt64(2305843009213693951)]:
        for n in [1, 2, 7, 8, 9, 15, 16, 17, 33]:
            var coefficients = List[UInt64]()
            var points = List[UInt64]()
            var values = List[UInt64]()
            for i in range(n):
                coefficients.append(rng.next() % p if i % 3 else UInt64(0))
                points.append(UInt64(3 * i + 2))
            for x in points:
                var y = UInt64(0)
                for i in range(n - 1, -1, -1):
                    y = (multiply(y, x, p) + coefficients[i]) % p
                values.append(y)
            var inverses = List[List[UInt64]]()
            for level in range(1, n):
                var row = List[UInt64]()
                for k in range(level, n):
                    row.append(power((points[k] + p - points[k - level]) % p, p - 2, p))
                inverses.append(row^)
            c.check(_interpolate(points, values, inverses, p) == coefficients, "interpolation recovers independent coefficients at nonconsecutive points, including wide fallback")
            if p < (UInt64(1) << 32):
                c.check(_interpolate32[False](points, values, inverses, _Prime32.of(p)) == coefficients, "scalar reciprocal interpolation control")


def bezout(mut c: Checker):
    var rng = SplitMix64(71981)
    for p in [UInt64(65521), UInt64(2147483659), UInt64(4294967291), UInt64(2305843009213693951)]:
        for n in [2, 7, 8, 9, 17, 32]:
            var f = List[UInt64](length=n + 1, fill=0)
            var g = f.copy()
            for i in range(n):
                f[i] = rng.next() % p
                g[i] = rng.next() % p
            f[n] = 1
            g[n] = 1
            var pair = _modular_bezout(f, g, p)
            c.check(Bool(pair), "deterministic coprime modular pair has Bezout cofactors")
            if pair:
                ref u = pair.value()[0]
                ref v = pair.value()[1]
                c.check(len(u) <= n and len(v) <= n, "Bezout degree bounds")
                var product = List[UInt64](length=2 * n + 1, fill=0)
                for i in range(len(u)):
                    for j in range(len(f)):
                        product[i + j] = (product[i + j] + multiply(u[i], f[j], p)) % p
                for i in range(len(v)):
                    for j in range(len(g)):
                        product[i + j] = (product[i + j] + multiply(v[i], g[j], p)) % p
                for i in range(len(product)):
                    c.check(product[i] == (UInt64(1) if i == 0 else UInt64(0)), "modular Bezout identity by independent schoolbook product")
            c.check(not _modular_bezout(f, f, p), "common factor is rejected")


def cofactors(mut c: Checker) raises:
    var work = Budget()
    var order = MonomialOrder()
    for f_text in ["x1*x0^9 + x0*x1 + 1", "x0^8 + x0*x1^2 + x1 - 1", "x0^2 + x1^2 - 1"]:
        var f = parse_polynomial(f_text, order, work)
        var g = parse_polynomial("x0^7 + x0*x1 + x1 - 1", order, work)
        var result = _bivariate_resultant(f, g, 0, [Int64(0), Int64(1)], work)
        c.check(Bool(result), "bivariate producer survives degree loss and skipped evaluation points")
        if result:
            ref r = result.value()
            c.check(r.resultant == resultant_value(f, g, 0, work), "cofactor producer agrees with resultant-value route")
            c.check(r.cofactors[0].mul(f, work).add(r.cofactors[1].mul(g, work), work) == r.resultant, "exact rational cofactor identity")
    var f = parse_polynomial("x0 + x1", order, work)
    c.check(not _bivariate_resultant(f, f, 0, [Int64(0), Int64(1)], work), "zero resultant remains a fallback")


def main() raises:
    var c = Checker("resultant-cofactor SIMD")
    products[1](c)
    products[2](c)
    products[4](c)
    products[8](c)
    small_arrays(c)
    interpolation(c)
    if c.failures:
        c.finish()
    bezout(c)
    cofactors(c)
    c.finish()
