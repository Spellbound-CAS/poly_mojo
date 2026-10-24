"""Truncated products across split boundaries and non-real exact series identities."""
from apn_mojo import Rational, ExactComplex
from polynomials import *
from support.check import Checker
from polynomials._toolkit_exact import _exact_product, _ExactOps
from polynomials._toolkit_ops import _divide_values, _convolve_values, _compose_values
from polynomials._series_ops import _exp_values
from polynomials._series_nmod import _native_series, _middle as _nmod_middle
from polynomials._nmod import _field, _to, _mulm, _addm
from polynomials._ntt import _modular_range
from polynomials._integer_product import _zz_mul
from polynomials._fntt import _fconvolve, _fpower, _FPRIMES, _FROOTS, _FCUBE_ROOTS, _FPRIME_COUNT, _FPRIME_BITS
from polynomials._ntt import _transform_length, _ntt_middle_product, _ks_range, _ks_core, _Segments, _digit_width, _zz_range, _bit_length, _PRIME_BITS, _transform_quotient, _shoup_quotient, _PRIMES, _PRIME_COUNT
from polynomials.coefficient_buffer import _view
from apn_mojo import Integer


def value[C: CoefficientRing](ctx: C.Context, k: Int, bits: Int) raises -> C:
    """A deterministic coefficient of about `bits` bits, with signs and denominators varying."""
    var x = Integer(k * k * 7919 + 13 * k + 5)
    var magnitude = x
    while magnitude.magnitude_bit_length() < bits:
        magnitude = magnitude * x + Integer(k + 1)
    var signed = -magnitude if k % 3 == 1 else magnitude
    comptime if C == ZZ:
        return rebind_var[C](ZZ(signed))
    elif C == QQ:
        return rebind_var[C](QQ(Rational(signed, Integer(k % 7 + 1))))
    else:
        var p = rebind[Modulus](ctx).value
        var r = UInt64(k * k * 7919 + 13 * k + 5) % p
        return rebind_var[C](ModInt(p - 1 - r if k % 2 else r, rebind[Modulus](ctx)))


def transform_products[C: CoefficientRing](ctx: C.Context, bits: Int, mut check: Checker) raises:
    """Transform truncated products return the work stack's exact lists; middle products agree with the full product."""
    var work = Budget()
    for sizes in [(64, 64), (100, 70), (130, 257), (300, 300)]:
        var a = List[C]()
        var b = List[C]()
        for k in range(sizes[0]):
            a.append(value[C](ctx, k, bits))
        for k in range(sizes[1]):
            b.append(value[C](ctx, k + 1000, bits))
        var full = DensePoly[C].from_coefficients(0, ctx, a.copy(), work).mul_schoolbook(DensePoly[C].from_coefficients(0, ctx, b.copy(), work), work)
        var whole = sizes[0] + sizes[1] - 1
        for length in [64, sizes[0] - 1, max(sizes[0], sizes[1]), max(sizes[0], sizes[1]) + 5, whole, whole + 9]:
            var got = _exact_product(_view(a), _view(b), length, ctx, work)
            var na = min(length, sizes[0])
            var nb = min(length, sizes[1])
            var n = min(length, na + nb - 1)
            # The reference shape: n coefficients, or the trimmed full product truncated to n.
            var expected = List[C]()
            var head_a = List[C]()
            var head_b = List[C]()
            for i in range(na):
                head_a.append(a[i].copy())
            for i in range(nb):
                head_b.append(b[i].copy())
            var full_low = DensePoly[C].from_coefficients(0, ctx, head_a^, work).mul_schoolbook(DensePoly[C].from_coefficients(0, ctx, head_b^, work), work)
            if n <= max(na, nb):
                for i in range(n):
                    expected.append(full_low.coefficients[i].copy() if i < len(full_low.coefficients) else C.zero(ctx))
            else:
                for i in range(min(n, len(full_low.coefficients))):
                    expected.append(full_low.coefficients[i].copy())
            check.check(got == expected, String("transform truncated product ", sizes[0], "x", sizes[1], " to ", length))
        for bounds in [(0, 64), (sizes[0] - 1, sizes[0] + 60), (whole - 70, whole + 5)]:
            var middle = _ntt_middle_product(_view(a), _view(b), bounds[0], bounds[1], work)
            var valid = Bool(middle)
            if middle:
                ref m = middle.value()
                valid = len(m) == bounds[1] - bounds[0]
                for i in range(len(m)):
                    var j = bounds[0] + i
                    var e = full.coefficients[j].copy() if j < len(full.coefficients) else C.zero(ctx)
                    valid = valid and m[i] == e
            check.check(valid, String("middle product ", bounds[0], " to ", bounds[1]))


def naive_zz(a: List[ZZ], b: List[ZZ], lo: Int, hi: Int) raises -> List[ZZ]:
    var out = List[ZZ]()
    for t in range(lo, hi):
        var acc = Integer(0)
        for i in range(max(0, t - len(b) + 1), min(t + 1, len(a))):
            acc = acc + a[i].value * b[t - i].value
        out.append(ZZ(acc))
    return out^


def segmented(mut check: Checker) raises:
    """Kronecker segmentation matches the schoolbook product on ragged digit counts, signs, zeros and ranges, for every count of float primes."""
    var work = Budget()
    for sizes in [(1, 1), (5, 3), (40, 70), (130, 129), (300, 280)]:
        var a = List[ZZ]()
        var b = List[ZZ]()
        var wide = sizes[0] == 300
        for k in range(sizes[0]):
            # Zero to four limbs (to thirty with `wide`), both signs, all-ones magnitudes.
            var bits = [0, 1, 64, 130, 250][k % 5] * (8 if wide else 1)
            if k % 11 == 5:
                a.append(ZZ((Integer(1) << Integer(bits + 1)) - Integer(1)))
            else:
                a.append(value[ZZ](NoContext(), k, bits) if k % 7 != 3 else ZZ(Integer(0)))
        for k in range(sizes[1]):
            b.append(value[ZZ](NoContext(), k + 500, [190, 3, 70, 128][k % 4] * (8 if wide else 1)))
        var whole = sizes[0] + sizes[1] - 1
        for bounds in [(0, whole), (0, 1), (whole // 2, whole + 3), (whole - 1, whole)]:
            var got = _ks_range(_view(a), _view(b), bounds[0], bounds[1], False, work)
            check.check(Bool(got) and got.value() == naive_zz(a, b, bounds[0], bounds[1]), String("segmented range ", sizes[0], "x", sizes[1], " ", bounds[0], "..", bounds[1]))
        # Every count of float primes, each with its widest digit.
        var la = sizes[0]
        var lb = sizes[1]
        var left = _Segments(_view(a), la, work)
        var right = _Segments(_view(b), lb, work)
        var expected = naive_zz(a, b, 0, whole)
        var agree = True
        comptime for k in range(2, 6):
            var width = _digit_width(k, min(la, lb), left, right)
            var stride = left.digits(width) + right.digits(width) - 1
            var n = _transform_length(whole * stride)
            agree = agree and width >= 1 and _ks_core[ZZ, k](left, right, 0, whole, width, stride, n, False, work) == expected
        check.check(agree, String("segmented cores of two to five primes, ", sizes[0], "x", sizes[1]))
        var square = _ks_range(_view(a), _view(a), 0, 2 * sizes[0] - 1, True, work)
        check.check(Bool(square) and square.value() == naive_zz(a, a, 0, 2 * sizes[0] - 1), "segmented square")


def segmented_columns(mut check: Checker) raises:
    """A side of one digit per coefficient convolves digit column by column, on either side, with signs, zeros, all-ones digits and ranges."""
    var work = Budget()
    for sizes in [(1, 1), (1, 40), (7, 300), (200, 150), (513, 64)]:
        var narrow = List[ZZ]()
        var wide = List[ZZ]()
        for k in range(sizes[0]):
            var bits = [0, 1, 17, 50, 61][k % 5]
            if k % 6 == 2:
                narrow.append(ZZ((Integer(1) << Integer(bits)) - Integer(1)))
            else:
                narrow.append(value[ZZ](NoContext(), k + 900, bits) if bits else ZZ(Integer(0)))
        for k in range(sizes[1]):
            var bits = [1500, 0, 64, 700, 2047][k % 5]
            if k % 7 == 4:
                wide.append(ZZ(-((Integer(1) << Integer(bits)) - Integer(1))))
            else:
                wide.append(value[ZZ](NoContext(), k + 40, bits) if bits else ZZ(Integer(0)))
        var whole = sizes[0] + sizes[1] - 1
        for bounds in [(0, whole), (0, 1), (whole // 2, whole + 3), (whole - 1, whole)]:
            var got = _ks_range(_view(narrow), _view(wide), bounds[0], bounds[1], False, work)
            var flipped = _ks_range(_view(wide), _view(narrow), bounds[0], bounds[1], False, work)
            var expected = naive_zz(narrow, wide, bounds[0], bounds[1])
            check.check(Bool(got) and got.value() == expected and Bool(flipped) and flipped.value() == expected, String("segmented columns ", sizes[0], "x", sizes[1], " ", bounds[0], "..", bounds[1]))
        # Every count of float primes; the narrow side takes one digit at some count.
        var left = _Segments(_view(narrow), sizes[0], work)
        var right = _Segments(_view(wide), sizes[1], work)
        var expected = naive_zz(narrow, wide, 0, whole)
        var agree = True
        var columns = False
        comptime for k in range(2, 6):
            var width = _digit_width(k, min(sizes[0], sizes[1]), left, right)
            var stride = left.digits(width) + right.digits(width) - 1
            columns = columns or (left.digits(width) == 1 and stride > 1)
            var n = _transform_length(whole * stride)
            agree = agree and width >= 1 and _ks_core[ZZ, k](left, right, 0, whole, width, stride, n, False, work) == expected
            agree = agree and _ks_core[ZZ, k](right, left, 0, whole, width, stride, n, False, work) == expected
        check.check(agree and columns, String("segmented column cores of two to five primes, ", sizes[0], "x", sizes[1]))


def exact_bits(k: Int, bits: Int) raises -> ZZ:
    """A coefficient of exactly `bits` bits, either sign."""
    var top = Integer(1) << Integer(bits - 1)
    var v = value[ZZ](NoContext(), k, bits + 8).value
    var m = (-v if v.sign() < 0 else v) % top + top
    return ZZ(-m if k % 3 == 1 else m)


def dispatch_bands(mut check: Checker) raises:
    """ZZ products match the schoolbook on both sides of every dispatch boundary: operand lengths 16, 32 and 256, and product bounds at 2, 3, 5, 6 and 7 float primes and 12, 13, 32 and 48 integer primes."""
    var work = Budget()
    for sizes in [(15, 20), (16, 16), (31, 40), (32, 32), (256, 260), (257, 257)]:
        var shorter = min(sizes[0], sizes[1])
        var log = _bit_length(shorter)
        var agree = True
        for bound in [2 * _FPRIME_BITS, 3 * _FPRIME_BITS, 5 * _FPRIME_BITS, 6 * _FPRIME_BITS, 7 * _FPRIME_BITS, 12 * _PRIME_BITS, 13 * _PRIME_BITS, 32 * _PRIME_BITS, 48 * _PRIME_BITS]:
            # Product bounds of 2 bits + log + 1 on either side of `bound`.
            var below = (bound - log - 1) // 2
            for bits in [below, below + 1]:
                var a = List[ZZ]()
                var b = List[ZZ]()
                for k in range(sizes[0]):
                    a.append(exact_bits(k, bits))
                for k in range(sizes[1]):
                    b.append(exact_bits(k + 1000, bits))
                var whole = sizes[0] + sizes[1] - 1
                var full = _zz_range(_view(a), _view(b), 0, whole, work)
                var low = _zz_range(_view(a), _view(b), 0, shorter, work)
                var expected = naive_zz(a, b, 0, whole)
                agree = agree and Bool(full) and full.value() == expected
                agree = agree and Bool(low) and low.value() == naive_zz(a, b, 0, shorter)
                var square = _zz_range(_view(a), _view(a), 0, 2 * sizes[0] - 1, work)
                agree = agree and Bool(square) and square.value() == naive_zz(a, a, 0, 2 * sizes[0] - 1)
        check.check(agree, String("ZZ products across dispatch bands, ", sizes[0], "x", sizes[1]))


def transform_quotients(mut check: Checker):
    """The division-free Shoup quotient of the twiddle tables matches the 128-bit division for every transform prime."""
    var primes = materialize[_PRIMES]()
    var valid = True
    for i in range(_PRIME_COUNT):
        var q = primes[i]
        var reciprocal = UInt64((UInt128(1) << 125) // UInt128(q))
        var x = UInt64(i + 1)
        for k in range(200):
            x = x * 6364136223846793005 + 1442695040888963407
            var y = (x >> 2) % q
            if k < 4:
                y = [UInt64(0), UInt64(1), q - 1, q - 2][k]
            valid = valid and _transform_quotient(y, q, reciprocal) == _shoup_quotient(y, q)
    check.check(valid, "transform quotients")


def float_prime_table(mut check: Checker):
    """Each float prime is c 2^32 + 1 with 3 | c in (2^49.99, 2^50), largest first, with roots of orders exactly 2^32 and 3 2^32.

    A root of order 2^32 makes every prime factor of p one more than a multiple
    of 2^32, above sqrt(p), so the order check also proves p prime.
    """
    var primes = materialize[_FPRIMES]()
    var roots = materialize[_FROOTS]()
    var cube_roots = materialize[_FCUBE_ROOTS]()
    var valid = True
    for i in range(_FPRIME_COUNT):
        var p = primes[i]
        var c = (p - 1) >> 32
        valid = valid and (c << 32) + 1 == p and c % 3 == 0
        valid = valid and p >= 1118122748098670 and p < (UInt64(1) << 50)
        valid = valid and (i == 0 or primes[i - 1] > p)
        var r = roots[i]
        valid = valid and _fpower(r, UInt64(1) << 32, p) == 1 and _fpower(r, UInt64(1) << 31, p) == p - 1
        var g = cube_roots[i]
        valid = valid and _fpower(g, UInt64(3) << 32, p) == 1 and _fpower(g, UInt64(3) << 31, p) != 1 and _fpower(g, UInt64(1) << 32, p) != 1
        valid = valid and _fpower(g, 3, p) == r
    check.check(valid, "float prime and root tables")


def words_integer(words: Int, seed: UInt64) raises -> Integer:
    """An integer of exactly `words` 64-bit words (zero for none), sign varying with the seed."""
    if words == 0:
        return Integer(0)
    var x = seed
    var v = Integer(1)
    while v.magnitude_bit_length() < 64 * words:
        x = x * 6364136223846793005 + 1442695040888963407
        v = (v << Integer(62)) + Integer(Int64(x >> 2))
    v = v >> Integer(v.magnitude_bit_length() - 64 * words)
    return -v if seed % 3 == 1 else v


def integer_products(mut check: Checker) raises:
    """_zz_mul equals apn's product on both sides of its crossover, for zero, one, squares and lopsided operands."""
    var work = Budget()
    var agree = True
    var sizes = [0, 1, 2, 249, 250, 400, 679, 680, 681, 1000, 2500]
    for i in range(len(sizes)):
        for j in range(i, len(sizes)):
            var a = ZZ(words_integer(sizes[i], UInt64(3 * i + 1)))
            var b = ZZ(words_integer(sizes[j], UInt64(5 * j + 2)))
            agree = agree and _zz_mul(a, b, work).value == a.value * b.value
            agree = agree and _zz_mul(b, a, work).value == a.value * b.value
        var a = ZZ(words_integer(sizes[i], UInt64(7 * i + 4)))
        var same = a.copy()
        agree = agree and _zz_mul(a, a, work).value == a.value * a.value
        agree = agree and _zz_mul(a, same, work).value == a.value * a.value
        for unit in [Integer(1), Integer(-1), Integer(0)]:
            agree = agree and _zz_mul(a, ZZ(unit.copy()), work).value == a.value * unit
    # Lopsided: past the area but under the shorter-operand bound, and over both.
    for shape in [(5000, 100), (5000, 249), (5000, 250), (2000, 300)]:
        var a = ZZ(words_integer(shape[0], 11))
        var b = ZZ(words_integer(shape[1], 13))
        agree = agree and _zz_mul(a, b, work).value == a.value * b.value
    check.check(agree, "segmented integer products")


def three_length(count: Int) -> Int:
    """The least 3 * 2^k >= count with 2^k >= 16."""
    var third = 16
    while 3 * third < count:
        third *= 2
    return 3 * third


def float_convolutions(mut check: Checker):
    """Double-precision transforms of 2^k and 3 * 2^k points reproduce exact cyclic-free convolutions for every float prime, below and above the vector width and the cache block."""
    var primes = materialize[_FPRIMES]()
    for index in range(_FPRIME_COUNT):
        var q = primes[index]
        var valid = True
        for la in [1, 3, 8, 9, 17, 100, 300]:
            for lb in [1, 7, 64]:
                var a = List[UInt64]()
                var b = List[UInt64]()
                var x = UInt64(la * 1000 + lb + index)
                for k in range(la):
                    x = x * 6364136223846793005 + 1442695040888963407
                    a.append(q - 1 if k % 5 == 0 else (x >> 14) % q)
                for k in range(lb):
                    x = x * 6364136223846793005 + 1442695040888963407
                    b.append(q - 1 if k % 3 == 0 else (x >> 14) % q)
                var count = la + lb - 1
                var got = _fconvolve(a, b, count, False, index, _transform_length(count))
                var square = _fconvolve(a, a, 2 * la - 1, True, index, _transform_length(2 * la - 1))
                # Lengths 3 * 2^k: a radix-3 layer over three power-of-two transforms.
                var got3 = _fconvolve(a, b, count, False, index, three_length(count))
                var square3 = _fconvolve(a, a, 2 * la - 1, True, index, three_length(2 * la - 1))
                for k in range(count):
                    var e = UInt128(0)
                    for i in range(max(0, k - lb + 1), min(k, la - 1) + 1):
                        e = (e + UInt128(a[i]) * UInt128(b[k - i])) % UInt128(q)
                    valid = valid and UInt64(e) == got[k] and UInt64(e) == got3[k]
                for k in range(2 * la - 1):
                    var e = UInt128(0)
                    for i in range(max(0, k - la + 1), min(k, la - 1) + 1):
                        e = (e + UInt128(a[i]) * UInt128(a[k - i])) % UInt128(q)
                    valid = valid and UInt64(e) == square[k] and UInt64(e) == square3[k]
        # Past the cache block: stages across the whole array, then block by
        # block; with 3 * 2^k points, inside thirds below and above the block.
        for la in [20000, 70000]:
            var a = List[UInt64]()
            var b = List[UInt64]()
            var x = UInt64(index + 77)
            for k in range(la):
                x = x * 6364136223846793005 + 1442695040888963407
                a.append(q - 1 if k % 5 == 0 else (x >> 14) % q)
            for k in range(5):
                b.append(q - 1 - UInt64(k))
            var count = la + 4
            var got = _fconvolve(a, b, count, False, index, _transform_length(count))
            var got3 = _fconvolve(a, b, count, False, index, three_length(count))
            for k in range(count):
                var e = UInt128(0)
                for i in range(max(0, k - 4), min(k, la - 1) + 1):
                    e = (e + UInt128(a[i]) * UInt128(b[k - i])) % UInt128(q)
                valid = valid and UInt64(e) == got[k] and UInt64(e) == got3[k]
        check.check(valid, String("float convolution modulo ", q))


def native_recurrences(mut check: Checker) raises:
    """Prime-field division, exponential and convolution recurrences on native words satisfy their defining identities."""
    var work = Budget()
    for p in [UInt64(101), (UInt64(1) << 62) - 57]:
        var ctx = Modulus(p, work)
        for n in [31, 32, 33, 100]:
            var a = List[ModInt]()
            var b = List[ModInt]()
            var x = UInt64(n)
            for k in range(n):
                x = x * 6364136223846793005 + 1442695040888963407
                a.append(ModInt(x >> 3, ctx))
                b.append(ModInt(UInt64(1) if k == 0 else p - 1 - (x >> 40) % p, ctx))
            # b (a / b) = a modulo x^n, by schoolbook.
            var q = _divide_values[_ExactOps[ModInt]](_view(a), _view(b), n, ctx, work)
            var back = DensePoly[ModInt].from_coefficients(0, ctx, b.copy(), work).mul_schoolbook(DensePoly[ModInt].from_coefficients(0, ctx, q.copy(), work), work)
            var valid = len(q) == n
            for k in range(n):
                var got = back.coefficients[k].value if k < len(back.coefficients) else UInt64(0)
                valid = valid and got == a[k].value
            check.check(valid, String("native series quotient modulo ", p, ", length ", n))
            # exp(f) with f(0) = 0: g' = f' g modulo x^(n-1).
            if UInt64(n) <= p:
                var f = a.copy()
                f[0] = ModInt(0, ctx)
                var g = _exp_values[_ExactOps[ModInt]](_view(f), n, ctx, work)
                var df = List[ModInt]()
                var dg = List[ModInt]()
                for k in range(1, n):
                    df.append(ModInt(UInt64((UInt128(f[k].value) * UInt128(k)) % UInt128(p)), ctx))
                    dg.append(ModInt(UInt64((UInt128(g[k].value) * UInt128(k)) % UInt128(p)), ctx))
                var rhs = DensePoly[ModInt].from_coefficients(0, ctx, df^, work).mul_schoolbook(DensePoly[ModInt].from_coefficients(0, ctx, g.copy(), work), work)
                var ok = len(g) == n and g[0].value == 1
                for k in range(n - 1):
                    var r = rhs.coefficients[k].value if k < len(rhs.coefficients) else UInt64(0)
                    ok = ok and r == dg[k].value
                check.check(ok, String("native series exponential modulo ", p, ", length ", n))
            var conv = _convolve_values[_ExactOps[ModInt]](_view(a), _view(b), n + 7, ctx, work)
            var full = DensePoly[ModInt].from_coefficients(0, ctx, a.copy(), work).mul_schoolbook(DensePoly[ModInt].from_coefficients(0, ctx, b.copy(), work), work)
            var same = len(conv) == n + 7
            for k in range(n + 7):
                var e = full.coefficients[k].value if k < len(full.coefficients) else UInt64(0)
                same = same and conv[k].value == e
            check.check(same, String("native convolution modulo ", p, ", length ", n))


def wide_newton(mut check: Checker) raises:
    """Wide-prime inverses, exponentials and square roots on both sides of their recurrence crossovers, and middle products by dot products or transforms."""
    var work = Budget()
    var p = (UInt64(1) << 62) - 57
    var ctx = Modulus(p, work)
    for n in [511, 512, 513, 895, 896, 897, 2559, 2560, 2561]:
        var a = List[ModInt]()
        var x = UInt64(n)
        for _ in range(n):
            x = x * 6364136223846793005 + 1442695040888963407
            a.append(ModInt(x >> 2, ctx))
        a[0] = ModInt(1, ctx)
        var series = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, ctx, a.copy(), work), n, work)
        var unit = series.mul(series.inverse(work), work)
        check.check(len(unit.polynomial.coefficients) == 1 and unit.polynomial.coefficients[0].value == 1, String("wide-prime inverse at order ", n))
        a[0] = ModInt(0, ctx)
        var nilpotent = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, ctx, a.copy(), work), n, work)
        var back = nilpotent.exp(work).log(work)
        check.check(back.polynomial == nilpotent.polynomial, String("wide-prime log(exp) at order ", n))
    # Square roots on both sides of the recurrence crossovers: 4096, and 2560
    # for primes from 2^28 below 2^59.
    for q in [p, UInt64(2147483647)]:
        var ring = Modulus(q, work)
        for n in [2559, 2560, 4095, 4096]:
            var a = List[ModInt]()
            var x = UInt64(n) + q
            for _ in range(n):
                x = x * 6364136223846793005 + 1442695040888963407
                a.append(ModInt(x >> 2, ring))
            a[0] = ModInt(1, ring)
            var series = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, ring, a.copy(), work), n, work)
            var root = series.sqrt(work)
            check.check(root.mul(root, work).polynomial == series.polynomial and root.polynomial.coefficients[0].value == 1, String("square root modulo ", q, " at order ", n))
    var f = _field(p)
    for shape in [(100, 100, 50, 100), (300, 200, 150, 300), (40, 900, 30, 600), (512, 256, 256, 512), (2000, 1000, 1000, 2000), (7, 9, 0, 20)]:
        var a = List[UInt64]()
        var b = List[UInt64]()
        var x = UInt64(shape[0] + shape[2])
        for _ in range(shape[0]):
            x = x * 6364136223846793005 + 1442695040888963407
            a.append(_to(x >> 2, f))
        for _ in range(shape[1]):
            x = x * 6364136223846793005 + 1442695040888963407
            b.append(_to(p - 1 if x % 3 == 0 else x >> 2, f))
        var got = _nmod_middle(a, b, shape[2], shape[3], f, work)
        var ok = len(got) == shape[3] - shape[2]
        for k in range(shape[2], shape[3]):
            var e = UInt64(0)
            for i in range(max(0, k - len(b) + 1), min(k + 1, len(a))):
                e = _addm(e, _mulm(a[i], b[k - i], f), f)
            ok = ok and got[k - shape[2]] == e
        check.check(ok, String("wide-prime middle product ", shape[0], "x", shape[1], " ", shape[2], "..", shape[3]))


def quadratic_exponentials(mut check: Checker) raises:
    """Rational exponentials by the word recurrence: word and wide denominators, and numerators too wide for it (Newton), satisfy exp(f) exp(-f) = 1, and log exp(f) = f at short orders."""
    var work = Budget()
    var one = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], work)
    for shape in [(8, Int64(1)), (8, Int64(6)), (20, Int64(2305843009213693951)), (60, Int64(7))]:
        for n in [2, 17, 120, 300]:
            var f = List[QQ]()
            var minus = List[QQ]()
            f.append(QQ(Rational(0)))
            minus.append(QQ(Rational(0)))
            var x = UInt64(n * 31 + shape[0])
            for _ in range(1, n):
                x = x * 6364136223846793005 + 1442695040888963407
                var top = Integer(1) << Integer(shape[0] - 1)
                var v = Integer(Int64(x >> 8)) % (top + top) - top
                f.append(QQ(Rational(v, Integer(shape[1]))))
                minus.append(QQ(Rational(-v, Integer(shape[1]))))
            var a = TruncatedSeries[QQ].from_polynomial(DensePoly[QQ].from_coefficients(0, NoContext(), f^, work), n, work)
            var b = TruncatedSeries[QQ].from_polynomial(DensePoly[QQ].from_coefficients(0, NoContext(), minus^, work), n, work)
            var e = a.exp(work)
            var ok = e.mul(b.exp(work), work).polynomial == one
            if n <= 17:
                # The logarithm of these wide outputs is slow; short orders suffice for it.
                ok = ok and e.log(work).polynomial == a.polynomial
            check.check(ok, String("rational exponential, ", shape[0], "-bit numerators over ", shape[1], ", order ", n))


def quadratic_inverses(mut check: Checker) raises:
    """Rational inverses and logarithms by the word recurrence, on both sides of its order crossover and with constant terms that leave it (Newton): f (1 / f) = 1 and exp(log f) = f."""
    var work = Budget()
    var one = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], work)
    for shape in [(8, Int64(1), Int64(1)), (62, Int64(1), Int64(-1)), (8, Int64(5), Int64(5)), (8, Int64(1), Int64(3)), (63, Int64(1), Int64(1))]:
        for n in [16, 40, 511, 512]:
            if n > 40 and shape[0] > 8:
                continue
            var f = List[QQ]()
            f.append(QQ(Rational(Integer(shape[2]), Integer(shape[1]))))
            var x = UInt64(n * 7 + shape[0])
            for _ in range(1, n):
                x = x * 6364136223846793005 + 1442695040888963407
                var top = Integer(1) << Integer(shape[0] - 1)
                f.append(QQ(Rational(Integer(Int64(x >> 1)) % (top + top) - top, Integer(shape[1]))))
            var a = TruncatedSeries[QQ].from_polynomial(DensePoly[QQ].from_coefficients(0, NoContext(), f^, work), n, work)
            var ok = a.mul(a.inverse(work), work).polynomial == one
            if shape[2] == shape[1] and n <= 40:
                # Constant term one: the logarithm reuses the inverse.
                ok = ok and a.log(work).exp(work).polynomial == a.polynomial
            check.check(ok, String("rational inverse, ", shape[0], "-bit numerators over ", shape[1], ", constant ", shape[2], ", order ", n))


def word_compositions(mut check: Checker) raises:
    """Rational compositions with word-sized scaled outer coefficients (word block sums) and wider ones (integer block sums) match the generic composition."""
    var work = Budget()
    for shape in [(8, Int64(1), Int64(1)), (63, Int64(1), Int64(1)), (70, Int64(1), Int64(1)), (8, Int64(3), Int64(1)), (8, Int64(1), Int64(5))]:
        for n in [16, 17, 40]:
            var f = List[QQ]()
            var g = List[QQ]()
            var x = UInt64(n * 13 + shape[0])
            for i in range(n):
                x = x * 6364136223846793005 + 1442695040888963407
                var top = Integer(1) << Integer(shape[0] - 1)
                var v = (Integer(x) << Integer(64) | Integer(x >> 3)) % (top + top) - top
                f.append(QQ(Rational(v, Integer(shape[1]))))
                x = x * 6364136223846793005 + 1442695040888963407
                g.append(QQ(Rational(Integer(Int64(x >> 57)) - Integer(64), Integer(shape[2]))) if i else QQ(Rational(0)))
            var fast = _compose_values[_ExactOps[QQ]](_view(f), _view(g), n, NoContext(), work)
            var slow = _compose_values[_ExactOps[QQ], False](_view(f), _view(g), n, NoContext(), work)
            var same = len(fast) == len(slow)
            for t in range(min(len(fast), len(slow))):
                same = same and fast[t] == slow[t]
            check.check(same, String("rational composition, ", shape[0], "-bit outer numerators over ", shape[1], ", inner over ", shape[2], ", order ", n))


def routed_series(mut check: Checker) raises:
    """Prime-field inverses and exponentials agree with their identities on both sides of the word/Newton crossovers."""
    var work = Budget()
    var p = (UInt64(1) << 62) - 57
    var ring = Modulus(p, work)
    for n in [1535, 1536, 1537, 3071, 3073]:
        var a = List[ModInt]()
        var x = UInt64(n)
        for k in range(n):
            x = x * 6364136223846793005 + 1442695040888963407
            a.append(ModInt(x >> 3, ring))
        var unit = a.copy()
        unit[0] = ModInt(1, ring)
        var g = _native_series(unit, ring, n, 0, work).value().copy()
        var bu = List[UInt64]()
        var gu = List[UInt64]()
        for k in range(n):
            bu.append(unit[k].value)
            gu.append(g[k].value)
        var product = _modular_range(bu, gu, 0, n, p, work)
        var ok = product[0] == 1
        for k in range(1, n):
            ok = ok and product[k] == 0
        check.check(ok, String("routed series inverse, order ", n))
        var zero = a.copy()
        zero[0] = ModInt(0, ring)
        var e = _native_series(zero, ring, n, 2, work).value().copy()
        var da = List[UInt64]()
        var eu = List[UInt64]()
        for k in range(n):
            eu.append(e[k].value)
            if k:
                da.append(UInt64((UInt128(zero[k].value) * UInt128(k)) % UInt128(p)))
        var rhs = _modular_range(da, eu, 0, n - 1, p, work)
        var fine = e[0].value == 1
        for k in range(n - 1):
            fine = fine and rhs[k] == UInt64((UInt128(e[k + 1].value) * UInt128(k + 1)) % UInt128(p))
        check.check(fine, String("routed series exponential, order ", n))


def main() raises:
    var check = Checker("toolkit_products")
    routed_series(check)
    native_recurrences(check)
    wide_newton(check)
    quadratic_exponentials(check)
    quadratic_inverses(check)
    word_compositions(check)
    float_prime_table(check)
    integer_products(check)
    float_convolutions(check)
    transform_quotients(check)
    segmented(check)
    segmented_columns(check)
    dispatch_bands(check)
    var work = Budget()
    transform_products[ZZ](NoContext(), 8, check)
    transform_products[ZZ](NoContext(), 200, check)
    # Eighteen limbs carry the residue sums past 128 bits and take 38 primes.
    transform_products[ZZ](NoContext(), 1100, check)
    # Past the 48-prime table: Kronecker segmentation.
    transform_products[ZZ](NoContext(), 1600, check)
    transform_products[QQ](NoContext(), 40, check)
    transform_products[ModInt](Modulus(65521, work), 0, check)
    transform_products[ModInt](Modulus((UInt64(1) << 62) - 57, work), 0, check)


    var ctx = Modulus(65521, work)
    for size in [63, 64, 65, 127, 128, 129, 257]:
        for short in [size, size - 3, size // 2]:
            var a = List[ModInt]()
            var b = List[ModInt]()
            for k in range(size):
                a.append(ModInt(UInt64((k * k + 3 * k + 1) % 17), ctx))
            for k in range(short):
                b.append(ModInt(UInt64((k * k + 5 * k + 2) % 19), ctx))
            var f = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(-2, ctx, a.copy(), work), size, work)
            var g = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(-2, ctx, b.copy(), work), size, work)
            var got = f.mul(g, work)
            # Independent native-residue oracle, with sums reduced each term.
            var valid = True
            for k in range(size):
                var expected = UInt64(0)
                for i in range(k + 1):
                    if i < len(a) and k - i < len(b):
                        expected = (expected + a[i].value * b[k - i].value) % ctx.value
                var actual = got.polynomial.coefficients[k].value if k < len(got.polynomial.coefficients) else UInt64(0)
                valid = valid and actual == expected
            check.check(valid, "truncated Karatsuba split agrees with modular oracle")
            var saved_f = f.polynomial.copy()
            var saved_g = g.polynomial.copy()
            var tiny = Budget(71)
            try:
                _ = f.mul(g, tiny)
                check.check(False, "truncated product ignored budget")
            except e:
                check.check(e.kind == PolyErrorKind.LIMIT and f.polynomial == saved_f and g.polynomial == saved_g, "product failure preserves operands")
    var imaginary = QQi(ExactComplex(0, 1))
    var long_values = List[ModInt](length=2049, fill=ModInt.one(ctx))
    var long_series = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, ctx, long_values^, work), 2049, work)
    var long_inverse = long_series.inverse(work)
    var long_expected = DensePoly[ModInt].from_coefficients(0, ctx, [ModInt.one(ctx), ModInt.one(ctx).neg(work)], work)
    check.check(long_inverse.polynomial == long_expected, "large modular Newton inverse retains truncation boundary")
    var coefficients = List[QQi](length=65, fill=imaginary)
    var f = TruncatedSeries[QQi].from_polynomial(DensePoly[QQi].from_coefficients(0, NoContext(), coefficients^, work), 65, work)
    var inverse = f.inverse(work)
    var expected = DensePoly[QQi].from_coefficients(0, NoContext(), [QQi(ExactComplex(0, -1)), QQi(ExactComplex(0, 1))], work)
    check.check(inverse.polynomial == expected, "Gaussian rational inverse constant is not real")
    coefficients = List[QQi]()
    coefficients.append(QQi.zero(NoContext()))
    var power = QQi.one(NoContext())
    for _ in range(1, 65):
        coefficients.append(power.copy())
        power = power.mul(imaginary, work)
    f = TruncatedSeries[QQi].from_polynomial(DensePoly[QQi].from_coefficients(0, NoContext(), coefficients^, work), 65, work)
    var reverted = f.revert(work)
    power = QQi.one(NoContext())
    var valid = reverted.polynomial.coefficients[0].is_zero()
    for k in range(1, 65):
        valid = valid and reverted.polynomial.coefficients[k] == power
        power = power.mul(imaginary.neg(work), work)
    check.check(valid, "x/(1-i*x) reverts to x/(1+i*x)")
    check.finish()
