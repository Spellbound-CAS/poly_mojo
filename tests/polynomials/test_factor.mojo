"""Exact integer/rational factors, Hensel identities, and corrupted products."""
from apn_mojo import Rational, Integer
from polynomials import Budget, PolyErrorKind, CoefficientRing, QQ, ZZ, NoContext, DensePoly, MonomialOrder, factor, factor_integer, check_factorization, Factorization, Modulus, factor_mod_p, prime
from polynomials.factor import _integer, _rational, _bound, _divides_word
from polynomials._hensel import _hensel, _hensel_start, _hensel_extend, _truncate, _divmod_monic, _long_divmod_monic, _build_tree, _image_inverse, _raise_inverse, _lift_image, _low_product, _reverse, _low, _block_divmod_with
from polynomials._dense_ops import _monic, _one
from polynomials._gcd import _mod_image
from support.corpus import read_cases, parse_polynomial, as_zz
from support.check import Checker
from support.rng import SplitMix64
from polynomials._kronecker import _kronecker_mul


def boundaries[C: CoefficientRing](mut c: Checker) raises:
    var work = Budget()
    for value in range(-2, 3):
        var qq = DensePoly[QQ].from_coefficients(19, NoContext(), [QQ(Rational(value))], work)
        comptime if C == QQ:
            var got = factor(qq, work)
            c.check(check_factorization(qq, got, work) and len(got.factors) == 0, "rational constants")
        else:
            var zz = DensePoly[ZZ].from_coefficients(19, NoContext(), [ZZ(Integer(value))], work)
            var got = factor_integer(zz, work)
            c.check(check_factorization(zz, got, work) and len(got.factors) == 0, "integer constants")
    var xq = DensePoly[QQ].from_coefficients(19, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], work)
    var f: DensePoly[C]
    comptime if C == QQ:
        f = rebind_var[DensePoly[C]](xq.copy())
    else:
        f = rebind_var[DensePoly[C]](_integer(xq, work))
    var before = f.copy()
    var tiny = Budget(25)
    try:
        comptime if C == QQ:
            _ = factor(rebind[DensePoly[QQ]](f), tiny)
        else:
            _ = factor_integer(rebind[DensePoly[ZZ]](f), tiny)
        c.check(False, "factor budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and f == before, "factor budget preserves input")


def zz(values: List[Int], mut work: Budget) raises -> DensePoly[ZZ]:
    var coefficients = List[ZZ]()
    for v in values:
        coefficients.append(ZZ(Integer(v)))
    return DensePoly[ZZ].from_coefficients(0, NoContext(), coefficients^, work)


def repeated_blocks(mut c: Checker) raises:
    # Integer Yun blocks: multiplicities, content and signs survive exactly.
    var work = Budget()
    var a = zz([1, 0, 1], work)
    var b = zz([3, 2], work)
    var d = zz([7, -1, 0, 0, 0, 1], work)
    var e = zz([-1, 1], work)
    var f = a.pow(3, work).mul(b.pow(2, work), work).mul(d, work).mul(e.pow(4, work), work).scale(ZZ(Integer(-6)), work)
    var got = factor_integer(f, work)
    c.check(check_factorization(f, got, work), "repeated blocks reconstruct")
    c.check(got.unit == ZZ(Integer(-6)), "repeated blocks keep signed content")
    c.check(len(got.factors) == 4, "four distinct repeated blocks")
    for entry in got.factors:
        if entry[0] == a:
            c.check(entry[1] == 3, "cube block")
        elif entry[0] == b:
            c.check(entry[1] == 2, "square block")
        elif entry[0] == d:
            c.check(entry[1] == 1, "simple block")
        else:
            c.check(entry[0] == e and entry[1] == 4, "fourth power block")
    var high = zz([-2, 0, 1], work).pow(7, work)
    var seven = factor_integer(high, work)
    c.check(len(seven.factors) == 1 and seven.factors[0][1] == 7, "seventh power")
    # Blocks of degree 17 and 20 take Yun's exact quotients through transforms.
    var g20 = List[Int](length=21, fill=0)
    g20[0] = 2
    g20[1] = 2
    g20[20] = 1
    var g17 = List[Int](length=18, fill=0)
    g17[0] = -3
    g17[5] = 3
    g17[17] = 1
    var p = zz(g20, work)
    var q = zz(g17, work)
    var big = p.pow(3, work).mul(q.pow(2, work), work)
    var large = factor_integer(big, work)
    c.check(check_factorization(big, large, work) and len(large.factors) == 2, "large repeated blocks reconstruct")
    for entry in large.factors:
        c.check((entry[0] == p and entry[1] == 3) or (entry[0] == q and entry[1] == 2), "large repeated multiplicities")


def binomials(mut c: Checker) raises:
    """The binomials x^n - 1 and x^n + 1 by the cyclotomic shortcut agree with the general algorithm on the same factors times x + 2."""
    var work = Budget()
    for n in [2, 6, 12, 15, 30, 105]:
        for sign in [-1, 1]:
            var values = List[Int](length=n + 1, fill=0)
            values[0] = sign
            values[n] = 1
            var f = zz(values, work)
            var shortcut = factor_integer(f, work)
            var general = factor_integer(f.mul(zz([2, 1], work), work), work)
            var rebuilt = DensePoly[ZZ].from_coefficients(0, NoContext(), [shortcut.unit.copy()], work)
            var expected = List[DensePoly[ZZ]]()
            for entry in shortcut.factors:
                c.check(entry[1] == 1, "binomial factors are simple")
                rebuilt = rebuilt.mul(entry[0], work)
                expected.append(entry[0].copy())
            c.check(rebuilt == f, "binomial factors rebuild the binomial")
            var others = List[DensePoly[ZZ]]()
            for entry in general.factors:
                if not (entry[0] == zz([2, 1], work)):
                    others.append(entry[0].copy())
            var same = len(others) == len(expected)
            if same:
                for i in range(len(expected)):
                    same = same and others[i] == expected[i]
            c.check(same, "cyclotomic shortcut matches the general algorithm")
    # A repeated binomial reaches the shortcut through its squarefree block.
    var square = zz([-1, 0, 0, 0, 0, 0, 1], work)
    var squared = factor_integer(square.mul(square, work), work)
    var all_double = len(squared.factors) == 4
    for entry in squared.factors:
        all_double = all_double and entry[1] == 2
    c.check(all_double, "(x^6 - 1)^2 has four double cyclotomic factors")


def swinnerton_dyer(mut c: Checker) raises:
    """S4 is irreducible but splits into quadratics modulo every prime; the Zassenhaus constant-term filter must reject every subset, and still accept x^2 + 3 beside it."""
    var work = Budget()
    var s4 = zz([46225, 0, -5596840, 0, 13950764, 0, -7453176, 0, 1513334, 0, -141912, 0, 6476, 0, -136, 0, 1], work)
    var alone = factor_integer(s4, work)
    c.check(len(alone.factors) == 1 and alone.factors[0][0] == s4 and alone.factors[0][1] == 1, "S4 is irreducible")
    var extra = zz([3, 0, 1], work)
    var both = factor_integer(s4.mul(extra, work), work)
    c.check(len(both.factors) == 2 and both.factors[0][0] == extra and both.factors[1][0] == s4, "S4 (x^2 + 3) splits into its two factors")
    # S5 has 16 quadratic modular factors: subsets of two find x^2 + 3 beside
    # it, and the lattice proves the rest irreducible.
    var s5 = zz([2000989041197056, 0, -44660812492570624, 0, 183876928237731840, 0, -255690851718529024, 0, 172580952324702208, 0, -65892492886671360, 0, 15459151516270592, 0, -2349014746136576, 0, 239210760462336, 0, -16665641517056, 0, 801918722048, 0, -26625650688, 0, 602397952, 0, -9028096, 0, 84864, 0, -448, 0, 1], work)
    var s5_alone = factor_integer(s5, work)
    c.check(len(s5_alone.factors) == 1 and s5_alone.factors[0][0] == s5 and s5_alone.factors[0][1] == 1, "S5 is irreducible")
    var s5_both = factor_integer(s5.mul(extra, work), work)
    c.check(len(s5_both.factors) == 2 and s5_both.factors[0][0] == extra and s5_both.factors[1][0] == s5, "S5 (x^2 + 3) splits into its two factors")


def divisibility_filter(mut c: Checker) raises:
    """The word-sized divisor test of recombination at its edges.

    Lifting stops once M exceeds the factor bound, so a true lc(q) h(0) is
    nonzero modulo M in factorization; the predicate itself must still keep
    a zero residue whenever the target reaches M.
    """
    var m = UInt128(1000003) * UInt128(1000003) * UInt128(1000003)
    c.check(not _divides_word(0, m - 1, m), "zero residue below the modulus is no divisor")
    c.check(_divides_word(0, m, m), "zero residue kept once the target reaches the modulus")
    c.check(_divides_word(0, 3 * m, m), "zero residue kept for a multiple of the modulus")
    c.check(_divides_word(m - 5, 10, m) and not _divides_word(m - 5, 4, m), "negative residues by symmetric magnitude")
    c.check(_divides_word(7, 7, m) and not _divides_word(8, 7, m), "positive residues by magnitude")
    c.check(_divides_word(m / 2, m / 2, m) and _divides_word(m / 2 + 1, m / 2, m), "the symmetric boundary")


def kronecker_products(mut c: Checker) raises:
    """Kronecker-packed products match schoolbook on signed operands of mixed sizes."""
    var work = Budget()
    var rng = SplitMix64(4021)
    var all = True
    for la in [1, 2, 3, 5, 9, 17, 31]:
        for lb in [1, 4, 8, 13]:
            for bits in [1, 31, 32, 33, 64, 100, 200, 384]:
                var sides = List[DensePoly[ZZ]]()
                for length in [la, lb]:
                    var coefficients = List[ZZ]()
                    for i in range(length):
                        var v = Integer(0)
                        for _ in range(bits // 32 + 1):
                            v = v * Integer(1 << 32) + Integer(Int(rng.next() >> 32))
                        v = v % (Integer(1) << bits)
                        if i == length - 1 and v.sign() == 0:
                            v = Integer(1)
                        # Extremes: all-negative, alternating and random signs.
                        var negative = (rng.next() & 1) == 1 if bits != 384 else (i % 2 == 0)
                        coefficients.append(ZZ(-v if negative else v))
                    sides.append(DensePoly[ZZ].from_coefficients(0, NoContext(), coefficients^, work))
                var expected = sides[0].mul_schoolbook(sides[1], work)
                all = all and _kronecker_mul(sides[0], sides[1], work) == expected
    c.check(all, "Kronecker products match schoolbook")
    # Borrow chains: every slot at -1 and at the extreme magnitude.
    var minus = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(-1)), ZZ(Integer(-1)), ZZ(Integer(-1))], work)
    c.check(_kronecker_mul(minus, minus, work) == minus.mul_schoolbook(minus, work), "Kronecker borrows through negative slots")
    var top = (Integer(1) << 200) - Integer(1)
    var wide = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(-top), ZZ(top.copy()), ZZ(-top), ZZ(top.copy())], work)
    c.check(_kronecker_mul(wide, wide, work) == wide.mul_schoolbook(wide, work), "Kronecker at full slot magnitude")


def carried_inverses(mut c: Checker) raises:
    """A divisor's series inverse from its image modulo p, raised p-adically,
    inverts its reversal modulo x^n and p^k, and divides as long division does."""
    var work = Budget()
    var rng = SplitMix64(9157)
    var values = List[Int]()
    for _ in range(60):
        values.append(Int(rng.next() % 2001) - 1000)
    values.append(1)
    var h = zz(values, work)
    var p = 10007
    var image = _monic(_mod_image(h, Modulus(UInt64(p), work), work), work)
    var inverse = _lift_image(_image_inverse(image, work), work)
    var inverse_modulus = ZZ(Integer(p))
    var target = ZZ(Integer(p) ** 4)
    _raise_inverse(h, inverse, inverse_modulus, target, work)
    var n = h.degree()
    var check = _truncate(_low_product(_low(_reverse(h, n + 1, work), n, work), inverse, n, work), target, work)
    var one = check.degree() == 0 and check.coefficients[0] == ZZ(Integer(1))
    var dividend = _truncate(h.mul(h, work).add(zz([3, -5, 7], work), work), target, work)
    var hs = _truncate(h, target, work)
    var long = _long_divmod_monic(dividend, hs, target, work)
    var block = _block_divmod_with(dividend, hs, inverse, target, work)
    c.check(one and inverse_modulus == target and long[0] == block[0] and long[1] == block[1], "carried series inverse divides as long division")


def main() raises:
    var work = Budget()
    var c = Checker("factor")
    carried_inverses(c)
    repeated_blocks(c)
    binomials(c)
    swinnerton_dyer(c)
    divisibility_filter(c)
    kronecker_products(c)
    var count = 0
    for record in read_cases():
        var op = record.field("op ")
        if op != "factor" and op != "factor_integer":
            continue
        var text = record.field("in f ")
        var sparse = parse_polynomial(text, MonomialOrder(), work)
        var f = DensePoly[QQ].from_sparse(sparse, work)
        var expected = List[Tuple[DensePoly[QQ], Int]]()
        for line in record.lines:
            if line.startswith("out factor "):
                var parts = line.split(" ")
                var exponent = Int(String(parts[2]))
                var factor_text = String(line[byte=12 + String(exponent).byte_length():])
                expected.append((DensePoly[QQ].from_sparse(parse_polynomial(factor_text, MonomialOrder(), work), work), exponent))
        var unit_text = record.field("out unit ")
        var expected_unit = QQ(Rational(unit_text))
        if op == "factor":
            var got = factor(f, work)
            c.check(got.unit == expected_unit and got.factors == expected, record.name)
            c.check(check_factorization(f, got, work), "rational product certificate")
            got.unit = got.unit.add(QQ.one(NoContext()), work)
            c.check(not check_factorization(f, got, work), "corrupt rational product")
        else:
            var zf = DensePoly[ZZ].from_sparse(as_zz(sparse, work), work)
            var got = factor_integer(zf, work)
            c.check(got.unit.value == expected_unit.value.to_integer_exact() and len(got.factors) == len(expected), record.name)
            for i in range(len(expected)):
                c.check(got.factors[i][1] == expected[i][1] and _rational(got.factors[i][0], work) == expected[i][0], "integer factor")
            c.check(check_factorization(zf, got, work), "integer product certificate")
            got.unit = got.unit.add(ZZ.one(NoContext()), work)
            c.check(not check_factorization(zf, got, work), "corrupt integer product")
        count += 1
    c.check(count == 21, "all twenty-one integer/rational factor cases")
    boundaries[QQ](c)
    boundaries[ZZ](c)
    # Planted primitive irreducibles: linear factors and positive quadratics.
    for n in range(1, 9):
        var a = DensePoly[ZZ].from_coefficients(-3, NoContext(), [ZZ(Integer(-n)), ZZ(Integer(2 * n + 1))], work)
        var b = DensePoly[ZZ].from_coefficients(-3, NoContext(), [ZZ(Integer(n)), ZZ(Integer(0)), ZZ(Integer(1))], work)
        var f = a.pow(2, work).mul(b.pow(UInt64(n % 3 + 1), work), work).scale(ZZ(Integer(-6)), work)
        var got = factor_integer(f, work)
        c.check(got.unit.value == Integer(-6) and len(got.factors) == 2, "planted unit and count")
        c.check(got.factors[0][0] == a and got.factors[0][1] == 2 and got.factors[1][0] == b and got.factors[1][1] == n % 3 + 1, "planted irreducibles and multiplicities")
        c.check(check_factorization(f, got, work), "planted product")
        got.factors[0][1] = 0
        c.check(not check_factorization(f, got, work), "zero exponent refused")
        got.factors[0][1] = -1
        c.check(not check_factorization(f, got, work), "negative exponent refused")
    # Force several base-p correction digits, including a nonmonic input.
    var h = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(35)), ZZ(Integer(31)), ZZ(Integer(6))], work)
    var ctx = Modulus(13, work)
    var fs = factor_mod_p(_monic(_mod_image(h, ctx, work), work), work, 42)
    var lifted = _hensel(h, fs, _bound(h, work), work)
    var product = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(6))], work)
    for g in lifted[0]:
        c.check(g.coefficients[g.degree()].value == Integer(1), "lift preserves monicity")
        product = product.mul(g, work)
    c.check(_truncate(h.sub(product, work), lifted[1], work).is_zero(), "Hensel product modulo lifted prime power")
    # The word-size lift equals the big-integer tree lift. A state with its
    # tree built extends on big integers.
    var reference = _hensel_start(fs, work)
    _build_tree(fs, reference, work)
    _hensel_extend(h, fs, _bound(h, work), reference, work)
    var same = reference.modulus == lifted[1] and len(reference.lifted) == len(lifted[0])
    for i in range(len(lifted[0])):
        same = same and _truncate(reference.lifted[i], reference.modulus, work) == lifted[0][i]
    c.check(same, "word-size lift matches the big-integer lift")
    # A bound between 2^62 and 2^124 takes the two-word Montgomery lift.
    var wide_left = List[ZZ]()
    var wide_right = List[ZZ]()
    for i in range(9):
        wide_left.append(ZZ(Integer((i * 7919 + 13) % 1000003 - 500001) * Integer(1021)))
        wide_right.append(ZZ(Integer((i * 104729 + 7) % 999983 - 499991) * Integer(1031)))
    wide_left[8] = ZZ(Integer(3))
    wide_right[8] = ZZ(Integer(5))
    var wide = DensePoly[ZZ].from_coefficients(0, NoContext(), wide_left^, work).mul(DensePoly[ZZ].from_coefficients(0, NoContext(), wide_right^, work), work)
    var wide_bound = _bound(wide, work)
    var wide_prime = Modulus(268435399, work)
    var wide_images = factor_mod_p(_monic(_mod_image(wide, wide_prime, work), work), work, 7)
    var wide_lift = _hensel(wide, wide_images, wide_bound, work)
    var wide_reference = _hensel_start(wide_images, work)
    _build_tree(wide_images, wide_reference, work)
    _hensel_extend(wide, wide_images, wide_bound, wide_reference, work)
    var wide_same = wide_bound.value.magnitude_bit_length() > 62 and wide_lift[1].value.magnitude_bit_length() < 124 and wide_reference.modulus == wide_lift[1]
    for i in range(len(wide_lift[0])):
        wide_same = wide_same and _truncate(wide_reference.lifted[i], wide_reference.modulus, work) == wide_lift[0][i]
    c.check(wide_same and len(wide_lift[0]) >= 2, "two-word lift matches the big-integer lift")
    # Bounds between 2^124 and 2^254 take the three- and four-limb lifts.
    for bits in [150, 189, 200, 230, 252]:
        var forced = ZZ(Integer(1) << bits)
        var limb_lift = _hensel(wide, wide_images, forced, work)
        var limb_reference = _hensel_start(wide_images, work)
        _build_tree(wide_images, limb_reference, work)
        _hensel_extend(wide, wide_images, forced, limb_reference, work)
        var limb_same = limb_reference.modulus == limb_lift[1] and limb_lift[1].value.magnitude_bit_length() >= 124
        for i in range(len(limb_lift[0])):
            limb_same = limb_same and _truncate(limb_reference.lifted[i], limb_reference.modulus, work) == limb_lift[0][i]
        c.check(limb_same, "multi-limb lift matches the big-integer lift")
    # Extending past 2^254 after a fixed-width lift builds the tree and lifts
    # again from the images, to the big-integer lift.
    var staged = _hensel_start(wide_images, work)
    _hensel_extend(wide, wide_images, ZZ(Integer(1) << 150), staged, work)
    _hensel_extend(wide, wide_images, ZZ(Integer(1) << 300), staged, work)
    var direct = _hensel_start(wide_images, work)
    _build_tree(wide_images, direct, work)
    _hensel_extend(wide, wide_images, ZZ(Integer(1) << 300), direct, work)
    var staged_same = staged.modulus == direct.modulus and staged.modulus.value.magnitude_bit_length() > 254
    for i in range(len(direct.lifted)):
        staged_same = staged_same and _truncate(staged.lifted[i], staged.modulus, work) == _truncate(direct.lifted[i], direct.modulus, work)
    c.check(staged_same, "a fixed-width lift extends to the big-integer lift")
    # A tree lift whose last step left the cofactors pending extends, after
    # lifting them, to the factors one longer lift gives.
    var deferred = _hensel_start(wide_images, work)
    _build_tree(wide_images, deferred, work)
    _hensel_extend(wide, wide_images, ZZ(Integer(1) << 150), deferred, work)
    _hensel_extend(wide, wide_images, ZZ(Integer(1) << 300), deferred, work)
    _hensel_extend(wide, wide_images, ZZ(Integer(1) << 700), deferred, work)
    var whole = _hensel_start(wide_images, work)
    _build_tree(wide_images, whole, work)
    _hensel_extend(wide, wide_images, ZZ(Integer(1) << 700), whole, work)
    var deferred_same = deferred.modulus == whole.modulus
    for i in range(len(whole.lifted)):
        deferred_same = deferred_same and _truncate(deferred.lifted[i], deferred.modulus, work) == _truncate(whole.lifted[i], whole.modulus, work)
    c.check(deferred_same, "deferred cofactors lift to the same factors")
    var first_prime = Integer(Int(prime(0, work)))
    for unlucky in range(2):
        var a = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(1)), ZZ(first_prime.copy() if unlucky == 0 else Integer(1))], work)
        var b = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(2) if unlucky == 0 else first_prime + Integer(1)), ZZ(Integer(1))], work)
        var input = a.mul(b, work)
        var got = factor_integer(input, work)
        c.check(len(got.factors) == 2 and ((got.factors[0][0] == a and got.factors[1][0] == b) or (got.factors[0][0] == b and got.factors[1][0] == a)), "skip first prime with degree drop or repeated modular factor")
        c.check(check_factorization(input, got, work), "unlucky factorization product")
        # Measure a warm deterministic call, then stop one work unit short.
        var measured = Budget()
        _ = factor_integer(input, measured)
        var late = Budget(UInt64.MAX - measured.remaining() - 1)
        var original = input.copy()
        try:
            _ = factor_integer(input, late)
            c.check(False, "late factorization budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT and input == original, "late budget fails closed without changing input")
    # Division by a power series inverse agrees with long division modulo M.
    var division_work = Budget()
    var big_modulus = Integer(1)
    for _ in range(40):
        big_modulus = big_modulus * Integer(1000003)
    var division_modulus = ZZ(big_modulus.copy())
    var state = Integer(12345)
    var dividend_values = List[ZZ]()
    var divisor_values = List[ZZ]()
    for i in range(301):
        state = (state * Integer(6364136223) + Integer(1442695)) % big_modulus
        dividend_values.append(ZZ(state - Integer(1000003)))
        if i < 120:
            divisor_values.append(ZZ(state % Integer(998244353) - Integer(499122176)))
    divisor_values.append(ZZ(Integer(1)))
    var numerator = _truncate(DensePoly[ZZ].from_coefficients(0, NoContext(), dividend_values^, division_work), division_modulus, division_work)
    var monic_divisor = DensePoly[ZZ].from_coefficients(0, NoContext(), divisor_values^, division_work)
    var fast = _divmod_monic(numerator, monic_divisor, division_modulus, division_work)
    var slow = _long_divmod_monic(numerator, monic_divisor, division_modulus, division_work)
    c.check(fast[0] == slow[0] and fast[1] == slow[1] and fast[1].degree() < monic_divisor.degree(), "series-inverse division matches long division")
    c.finish()
