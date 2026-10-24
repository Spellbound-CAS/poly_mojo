"""Closed-form sparse products across rings, orders, collisions and limits."""
from apn_mojo import Integer, Rational
from polynomials import Budget, CoefficientRing, ZZ, QQ, ModInt, Modulus, NoContext, SparsePoly, Term, Monomial, MonomialOrder, PolyErrorKind
from support.products import geometric, expected_product
from support.check import Checker
from support.gen import polynomial
from support.rng import SplitMix64
from polynomials.sparse import _heap_product
from polynomials._packed_products import _integer, _packed_power, _power_recurrence, _packed_poly_product, _packed_poly_power, _packed_poly_divide, _packed_divide
from polynomials.division import divide_exact
from polynomials._packed_poly import _PackedPoly


def products[C: CoefficientRing](base: C, mut c: Checker) raises:
    var work = Budget()
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(1), Int64(0)], work)
        for n in [2, 8, 32, 34]:
            for shape in range(3):
                var a = geometric(n, 0, False, base, order, work)
                var b = geometric(n, Int64(1 if shape == 1 else 0), shape == 2, base, order, work)
                var expected = expected_product(n, shape == 1, shape == 2, base, order, work)
                var product = a.mul(b, work)
                product.check_invariants(work)
                c.check(product == expected and b.mul(a, work) == expected, "closed-form product and reversed operands")
        var a = geometric(2, 0, False, base, order, work)
        var b = geometric(34, 0, False, base, order, work)
        var square = base.mul(base, work)
        var terms = List[Term[C]]()
        for k in range(35):
            var value = square.copy() if k == 0 or k == 34 else square.add(square, work)
            terms.append(Term[C](Monomial([(Int64(0), UInt32(k))], work) if k else Monomial(), value^))
        var expected = SparsePoly[C].from_terms(terms^, order, base.context(), work)
        c.check(a.mul(b, work) == expected and b.mul(a, work) == expected, "asymmetric product")
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0)], work)
    var a = geometric(34, 0, False, base, order, work)
    var b = geometric(34, 0, True, base, order, work)
    var saved_a = a.copy()
    var saved_b = b.copy()
    var expected = expected_product(34, False, True, base, order, work)
    for limit in [UInt64(0), UInt64(100), UInt64(1000), UInt64(10000)]:
        var tiny = Budget(limit)
        try:
            var result = a.mul(b, tiny)
            c.check(result == expected, "limited product is complete")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "multiplication budget fails closed")
        c.check(a == saved_a and b == saved_b, "failed multiplication preserves operands")
    var rng = SplitMix64(118)
    var kernel = MonomialOrder()
    for kind in [MonomialOrder.LEX, MonomialOrder.GRLEX, MonomialOrder.GREVLEX]:
        var a = polynomial[C](base.context(), kernel, rng, work, terms=48, degree=7)
        var b = polynomial[C](base.context(), kernel, rng, work, terms=35, degree=7)
        var ordered = MonomialOrder(kind, [Int64(2), Int64(0), Int64(1)], work)
        var expected = a.mul(b, work).reorder(ordered, work)
        var left = a.reorder(ordered, work)
        var right = b.reorder(ordered, work)
        c.check(left.mul(right, work) == expected and _heap_product(left, right, work) == expected, "changing order commutes with multiplication on generated multivariate inputs")
    var overflow_order = MonomialOrder(MonomialOrder.LEX, [Int64(0)], work)
    var large = geometric(8, 0, False, base, overflow_order, work)
    large.terms[len(large.terms) - 1].monomial.pairs[0] = (Int64(0), UInt32.MAX)
    var saved = large.copy()
    try:
        _ = _heap_product(large, large, work)
        c.check(False, "overflowing monomial product accepted")
    except e:
        c.check(e.kind == PolyErrorKind.CAPACITY, "heap exponent overflow is CAPACITY")
    c.check(large == saved, "exponent overflow preserves operands")


def word[C: CoefficientRing](ctx: C.Context, x: Int64) -> C:
    comptime if C == ZZ:
        return rebind_var[C](ZZ(Integer(x)))
    else:
        var p = rebind[Modulus](ctx).value
        var r = UInt64(x) % p if x >= 0 else p - UInt64(-x) % p
        return rebind_var[C](ModInt(r, rebind[Modulus](ctx)))


def dense_block[C: CoefficientRing](ctx: C.Context, order: MonomialOrder, mut rng: SplitMix64, n: Int, degree: Int, scale: Int, bits: Int, mut work: Budget) raises -> SparsePoly[C]:
    """Random terms of degree below `degree` in each of n variables, exponents scaled by `scale`, coefficients of `bits` bits."""
    var terms = List[Term[C]]()
    var count = 1
    for _ in range(n):
        count *= degree
    for index in range(count):
        if rng.below(5) == 0:
            continue
        var pairs = List[Tuple[Int64, UInt32]]()
        var rest = index
        for i in range(n):
            var e = rest % degree
            rest //= degree
            if e:
                pairs.append((Int64(i), UInt32(e * scale)))
        var magnitude = Int64(rng.next() >> UInt64(64 - bits)) | 1
        var value = -magnitude if rng.below(2) else magnitude
        terms.append(Term[C](Monomial(pairs^, work), word[C](ctx, value)))
    return SparsePoly[C].from_terms(terms^, order, ctx, work)


def packed[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    """Packed-exponent products agree with the heap product, including window splits and wide accumulators."""
    var work = Budget()
    var rng = SplitMix64(4242)
    for n in [1, 2, 3, 4]:
        var variables = List[Int64]()
        for i in range(n):
            variables.append(Int64(i))
        var order = MonomialOrder(MonomialOrder.LEX, variables^, work)
        for shape in range(4):
            # Shapes: small exponents; large exponents that split the window; 63-bit
            # coefficients; and a product with cancellation against its negation.
            var degree = [200, 14, 6, 4][n - 1]
            var scale = 1 if shape != 1 else 97
            var bits = 63 if shape == 2 else 20
            var a = dense_block[C](ctx, order, rng, n, degree, scale, bits, work)
            var b = dense_block[C](ctx, order, rng, n, degree, scale, bits, work)
            if shape == 3:
                b = b.sub(a, work)
            var expected = _heap_product(a, b, work)
            var product = a.mul(b, work)
            product.check_invariants(work)
            c.check(product == expected, String("packed product, ", n, " variables, shape ", shape))
            c.check(b.mul(a, work) == expected, "packed product commutes")
            c.check(a.mul(a, work) == _heap_product(a, a, work), "packed square")
    # Exponents beyond the packed key range use the heap.
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1), Int64(2)], work)
    var a = dense_block[C](ctx, order, rng, 3, 8, 1, 10, work)
    var wide = Monomial([(Int64(0), UInt32(1 << 21)), (Int64(1), UInt32(1 << 21)), (Int64(2), UInt32(1 << 21))], work)
    var terms: List[Term[C]] = [Term[C](wide^, C.one(ctx))]
    a = a.add(SparsePoly[C].from_terms(terms^, order, ctx, work), work)
    c.check(a.mul(a, work) == _heap_product(a, a, work), "wide exponents fall back to the heap")


def packed_poly[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    """Packed storage round-trips, and products of packed operands match the heap product term for term."""
    var work = Budget()
    var rng = SplitMix64(9191)
    for n in [1, 2, 3, 4]:
        var variables = List[Int64]()
        for i in range(n):
            variables.append(Int64(i))
        var order = MonomialOrder(MonomialOrder.LEX, variables^, work)
        for shape in range(4):
            var degree = [200, 14, 6, 4][n - 1]
            var scale = 1 if shape != 1 else 97
            var bits = 63 if shape == 2 else 20
            var a = dense_block[C](ctx, order, rng, n, degree, scale, bits, work)
            var b = dense_block[C](ctx, order, rng, n, degree, scale, bits, work)
            if shape == 3:
                b = b.sub(a, work)
            var pa = _PackedPoly[C].from_sparse(a, work)
            var pb = _PackedPoly[C].from_sparse(b, work)
            c.check(Bool(pa) and Bool(pb), "LEX operands pack")
            c.check(pa.value().to_sparse(work) == a, "packed storage round-trips")
            var product = _packed_poly_product(pa.value(), pb.value(), work)
            var expected = _heap_product(a, b, work)
            c.check(Bool(product) and product.value().to_sparse(work) == expected, String("packed operands multiply, ", n, " variables, shape ", shape))
            var square = _packed_poly_product(pa.value(), pa.value(), work)
            c.check(Bool(square) and square.value().to_sparse(work) == _heap_product(a, a, work), "packed square")
            if product:
                # The product's radix is its exact exponent bound, so packing the expected product reproduces it.
                var repacked = _PackedPoly[C].from_sparse(expected, work)
                c.check(Bool(repacked) and repacked.value() == product.value(), "product keys sit in the product's own radix")
                if shape != 3 and len(b.terms):
                    var quotient = _packed_poly_divide(product.value(), pb.value(), work)
                    if quotient:
                        c.check(Bool(quotient.value()) and quotient.value().value().to_sparse(work) == a, String("packed exact quotient, ", n, " variables, shape ", shape))
                    if quotient and Bool(quotient.value()):
                        c.check(quotient.value().value() == _PackedPoly[C].from_sparse(a, work).value(), "quotient keys sit in the quotient's own radix")
                    var bumped = expected.add(SparsePoly[C].constant(C.one(ctx), order, work), work)
                    var pbumped = _PackedPoly[C].from_sparse(bumped, work)
                    var inexact = _packed_poly_divide(pbumped.value(), pb.value(), work)
                    var reference = divide_exact(bumped, b, work)
                    if inexact:
                        c.check(Bool(inexact.value()) == Bool(reference), "packed inexact division agrees with divide_exact")
    var grevlex = MonomialOrder(MonomialOrder.GREVLEX, [Int64(0), Int64(1)], work)
    var other = dense_block[C](ctx, grevlex, rng, 2, 4, 1, 10, work)
    c.check(not _PackedPoly[C].from_sparse(other, work), "only LEX over variables 0..n-1 packs")


def packed_powers(mut c: Checker) raises:
    """Packed powers equal SparsePoly.pow: the recurrence, repeated squaring, and steps past word coefficients."""
    var work = Budget()
    var specs: List[List[Tuple[Int, Int, Int, Int]]] = [
        [(1, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1)],
        [(3, 0, 0, 0), (-2, 1, 0, 0), (5, 0, 2, 0), (7, 0, 0, 1)],
        [(1, 0, 0, 0), (1, 1, 0, 0), (2, 0, 2, 0), (3, 0, 0, 3)],
        [(1, 0, 0, 0), (97, 5, 0, 0), (-89, 0, 5, 0), (83, 0, 0, 4)],
    ]
    var three = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1), Int64(2)], work)
    var divided = 0
    var remainders = 0
    for spec in specs:
        var f = few_terms(spec, three, work)
        var packed = _PackedPoly[ZZ].from_sparse(f, work)
        for k in [1, 2, 3, 5, 8, 12]:
            var power = _packed_poly_power(packed.value(), k, work)
            c.check(Bool(power) and power.value().to_sparse(work) == f.pow(UInt64(k), work), String("packed power k = ", k))
        # Exact and inexact quotients of powers, against the SparsePoly kernel and divide_exact.
        var big = f.pow(UInt64(10), work)
        var pbig = _PackedPoly[ZZ].from_sparse(big, work)
        var bumped = big.add(SparsePoly[ZZ].constant(ZZ(Integer(1)), three, work), work)
        var pbumped = _PackedPoly[ZZ].from_sparse(bumped, work)
        for j in [1, 3, 4]:
            var g = f.pow(UInt64(j), work)
            var pg = _PackedPoly[ZZ].from_sparse(g, work)
            var quotient = _packed_poly_divide(pbig.value(), pg.value(), work)
            var reference = _packed_divide(big, g, work)
            c.check(Bool(quotient) == Bool(reference), "packed and sparse division take the same inputs")
            if quotient:
                divided += 1
                c.check(Bool(quotient.value()) and quotient.value().value().to_sparse(work) == f.pow(UInt64(10 - j), work), String("packed exact quotient, j = ", j))
                c.check(Bool(quotient.value()) and quotient.value().value() == _PackedPoly[ZZ].from_sparse(f.pow(UInt64(10 - j), work), work).value(), "quotient keys sit in the quotient's own radix")
            var inexact = _packed_poly_divide(pbumped.value(), pg.value(), work)
            if inexact:
                remainders += 1
                c.check(not inexact.value() and not divide_exact(bumped, g, work), "packed division detects a remainder")
    var p = Modulus(101, work)
    var terms = List[Term[ModInt]]()
    terms.append(Term[ModInt](Monomial(), ModInt(3, p)))
    terms.append(Term[ModInt](Monomial([(Int64(0), UInt32(2))], work), ModInt(7, p)))
    terms.append(Term[ModInt](Monomial([(Int64(1), UInt32(1)), (Int64(2), UInt32(3))], work), ModInt(100, p)))
    var g = SparsePoly[ModInt].from_terms(terms^, three, p, work)
    var pg = _PackedPoly[ModInt].from_sparse(g, work)
    for k in [2, 7, 16]:
        var power = _packed_poly_power(pg.value(), k, work)
        c.check(Bool(power) and power.value().to_sparse(work) == g.pow(UInt64(k), work), "packed power over GF(101)")
    var g16 = _PackedPoly[ModInt].from_sparse(g.pow(UInt64(16), work), work)
    var g5 = _PackedPoly[ModInt].from_sparse(g.pow(UInt64(5), work), work)
    var modular = _packed_poly_divide(g16.value(), g5.value(), work)
    if modular:
        divided += 1
        c.check(Bool(modular.value()) and modular.value().value().to_sparse(work) == g.pow(UInt64(11), work), "packed exact quotient over GF(101)")
    c.check(divided >= 10 and remainders >= 8, "packed division takes word-sized divisors")


def few_terms(spec: List[Tuple[Int, Int, Int, Int]], order: MonomialOrder, mut work: Budget) raises -> SparsePoly[ZZ]:
    """A polynomial from (coefficient, e0, e1, e2) terms."""
    var terms = List[Term[ZZ]]()
    for t in spec:
        var pairs = List[Tuple[Int64, UInt32]]()
        if t[1]:
            pairs.append((Int64(0), UInt32(t[1])))
        if t[2]:
            pairs.append((Int64(1), UInt32(t[2])))
        if t[3]:
            pairs.append((Int64(2), UInt32(t[3])))
        terms.append(Term[ZZ](Monomial(pairs^, work), ZZ(Integer(t[0]))))
    return SparsePoly[ZZ].from_terms(terms^, order, NoContext(), work)


def powers(mut c: Checker) raises:
    """Powers of few-term bases by the power recurrence agree with repeated heap products."""
    var work = Budget()
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1), Int64(2)], work)
    # Fateman-like; a negative non-unit leading coefficient; a leading term that
    # ties in the first variable, so the recurrence grades by the second; and
    # sums of about 105 bits against the 124-bit bound.
    var bases = [
        few_terms([(1, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1)], order, work),
        few_terms([(-1, 0, 0, 0), (3, 0, 1, 0), (-2, 1, 0, 0), (5, 0, 2, 3)], order, work),
        few_terms([(1, 1, 0, 0), (-4, 1, 1, 0), (7, 0, 0, 2), (1, 0, 0, 0)], order, work),
        few_terms([(3, 0, 0, 0), (-5, 0, 0, 1), (7, 0, 1, 0), (9, 1, 0, 0)], order, work),
    ]
    var exponents = [[3, 7, 16, 30], [3, 5, 12, 2], [4, 9, 13, 2], [3, 8, 20, 6]]
    var taken = 0
    for b in range(len(bases)):
        ref f = bases[b]
        for k in exponents[b]:
            var expected = f.copy()
            for _ in range(k - 1):
                expected = _heap_product(expected, f, work)
            var direct = _power_recurrence(f, k, work)
            c.check(Bool(direct), String("power recurrence applies, base ", b, ", exponent ", k))
            if direct:
                direct.value().check_invariants(work)
                c.check(direct.value() == expected, String("power recurrence, base ", b, ", exponent ", k))
            if _packed_power(f, k, work):
                taken += 1
            var got = f.pow(UInt64(k), work)
            c.check(got == expected, String("power, base ", b, ", exponent ", k))
    c.check(taken >= 3, "large powers of few-term bases take the recurrence")
    c.check(not _packed_power(bases[0], 2, work), "squares use the product kernel")
    var spread = few_terms([(1, 0, 0, 0), (1, 1, 0, 0), (2, 0, 2, 0), (5, 0, 0, 5)], order, work)
    c.check(not _packed_power(spread, 12, work), "wide exponent spreads keep repeated squaring")
    c.check(Bool(_power_recurrence(spread, 12, work)), "the recurrence itself still applies")


def integer_words(mut c: Checker) raises:
    """192-bit two's complement words convert exactly at every width boundary."""
    var two64 = Integer(1) << 64
    var values = List[Integer]()
    for base in [Integer(0), Integer(1) << 63, two64, Integer(1) << 127, Integer(1) << 128, Integer(1) << 190]:
        for delta in [-1, 0, 1]:
            values.append(base + delta)
            values.append(-(base + delta))
    values.append(Integer(12345678901234567))
    for v in values:
        # Two's complement words of v modulo 2^192.
        var m = v % (Integer(1) << 192)
        var w0 = UInt64(Int(m % two64 >> 32) << 32) | UInt64(Int(m % (Integer(1) << 32)))
        var w1 = UInt64(Int((m >> 64) % two64 >> 32) << 32) | UInt64(Int((m >> 64) % (Integer(1) << 32)))
        var w2 = UInt64(Int((m >> 128) % two64 >> 32) << 32) | UInt64(Int((m >> 128) % (Integer(1) << 32)))
        if v >= Integer(1) << 191 or v < -(Integer(1) << 191):
            continue
        c.check(_integer((w0, w1, w2)) == v, String("words convert ", v))


def narrow_bound(mut c: Checker) raises:
    """Products whose coefficient bound just fits 64-bit slots (and just misses) match the heap product."""
    var work = Budget()
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1), Int64(2)], work)
    # max|f| * |g|_1 < 2^(29 + 31) = 2^60, plus one for a square: the 64-bit slots.
    # 2^(31 + 33) = 2^64: the 128-bit slots.
    for bits in [29, 31]:
        var m = (1 << bits) - 1
        var f = few_terms([(m, 0, 0, 0), (-m, 1, 0, 0), (m, 0, 1, 0), (-m, 1, 1, 1)], order, work)
        var g = few_terms([(m, 0, 0, 0), (m, 1, 0, 0), (-m, 0, 1, 0), (m, 0, 0, 1)], order, work)
        c.check(f.mul(g, work) == _heap_product(f, g, work), String("product at the 64-bit bound, ", bits, " bits"))
        c.check(f.mul(f, work) == _heap_product(f, f, work), String("square at the 64-bit bound, ", bits, " bits"))
        var pf = _PackedPoly[ZZ].from_sparse(f, work)
        var pg = _PackedPoly[ZZ].from_sparse(g, work)
        var product = _packed_poly_product(pf.value(), pg.value(), work)
        c.check(Bool(product) and product.value().to_sparse(work) == _heap_product(f, g, work), "packed product at the 64-bit bound")


def main() raises:
    var c = Checker("sparse products")
    integer_words(c)
    narrow_bound(c)
    powers(c)
    var work = Budget()
    products(ZZ(Integer(1) << 127), c)
    products(QQ(Rational(1, 1009)), c)
    products(ModInt(1, Modulus(2, work)), c)
    products(ModInt(1, Modulus(101, work)), c)
    packed[ZZ](NoContext(), c)
    packed[ModInt](Modulus(101, work), c)
    packed[ModInt](Modulus((UInt64(1) << 62) - 57, work), c)
    packed_poly[ZZ](NoContext(), c)
    packed_poly[ModInt](Modulus(101, work), c)
    packed_powers(c)
    c.finish()
