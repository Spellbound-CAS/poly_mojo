"""Borrowed proof checks, prepared state and fused construction regressions."""
from apn_mojo import Integer, Rational
from polynomials import *
from polynomials.algebraic_types import _check_real
from polynomials.algebraic_specialization import _prepare_point
from polynomials._sparse_ops import _subtract_multiples
from polynomials._dense_reduction import _dense_reducer
from polynomials._subsets import _Recombination
from polynomials._lattice import _orthogonal
from support.corpus import parse_polynomial, as_mod, as_zz
from support.check import Checker


def witnesses(mut c: Checker, mut work: Budget) raises:
    var f = DensePoly[QQ].from_sparse(parse_polynomial("x0^2 - 2", MonomialOrder(), work), work)
    var value = real_root(f, QQ(Rational(1)), QQ(Rational(2)), work)
    c.check(_check_real(value, work), "borrowed Sturm evidence accepted")
    var before = value.copy()
    for kind in range(5):
        var bad = value.copy()
        if kind == 0:
            bad.sequence.quotients[0].coefficients[0] = QQ(Rational(3))
        elif kind == 1:
            bad.low = QQ(Rational(2))
        elif kind == 2:
            bad.high = QQ(Rational(3))
            bad.low = QQ(Rational(2))
        elif kind == 3:
            _ = bad.sequence.polynomials.pop()
        else:
            bad.sequence.polynomials[0].coefficients[0] = QQ(Rational(-3))
        c.check(not _check_real(bad, work), "borrowed checker rejects malformed witness")
    var tiny = Budget(0)
    try:
        _ = _check_real(value, tiny)
        c.check(False, "borrowed check ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and value.polynomial == before.polynomial and value.sequence.polynomials == before.sequence.polynomials and value.sequence.quotients == before.sequence.quotients, "budget failure preserves witness")
    var point: List[Tuple[Int64, RealAlgebraic]] = [(Int64(9), value.copy()), (Int64(2), value.copy()), (Int64(5), real_rational(QQ(Rational(3)), work))]
    var prepared = _prepare_point(point, work)
    c.check(prepared.index(2, work) == 1 and prepared.index(9, work) == 0 and prepared.index(5, work) == 2 and prepared.index(7, work) == -1, "sorted index retains original coordinate bindings")
    var p = parse_polynomial("x2*x10^2 - x9*x10^2 + x5*x10 + 1", MonomialOrder(), work)
    var proof = algebraic_specialize(p, point, work)
    c.check(check_algebraic_specialization(p, point, proof, work), "prepared point preserves dependent and zero coefficients")
    var swapped: List[Tuple[Int64, RealAlgebraic]] = [point[2].copy(), point[1].copy(), point[0].copy()]
    var reordered = algebraic_specialize(p, swapped, work)
    c.check(check_algebraic_specialization(p, point, reordered, work), "assignment permutation preserves transcript schedule")
    point.append((Int64(2), value.copy()))
    try:
        _ = _prepare_point(point, work)
        c.check(False, "duplicate prepared coordinate accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "duplicate prepared coordinate rejected")


def sparse[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], expected: SparsePoly[C], scaled_expected: SparsePoly[C], mut c: Checker, mut work: Budget) raises:
    var one = C.one(f.ring)
    var left = Term[C](Monomial([(Int64(0), UInt32(1))], work), one.copy())
    var right = Term[C](Monomial([(Int64(1), UInt32(1))], work), one.copy())
    var got = _subtract_multiples(f, left, g, right, work)
    got.check_invariants(work)
    c.check(got == expected, "shifted streams cancel their common leading term")
    var two = one.add(one, work)
    var three = two.add(one, work)
    var scaled_left = Term[C](left.monomial.copy(), three^)
    var scaled_right = Term[C](right.monomial.copy(), two^)
    c.check(_subtract_multiples(f, scaled_left, g, scaled_right, work) == scaled_expected, "nonunit multipliers retain exact coefficient arithmetic")
    c.check(_subtract_multiples(f, left, f, left, work).is_zero(), "identical translated streams cancel completely")
    var zero = SparsePoly[C].zero(f.order, f.ring, work)
    c.check(_subtract_multiples(zero, left, zero, right, work).is_zero(), "empty shifted streams")
    var zero_term = Term[C](left.monomial.copy(), C.zero(f.ring))
    c.check(_subtract_multiples(f, zero_term, g, zero_term, work).is_zero(), "zero multipliers skip streams")
    var before = f.copy()
    var tiny = Budget(0)
    try:
        _ = _subtract_multiples(f, left, g, right, tiny)
        c.check(False, "fused merge ignored work budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and before == f, "fused merge preserves borrowed input on exhaustion")


def sparse_cases(mut c: Checker, mut work: Budget) raises:
    var orders = [MonomialOrder(), MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1)], work), MonomialOrder(MonomialOrder.GRLEX, [Int64(0), Int64(1)], work), MonomialOrder(MonomialOrder.GREVLEX, [Int64(1), Int64(0)], work)]
    for order in orders:
        var f = parse_polynomial("x1^2 + 2*x0 + 1", order, work)
        var g = parse_polynomial("x0*x1 + x0 - 1", order, work)
        var expected = parse_polynomial("2*x0^2 - x0*x1 + x0 + x1", order, work)
        var scaled_expected = parse_polynomial("x0*x1^2 + 6*x0^2 - 2*x0*x1 + 3*x0 + 2*x1", order, work)
        sparse(f, g, expected, scaled_expected, c, work)
        sparse(as_zz(f, work), as_zz(g, work), as_zz(expected, work), as_zz(scaled_expected, work), c, work)
        var ctx = Modulus(2, work)
        sparse(as_mod(f, ctx, work), as_mod(g, ctx, work), as_mod(expected, ctx, work), as_mod(scaled_expected, ctx, work), c, work)
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0)], work)
    var huge = SparsePoly[QQ].from_terms([Term[QQ](Monomial([(Int64(0), UInt32.MAX)], work), QQ.one(NoContext()))], order, NoContext(), work)
    var shift = Term[QQ](Monomial([(Int64(0), UInt32(1))], work), QQ.one(NoContext()))
    try:
        _ = _subtract_multiples(huge, shift, huge, shift, work)
        c.check(False, "overflowing shifted streams accepted")
    except e:
        c.check(e.kind == PolyErrorKind.CAPACITY, "fused merge checks exponent capacity before cancellation")


def reduction[C: CoefficientField](a: DensePoly[C], b: DensePoly[C], g: DensePoly[C], expected: DensePoly[C], mut c: Checker, mut work: Budget) raises:
    var reducer = _dense_reducer(g, work)
    var before = a.copy()
    var before_b = b.copy()
    for limit in [UInt64(0), UInt64(4), UInt64(15), UInt64(35)]:
        var tiny = Budget(limit)
        try:
            var got = reducer.product(a, b, tiny)
            c.check(got == expected, "limited prepared reduction returns complete result")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "prepared reduction preserves typed exhaustion")
        c.check(a == before and b == before_b, "prepared reduction preserves both borrowed coefficient vectors")
        c.check(reducer.product(a, b, work) == expected, "prepared modulus reusable after interrupted reduction")
    var zero = DensePoly[C](a.variable, a.ring.copy())
    c.check(reducer.product(a, zero, work).is_zero(), "zero product needs no reduction")
    c.check(reducer.product(a, b, work) == expected, "prepared modulus reusable after zero product")
    var unit = DensePoly[C].from_coefficients(a.variable, a.ring, [C.one(a.ring)], work)
    var constant = _dense_reducer(unit, work)
    c.check(constant.product(a, b, work).is_zero(), "constant modulus has zero remainder")
    var large_values = List[C](length=65, fill=C.zero(a.ring))
    large_values[0] = C.one(a.ring)
    large_values[64] = C.one(a.ring)
    var large = DensePoly[C].from_coefficients(a.variable, a.ring, large_values^, work)
    var large_modulus = large.copy()
    large_modulus.coefficients[0] = C.one(a.ring).neg(work)
    var prepared_large = _dense_reducer(large_modulus, work)
    var two = C.one(a.ring).add(C.one(a.ring), work)
    var four = DensePoly[C].from_coefficients(a.variable, a.ring, [two.mul(two, work)], work)
    c.check(prepared_large.product(large, large, work) == four, "Karatsuba product reduced in place: (x^64+1)^2 modulo x^64-1 is four")


def reductions(mut c: Checker, mut work: Budget) raises:
    var order = MonomialOrder()
    var a = parse_polynomial("x0^2 + 2*x0 + 3", order, work)
    var b = parse_polynomial("x0^2 - x0 + 1", order, work)
    var modulus = parse_polynomial("2*x0^2 - 2", order, work)
    var expected = parse_polynomial("6", order, work)
    reduction(DensePoly[QQ].from_sparse(a, work), DensePoly[QQ].from_sparse(b, work), DensePoly[QQ].from_sparse(modulus, work), DensePoly[QQ].from_sparse(expected, work), c, work)
    var ctx = Modulus(7, work)
    reduction(DensePoly[ModInt].from_sparse(as_mod(a, ctx, work), work), DensePoly[ModInt].from_sparse(as_mod(b, ctx, work), work), DensePoly[ModInt].from_sparse(as_mod(modulus, ctx, work), work), DensePoly[ModInt].from_sparse(as_mod(expected, ctx, work), work), c, work)
    var g = DensePoly[QQ].from_sparse(modulus, work)
    var reducer = _dense_reducer(g, work)
    g.coefficients[0] = QQ(Rational(100))
    c.check(reducer.product(DensePoly[QQ].from_sparse(a, work), DensePoly[QQ].from_sparse(b, work), work) == DensePoly[QQ].from_sparse(expected, work), "prepared modulus owns its coefficients")
    try:
        _ = _dense_reducer(DensePoly[QQ](0, NoContext()), work)
        c.check(False, "zero prepared modulus accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero modulus rejected")
    var wrong = DensePoly[QQ].from_sparse(a, work)
    wrong.variable = 1
    try:
        _ = reducer.product(wrong, wrong, work)
        c.check(False, "prepared reduction accepts different variable")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "prepared variable context enforced")
    var modular = _dense_reducer(DensePoly[ModInt].from_sparse(as_mod(modulus, ctx, work), work), work)
    var other = DensePoly[ModInt].from_sparse(as_mod(a, Modulus(11, work), work), work)
    try:
        _ = modular.product(other, other, work)
        c.check(False, "prepared reduction accepts different prime")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "prepared prime context enforced")


def selections(mut c: Checker, mut work: Budget) raises:
    var search = _Recombination(5, work)
    var count = 0
    while search.ready():
        count += 1
        search.advance(work)
    c.check(count == 15, "five singleton and ten pair candidates")
    var pairs = _Recombination(4, work)
    for _ in range(4):
        pairs.advance(work)
    var expected: List[List[Int]] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
    for selection in expected:
        c.check(pairs.selection.indices == selection, "pair traversal remains lexicographic")
        pairs.advance(work)
    c.check(not pairs.ready(), "pair exhaustion terminates traversal")
    var remaining = _Recombination(5, work)
    remaining.advance(work)
    remaining.accept(work)
    c.check(remaining.active == [0, 2, 3, 4] and remaining.selection.indices == [0], "accept removes selected factor and restarts at singleton")
    for limit in [UInt64(0), UInt64(5), UInt64(8)]:
        var tiny = Budget(limit)
        try:
            remaining.accept(tiny)
            c.check(False, "selection removal ignored budget")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT and remaining.active == [0, 2, 3, 4] and remaining.selection.indices == [0], "failed removal preserves paired selection state")
    var pair_removal = _Recombination(4, work)
    for _ in range(5):
        pair_removal.advance(work)
    pair_removal.accept(work)
    c.check(pair_removal.active == [1, 3] and pair_removal.selection.indices == [0], "noncontiguous pair removal preserves active factor order")
    for size in range(2):
        var empty = _Recombination(size, work)
        c.check(not empty.ready(), "no recombination candidate for fewer than two factors")
        try:
            empty.accept(work)
            c.check(False, "accepted empty selection")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "empty selection acceptance rejected")


def cyclic_reduction(mut c: Checker, mut work: Budget) raises:
    var ctx = Modulus(65521, work)
    for size in [8, 64]:
        var a = List[ModInt]()
        var b = List[ModInt]()
        var expected = List[ModInt]()
        for i in range(size):
            a.append(ModInt(UInt64(3 * i + 2), ctx))
            b.append(ModInt(UInt64(i * i + 1), ctx))
        # Independent cyclic convolution, using native integer residues.
        for k in range(size):
            var value = UInt64(0)
            for i in range(size):
                var j = (k - i + size) % size
                value = (value + UInt64(3 * i + 2) * UInt64(j * j + 1)) % ctx.value
            expected.append(ModInt(value, ctx))
        var coefficients = List[ModInt](length=size + 1, fill=ModInt.zero(ctx))
        coefficients[0] = ModInt(2, ctx).neg(work)
        coefficients[size] = ModInt(2, ctx)
        var modulus = DensePoly[ModInt].from_coefficients(0, ctx, coefficients^, work)
        var prepared = _dense_reducer(modulus, work)
        var actual = prepared.product(DensePoly[ModInt].from_coefficients(0, ctx, a^, work), DensePoly[ModInt].from_coefficients(0, ctx, b^, work), work)
        c.check(actual == DensePoly[ModInt].from_coefficients(0, ctx, expected^, work), "dense modular product matches independent cyclic convolution")


def lattice(mut c: Checker, mut work: Budget) raises:
    var input: List[List[ZZ]] = [[ZZ(Integer(1)), ZZ(Integer(1)), ZZ(Integer(1))], [ZZ(Integer(-1)), ZZ(Integer(0)), ZZ(Integer(2))], [ZZ(Integer(3)), ZZ(Integer(5)), ZZ(Integer(6))]]
    var result = _orthogonal(input, work)
    c.check(result.norms == [QQ(Rational(3)), QQ(Rational(14, 3)), QQ(Rational(9, 14))], "hand-derived orthogonal squared norms")
    c.check(result.mu[1][0] == QQ(Rational(1, 3)) and result.mu[2][0] == QQ(Rational(14, 3)) and result.mu[2][1] == QQ(Rational(13, 14)), "hand-derived rational projection coefficients")


def main() raises:
    var work = Budget()
    var c = Checker("reused operation state")
    witnesses(c, work)
    sparse_cases(c, work)
    reductions(c, work)
    cyclic_reduction(c, work)
    selections(c, work)
    lattice(c, work)
    c.finish()
