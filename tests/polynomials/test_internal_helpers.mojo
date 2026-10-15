"""Hand-derived coefficient images, subset traversal and sparse construction."""
from apn_mojo import Integer, Rational
from polynomials import Budget, PolyErrorKind, QQ, ZZ, ModInt, Modulus, NoContext, DensePoly, SparsePoly, Monomial, MonomialOrder, Term, RealAlgebraic, real_rational
from polynomials.factor import _integer, _rational
from polynomials._images import _integer_image, _rational_image, _mod_image, _sparse_mod_image
from polynomials._subsets import _Subsets
from polynomials._views import _join
from polynomials.algebraic_specialization import _groups
from support.corpus import parse_polynomial
from support.check import Checker


def conversions(mut c: Checker) raises:
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(17, NoContext(), [QQ(Rational(1, 6)), QQ.zero(NoContext()), QQ(Rational(-3, 10)), QQ(Rational(7, 15))], work)
    var expected = DensePoly[ZZ].from_coefficients(17, NoContext(), [ZZ(Integer(5)), ZZ(Integer(0)), ZZ(Integer(-9)), ZZ(Integer(14))], work)
    var embedded = DensePoly[QQ].from_coefficients(17, NoContext(), [QQ(Rational(5)), QQ(Rational(0)), QQ(Rational(-9)), QQ(Rational(14))], work)
    c.check(_integer(f, work) == expected, "LCM 30 clears denominators without primitive normalization")
    c.check(_rational(expected, work) == embedded, "integer embedding preserves variable and internal zeros")
    var ctx = Modulus(7, work)
    var modular = DensePoly[ModInt].from_coefficients(17, ctx, [ModInt(5, ctx), ModInt(0, ctx), ModInt(5, ctx)], work)
    c.check(_mod_image(expected, ctx, work) == modular, "modular image removes vanishing leading coefficient")
    var zero = DensePoly[QQ](17, NoContext())
    c.check(_integer(zero, work).is_zero() and _rational(_integer(zero, work), work) == zero, "zero coefficient images")
    var constant = DensePoly[QQ].from_coefficients(17, NoContext(), [QQ(Rational(-4, 9))], work)
    c.check(_integer(constant, work).coefficients[0] == ZZ(Integer(-4)), "negative constant retains numerator content")
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(17)], work)
        c.check(_integer_image(f.to_sparse(order, work), work) == expected.to_sparse(order, work), "dense and sparse denominator clearing commute")
        c.check(_rational_image(expected.to_sparse(order, work), work) == embedded.to_sparse(order, work), "dense and sparse integer embedding commute")
        c.check(_sparse_mod_image(expected.to_sparse(order, work), ctx, work) == modular.to_sparse(order, work), "sparse modular image removes zero residues")
    var saved = f.copy()
    for limit in [UInt64(0), UInt64(10), UInt64(100), UInt64(1000)]:
        var tiny = Budget(limit)
        try:
            c.check(_integer(f, tiny) == expected, "limited conversion is complete")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "coefficient conversion fails closed")
        c.check(f == saved, "conversion failure preserves input")


def subsets(mut c: Checker) raises:
    var work = Budget()
    var counts: List[Int] = [1, 6, 15, 20, 15, 6, 1]
    for size in range(7):
        var cursor = _Subsets(6, size, work)
        var count = 1
        while cursor.advance(work):
            count += 1
        c.check(count == counts[size], "binomial number of subsets")
    var cursor = _Subsets(4, 2, work)
    var expected: List[List[Int]] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
    for i in range(len(expected)):
        c.check(cursor.indices == expected[i], "lexicographic combination order")
        c.check(cursor.advance(work) == (i + 1 < len(expected)), "exact subset exhaustion")
    var wide = _Subsets(70, 1, work)
    var count = 1
    while wide.advance(work):
        count += 1
    c.check(count == 70, "subsets have no machine-word bitmask limit")
    for size in [-1, 5]:
        try:
            _ = _Subsets(4, size, work)
            c.check(False, "invalid subset size accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "invalid subset size is DOMAIN")
    var limited = _Subsets(5, 3, work)
    var saved = limited.indices.copy()
    var tiny = Budget(0)
    try:
        _ = limited.advance(tiny)
        c.check(False, "unbudgeted subset advance accepted")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and limited.indices == saved, "subset failure preserves cursor")


def builders(mut c: Checker) raises:
    var work = Budget()
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(1), Int64(0)], work)
        var a = parse_polynomial("x0 + 1", order, work)
        var b = parse_polynomial("-x0 - 1", order, work)
        var expected = parse_polynomial("x0 + 1 - x0*x1^2 - x1^2", order, work)
        var coefficients: List[SparsePoly[QQ]] = [a.copy(), SparsePoly[QQ].zero(order, NoContext(), work), b.copy()]
        c.check(_join(coefficients, 1, a, work) == expected, "joining preserves gaps and signs")
        c.check(_join(List[SparsePoly[QQ]](), 1, a, work).is_zero(), "empty join")
        var overlap: List[SparsePoly[QQ]] = [parse_polynomial("x1", order, work), parse_polynomial("-1", order, work)]
        c.check(_join(overlap, 1, a, work).is_zero(), "joining combines and cancels overlapping terms")
        var point: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), real_rational(QQ(Rational(1)), work))]
        var grouped = _groups(expected, point, work)
        c.check(len(grouped[0]) == 2, "two remaining monomial groups")
        for i in range(len(grouped[0])):
            var degree = grouped[0][i].exponent(1, work)
            c.check((degree == 0 or degree == 2) and grouped[1][i] == (a if degree == 0 else b), "group coefficients follow hand decomposition")
        var saved = expected.copy()
        for limit in [UInt64(0), UInt64(30), UInt64(300)]:
            var tiny = Budget(limit)
            try:
                _ = _groups(expected, point, tiny)
                c.check(True, "limited grouping completes")
            except e:
                c.check(e.kind == PolyErrorKind.LIMIT, "grouping fails closed")
            c.check(expected == saved, "grouping preserves input")
    var order = MonomialOrder()
    var interleaved = parse_polynomial("x0*x1^2 + x1 + 1", order, work)
    var point: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), real_rational(QQ(Rational(1)), work))]
    var grouped = _groups(interleaved, point, work)
    c.check(grouped[0][0].exponent(1, work) == 0 and grouped[0][1].exponent(1, work) == 2 and grouped[0][2].exponent(1, work) == 1, "grouping preserves first-occurrence certificate order")
    point.append((Int64(1), real_rational(QQ(Rational(2)), work)))
    grouped = _groups(interleaved, point, work)
    c.check(len(grouped[0]) == 1 and len(grouped[0][0].pairs) == 0 and grouped[1][0] == interleaved, "complete assignment forms one coefficient group")
    var zero = SparsePoly[QQ].zero(order, NoContext(), work)
    c.check(len(_groups(zero, point, work)[0]) == 0, "zero has no specialization groups")
    var overflowing = SparsePoly[QQ].from_terms([Term[QQ](Monomial([(Int64(1), UInt32.MAX)], work), QQ.one(NoContext()))], order, NoContext(), work)
    var values: List[SparsePoly[QQ]] = [zero.copy(), overflowing.copy()]
    try:
        _ = _join(values, 1, zero, work)
        c.check(False, "join exponent overflow accepted")
    except e:
        c.check(e.kind == PolyErrorKind.CAPACITY and values[1] == overflowing, "join overflow preserves borrowed coefficients")


def main() raises:
    var c = Checker("internal helpers")
    conversions(c)
    subsets(c)
    builders(c)
    c.finish()
