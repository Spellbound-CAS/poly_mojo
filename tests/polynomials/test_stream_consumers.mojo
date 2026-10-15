"""Independent identities for consumed tree nodes and ordered term streams."""
from apn_mojo import Rational, Integer, ExactComplex
from polynomials import *
from support.check import Checker
from support.products import geometric, expected_product
from polynomials.sparse import _ProductCursor, _ProductOrder


def product_checks[C: CoefficientRing](base: C, mut c: Checker, mut work: Budget) raises:
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(0)], work)
        var a = geometric(34, 0, False, base, order, work)
        var b = geometric(34, 0, True, base, order, work)
        var expected = expected_product(34, False, True, base, order, work)
        var generators: List[SparsePoly[C]] = [a.copy()]
        var certificate = MembershipCertificate[C]([b.copy()])
        c.check(check_membership(expected, generators, certificate, work), "independent alternating convolution")
        var zero = SparsePoly[C].zero(order, base.context(), work)
        var cancel = MembershipCertificate[C]([b.copy(), b.neg(work)])
        var repeated: List[SparsePoly[C]] = [a.copy(), a.copy()]
        c.check(check_membership(zero, repeated, cancel, work), "complete cancellation across products")
        for units in [0, 100, 1000]:
            var limited = Budget(UInt64(units))
            try:
                c.check(check_membership(expected, generators, certificate, limited), "limited checker returns only complete identity")
            except e:
                c.check(e.kind == PolyErrorKind.LIMIT, "checker fails closed during consumption")
        repeated[1].terms[0].coefficient = C.zero(base.context())
        try:
            _ = check_membership(a, repeated, cancel, work)
            c.check(False, "early mismatch hid malformed later generator")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "all streamed operands validated")


def main() raises:
    var c = Checker("stream consumers")
    var work = Budget()
    var prime = Modulus(65537, work)
    # Base-three digits are lexicographic exponents on (-7, 0, 9).
    # Numeric index order is therefore an independent lex comparator oracle.
    var cursors = List[_ProductCursor[ModInt]]()
    for value in range(27):
        var pairs = List[Tuple[Int64, UInt32]]()
        for identifier, place in [(Int64(-7), 9), (Int64(0), 3), (Int64(9), 1)]:
            var exponent = value // place % 3
            if exponent:
                pairs.append((identifier, UInt32(exponent)))
        var cursor = _ProductCursor[ModInt](value, prime)
        cursor.term = Term[ModInt](Monomial(pairs^, work), ModInt(1, prime))
        cursors.append(cursor^)
    var sorted_order = MonomialOrder(MonomialOrder.LEX, [Int64(-7), Int64(0), Int64(9)], work)
    for a in range(27):
        for b in range(27):
            c.check(_ProductOrder[ModInt, True].compare(cursors, a, b, sorted_order, work) == Int(a > b) - Int(a < b), "internal ascending-variable lex order")
    product_checks(ZZ(Integer(1) << 127), c, work)
    product_checks(QQ(Rational(1, 1009)), c, work)
    product_checks(QQi(ExactComplex(Rational(2, 7), Rational(3, 11))), c, work)
    product_checks(ModInt(1, Modulus(2, work)), c, work)
    for n in [1, 7, 32, 65]:
        var points = List[ModInt]()
        var coefficients = List[ModInt]()
        for i in range(n):
            points.append(ModInt(UInt64(i), prime))
            coefficients.append(ModInt(UInt64(i + 1), prime))
        var f = DensePoly[ModInt].from_coefficients(-17, prime, coefficients^, work)
        var plan = MultipointPlan[ModInt].prepare_interpolation(-17, prime, points, work)
        var values = List[ModInt]()
        for point in points:
            values.append(evaluate_dense(f, point, work))
        for _ in range(2):
            c.check(plan.evaluate(f, work) == values, "tree agrees with Horner including odd frontiers")
            c.check(plan.interpolate(values, work) == f, "unique degree-bounded interpolant and reusable plan")
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(0), Int64(1)], work)
        var terms = List[Term[QQ]]()
        for i in range(65):
            terms.append(Term[QQ](Monomial([(Int64(0), UInt32(i + 1)), (Int64(1), UInt32(65 - i))], work), QQ(Rational(1))))
        var p = SparsePoly[QQ].from_terms(terms^, order, NoContext(), work)
        var one = SparsePoly[QQ].from_terms([Term[QQ](Monomial(), QQ(Rational(1)))], order, NoContext(), work)
        # Two independently planted identities: (p+1)(p-1)=p^2-1 and
        # p*(p+1) + (-p)*p = p, with cancellation across product streams.
        var a = p.add(one, work)
        var b = p.sub(one, work)
        c.check(a.mul(b, work) == p.mul(p, work).sub(one, work), "ordered product distributivity")
        var generators: List[SparsePoly[QQ]] = [a.copy(), p.copy()]
        var cofactors: List[SparsePoly[QQ]] = [p.copy(), p.neg(work)]
        var cert = MembershipCertificate[QQ](cofactors^)
        c.check(check_membership(p, generators, cert, work), "cancelling membership products")
        c.check(not check_membership(a, generators, cert, work), "corrupted membership target rejected")
        c.check(deserialize_sparse[QQ](serialize_sparse(p, work), NoContext(), work) == p, "streamed sparse interchange")
        if kind:
            var dividend = a.mul(b, work).add(one, work)
            var divisors: List[SparsePoly[QQ]] = [b.copy(), a.copy()]
            var division = divide_by_list(dividend, divisors, work)
            c.check(division.quotients[0] == a and division.quotients[1].is_zero() and division.remainder == one, "first eligible divisor and planted remainder")
            c.check(check_list_division(dividend, divisors, division, work), "division identity consumer")
        var snapshot = serialize_sparse(p, work)
        var tiny = Budget(100)
        try:
            _ = p.mul(p, tiny)
            c.check(False, "limited multiplication accepted")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "stream budget exhaustion")
        c.check(snapshot == serialize_sparse(p, work), "borrowed input preserved")
    var hand = "POLYNOMIALS 1\nDENSE\nQQ\n-17\n3\n1\n-2\n3\n"
    var dense = deserialize_dense[QQ](hand, NoContext(), work)
    c.check(serialize_dense(dense, work) == hand, "canonical bytes match handwritten fixture")
    # (1+...+x^128)^2 has triangular coefficients. Its long reduction
    # exercises bounded-frontier compaction and slot recycling.
    var triangle = List[ModInt]()
    var ones = List[ModInt]()
    for i in range(257):
        triangle.append(ModInt(UInt64(min(i + 1, 257 - i)), prime))
    for _ in range(129):
        ones.append(ModInt(1, prime))
    var order = MonomialOrder(MonomialOrder.GREVLEX, [Int64(0)], work)
    var square = DensePoly[ModInt].from_coefficients(0, prime, triangle^, work).to_sparse(order, work)
    var root = DensePoly[ModInt].from_coefficients(0, prime, ones^, work).to_sparse(order, work)
    var divisors: List[SparsePoly[ModInt]] = [root.copy()]
    var reduced = divide_by_list(square, divisors, work)
    c.check(reduced.quotients[0] == root and reduced.remainder.is_zero(), "bounded frontier triangular convolution")
    c.check(check_list_division(square, divisors, reduced, work), "large modular division checker")
    var grid = List[Term[ModInt]]()
    var grid_square = List[Term[ModInt]]()
    for i in range(33):
        for j in range(33):
            var pairs = List[Tuple[Int64, UInt32]]()
            if i:
                pairs.append((Int64(0), UInt32(i)))
            if j:
                pairs.append((Int64(1), UInt32(j)))
            var monomial = Monomial(pairs^, work)
            grid_square.append(Term[ModInt](monomial.copy(), ModInt(UInt64(min(i + 1, 33 - i) * min(j + 1, 33 - j)), prime)))
            if i < 17 and j < 17:
                grid.append(Term[ModInt](monomial^, ModInt(1, prime)))
    var grid_poly = SparsePoly[ModInt].from_terms(grid^, MonomialOrder(), prime, work)
    var grid_expected = SparsePoly[ModInt].from_terms(grid_square^, MonomialOrder(), prime, work)
    c.check(grid_poly.mul(grid_poly, work) == grid_expected, "mixed-support KERNEL product has tensor triangular coefficients")
    for units in [0, 100, 10000, 100000]:
        var limited = Budget(UInt64(units))
        try:
            _ = divide_by_list(square, divisors, limited)
            c.check(False, "reduction accepted insufficient budget")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "partial reduction is not returned")
        c.check(check_list_division(square, divisors, reduced, work), "failed reduction preserves borrowed data")
    for bad in ["POLYNOMIALS 1\nSPARSE\nQQ\n0\n0\n2\n1\n0\n1\n1\n1\n0\n1\n2\n", "POLYNOMIALS 1\nSPARSE\nQQ\n0\n0\n1\n1\n0\n0\n1\n", "POLYNOMIALS 1\nSPARSE\nQQ\n0\n0\n1\n0\n0\n"]:
        try:
            _ = deserialize_sparse[QQ](bad, NoContext(), work)
            c.check(False, "malformed streamed term accepted")
        except e:
            c.check(e.kind == PolyErrorKind.ENCODING, "malformed streamed term rejected")
    c.finish()
