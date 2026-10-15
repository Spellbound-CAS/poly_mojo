"""Closed-form sparse products across rings, orders, collisions and limits."""
from apn_mojo import Integer, Rational
from polynomials import Budget, CoefficientRing, ZZ, QQ, ModInt, Modulus, NoContext, SparsePoly, Term, Monomial, MonomialOrder, PolyErrorKind
from support.products import geometric, expected_product
from support.check import Checker
from support.gen import polynomial
from support.rng import SplitMix64
from polynomials.sparse import _heap_product


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


def main() raises:
    var c = Checker("sparse products")
    var work = Budget()
    products(ZZ(Integer(1) << 127), c)
    products(QQ(Rational(1, 1009)), c)
    products(ModInt(1, Modulus(2, work)), c)
    products(ModInt(1, Modulus(101, work)), c)
    c.finish()
