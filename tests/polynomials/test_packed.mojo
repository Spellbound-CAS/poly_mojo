"""PackedPoly through the public surface: round trips, products, powers and exact quotients against SparsePoly."""
from apn_mojo import Integer, Rational
from polynomials import Budget, CoefficientRing, ZZ, QQ, ModInt, Modulus, NoContext, SparsePoly, Term, Monomial, MonomialOrder, PolyError, PolyErrorKind, divide_exact
from polynomials import PackedPoly, packed_mul, packed_pow, packed_divide_exact
from support.check import Checker
from support.gen import polynomial
from support.rng import SplitMix64


def pack[C: CoefficientRing](p: SparsePoly[C], mut work: Budget) raises -> PackedPoly[C]:
    var packed = PackedPoly[C].from_sparse(p, work)
    if not packed:
        raise Error("expected a packable polynomial")
    return packed.take()


def sum_of_variables[C: CoefficientRing](ctx: C.Context, order: MonomialOrder, n: Int, mut work: Budget) raises -> SparsePoly[C]:
    """1 + x0 + ... + x(n-1)."""
    var terms = List[Term[C]]()
    terms.append(Term[C](Monomial(), C.one(ctx)))
    for i in range(n):
        terms.append(Term[C](Monomial([(Int64(i), UInt32(1))], work), C.one(ctx)))
    return SparsePoly[C].from_terms(terms^, order, ctx, work)


def agree[C: CoefficientRing](ctx: C.Context, name: String, mut c: Checker) raises:
    """Random operands: every packed operation equals its SparsePoly counterpart."""
    var work = Budget()
    var rng = SplitMix64(4242)
    for n in [1, 2, 3, 4]:
        var variables = List[Int64]()
        for i in range(n):
            variables.append(Int64(i))
        var order = MonomialOrder(MonomialOrder.LEX, variables^, work)
        for round in range(6):
            var bits = 90 if round == 5 else 0
            var a = polynomial[C](ctx, order, rng, work, terms=12, variables=n, degree=5, coefficient_bound=50, coefficient_bits=bits)
            var b = polynomial[C](ctx, order, rng, work, terms=9, variables=n, degree=4, coefficient_bound=50)
            var pa = pack(a, work)
            var pb = pack(b, work)
            c.check(pa.to_sparse(work) == a and len(pa) == len(a.terms), name + ": round trip")
            var product = packed_mul(pa, pb, work)
            var expected = a.mul(b, work)
            c.check(Bool(product) and product.value().to_sparse(work) == expected, name + ": product")
            if product:
                c.check(product.value() == pack(expected, work), name + ": product radix is exact")
                if not b.is_zero():
                    var quotient = packed_divide_exact(product.value(), pb, work)
                    c.check(Bool(quotient) and quotient.value().to_sparse(work) == a, name + ": exact quotient")
                    if quotient:
                        c.check(quotient.value() == pack(a, work), name + ": quotient radix is exact")
                    var bumped = expected.add(SparsePoly[C].constant(C.one(ctx), order, work), work)
                    var inexact = packed_divide_exact(pack(bumped, work), pb, work)
                    c.check(Bool(inexact) == Bool(divide_exact(bumped, b, work)), name + ": inexact quotient agrees with divide_exact")
            for k in [0, 1, 2, 3, 5]:
                var power = packed_pow(pb, UInt64(k), work)
                c.check(Bool(power) and power.value().to_sparse(work) == b.pow(UInt64(k), work), String(name, ": power ", k))


def powers_and_quotients(mut c: Checker) raises:
    """Fateman-style powers take the packed kernels; quotients of powers are exact."""
    var work = Budget()
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1), Int64(2), Int64(3)], work)
    var f = sum_of_variables[ZZ](NoContext(), order, 4, work)
    var pf = pack(f, work)
    var p10 = packed_pow(pf, 10, work)
    c.check(Bool(p10) and p10.value().to_sparse(work) == f.pow(10, work), "fateman power")
    var p4 = packed_pow(pf, 4, work)
    var q = packed_divide_exact(p10.value(), p4.value(), work)
    c.check(Bool(q) and q.value() == packed_pow(pf, 6, work).value(), "quotient of powers")
    var one = SparsePoly[ZZ].constant(ZZ(Integer(1)), order, work)
    var shifted = pack(f.pow(10, work).add(one, work), work)
    c.check(not packed_divide_exact(shifted, p4.value(), work), "a remainder means no quotient")
    var square = packed_mul(p10.value(), p10.value(), work)
    c.check(Bool(square) and square.value() == packed_pow(pf, 20, work).value(), "square of a power")
    var p = Modulus(101, work)
    var g = sum_of_variables[ModInt](p, order, 4, work)
    var pg = pack(g, work)
    var g12 = packed_pow(pg, 12, work)
    c.check(Bool(g12) and g12.value().to_sparse(work) == g.pow(12, work), "power over GF(101)")
    var g5 = packed_pow(pg, 5, work)
    var gq = packed_divide_exact(g12.value(), g5.value(), work)
    c.check(Bool(gq) and gq.value() == packed_pow(pg, 7, work).value(), "quotient over GF(101)")


def limits(mut c: Checker) raises:
    """Orders that do not pack, key ranges past 2**62, zero, and typed errors."""
    var work = Budget()
    var lex = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1)], work)
    var grevlex = MonomialOrder(MonomialOrder.GREVLEX, [Int64(0), Int64(1)], work)
    var other = sum_of_variables[ZZ](NoContext(), grevlex, 2, work)
    c.check(not PackedPoly[ZZ].from_sparse(other, work), "GREVLEX does not pack")
    var reversed = MonomialOrder(MonomialOrder.LEX, [Int64(1), Int64(0)], work)
    c.check(not PackedPoly[ZZ].from_sparse(sum_of_variables[ZZ](NoContext(), reversed, 2, work), work), "LEX must list variables 0..n-1 in order")
    var terms = List[Term[ZZ]]()
    terms.append(Term[ZZ](Monomial(), ZZ(Integer(1))))
    terms.append(Term[ZZ](Monomial([(Int64(0), UInt32(1 << 30)), (Int64(1), UInt32(1 << 30))], work), ZZ(Integer(1))))
    var wide = SparsePoly[ZZ].from_terms(terms^, lex, NoContext(), work)
    var pw = pack(wide, work)
    c.check(not packed_mul(pw, pw, work), "a product past 2**62 has no packed form")
    c.check(not packed_pow(pw, 2, work), "a power past 2**62 has no packed form")
    var p1 = packed_pow(pw, 1, work)
    c.check(Bool(p1) and p1.value() == pw, "the first power fits")
    var zero = pack(SparsePoly[ZZ].zero(lex, NoContext(), work), work)
    c.check(zero.is_zero() and len(zero) == 0, "zero packs")
    var x = pack(sum_of_variables[ZZ](NoContext(), lex, 2, work), work)
    var zx = packed_mul(zero, x, work)
    c.check(Bool(zx) and zx.value().is_zero(), "zero product")
    var zq = packed_divide_exact(zero, x, work)
    c.check(Bool(zq) and zq.value().is_zero(), "zero quotient")
    var z0 = packed_pow(zero, 0, work)
    c.check(Bool(z0) and z0.value().to_sparse(work) == SparsePoly[ZZ].constant(ZZ(Integer(1)), lex, work), "zero to the zero is one")
    try:
        _ = packed_divide_exact(x, zero, work)
        c.check(False, "division by zero raises")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "division by zero raises DOMAIN")
    var three = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1), Int64(2)], work)
    var y = pack(sum_of_variables[ZZ](NoContext(), three, 3, work), work)
    try:
        _ = packed_mul(x, y, work)
        c.check(False, "different orders raise")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "different orders raise DOMAIN")
    var constant_order = MonomialOrder(MonomialOrder.LEX, List[Int64](), work)
    var seven = pack(SparsePoly[ZZ].constant(ZZ(Integer(7)), constant_order, work), work)
    var s3 = packed_pow(seven, 3, work)
    c.check(Bool(s3) and s3.value().to_sparse(work) == SparsePoly[ZZ].constant(ZZ(Integer(343)), constant_order, work), "constants with no variables")
    var s1 = packed_divide_exact(s3.value(), seven, work)
    c.check(Bool(s1) and s1.value().to_sparse(work) == SparsePoly[ZZ].constant(ZZ(Integer(49)), constant_order, work), "constant quotient")


def main() raises:
    var c = Checker("packed")
    var work = Budget()
    agree[ZZ](NoContext(), "ZZ", c)
    agree[QQ](NoContext(), "QQ", c)
    agree[ModInt](Modulus(101, work), "GF(101)", c)
    agree[ModInt](Modulus((UInt64(1) << 62) - 57, work), "GF(2^62 - 57)", c)
    powers_and_quotients(c)
    limits(c)
    c.finish()
