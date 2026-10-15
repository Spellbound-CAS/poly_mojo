"""Division identities, independent corpus cases, and malformed certificates."""
from apn_mojo import Integer
from polynomials import Budget, PolyErrorKind, CoefficientRing, CoefficientField, ZZ, QQ, ModInt, Modulus, NoContext, MonomialOrder, DensePoly, SparsePoly, Term
from polynomials.division import divmod, pseudo_divmod, divide_exact, divide_by_list
from polynomials.certificates import check_division, check_list_division
from polynomials.coefficients import _pow
from support.check import Checker
from support.gen import coefficient, polynomial
from support.rng import SplitMix64
from support.corpus import read_cases, parse_polynomial, as_zz

def check_ring[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    var work = Budget()
    var rng = SplitMix64(772)
    var order = MonomialOrder(MonomialOrder.GREVLEX, [0, 1, 2], work)
    for i in range(35):
        var a = polynomial[C](ctx, order, rng, work, terms=3, variables=3, degree=2)
        var b = polynomial[C](ctx, order, rng, work, terms=3, variables=3, degree=2)
        if b.is_zero():
            continue
        var product = a.mul(b, work)
        var exact = divide_exact(product, b, work)
        c.check(Bool(exact) and exact.value() == a, "exact generated product")
        var stored = divide_exact(product.reorder(MonomialOrder(), work), b.reorder(MonomialOrder(), work), work)
        c.check(Bool(stored) and stored.value() == a.reorder(MonomialOrder(), work), "exact division in storage order")
        var da = List[C]()
        var db = List[C]()
        for _ in range(i % 9):
            da.append(coefficient[C](ctx, rng.below(11) - 5, work))
        for _ in range(i % 5):
            db.append(coefficient[C](ctx, rng.below(11) - 5, work))
        db.append(C.one(ctx).add(C.one(ctx), work))
        var f = DensePoly[C].from_coefficients(0, ctx, da^, work)
        var g = DensePoly[C].from_coefficients(0, ctx, db^, work)
        var result = pseudo_divmod(f, g, work)
        var scaled = f.scale(_pow(g.coefficients[len(g.coefficients) - 1], UInt64(max(f.degree() - g.degree() + 1, 0)), work), work)
        c.check(check_division(scaled, g, result, work), "pseudo division identity and degree")
        comptime if conforms_to(C, CoefficientField):
            var ordinary = divmod(f, g, work)
            c.check(check_division(f, g, ordinary, work), "field division identity and degree")
            ordinary.remainder = ordinary.remainder.add(DensePoly[C].from_coefficients(0, ctx, [C.one(ctx)], work), work)
            c.check(not check_division(f, g, ordinary, work), "corrupt division rejected")
            var divisors: List[SparsePoly[C]] = [a.copy(), b.copy()]
            if not a.is_zero():
                var listing = divide_by_list(product.add(a, work), divisors, work)
                c.check(check_list_division(product.add(a, work), divisors, listing, work), "list division identity")
                for term in listing.remainder.terms:
                    for divisor in divisors:
                        c.check(not divisor.leading(work).monomial.divides(term.monomial, work), "list remainder irreducible")
                listing.remainder = listing.remainder.add(SparsePoly[C].constant(C.one(ctx), order, work), work)
                c.check(not check_list_division(product.add(a, work), divisors, listing, work), "corrupt list division rejected")
    var zero = SparsePoly[C].zero(order, ctx, work)
    var one = SparsePoly[C].constant(C.one(ctx), order, work)
    var x = SparsePoly[C].variable(0, order, ctx, work)
    c.check(not divide_exact(one, x, work), "nonmultiple refused")
    c.check(divide_exact(zero, x, work).value().is_zero(), "zero quotient")
    try:
        _ = divide_exact(one, zero, work)
        c.check(False, "zero exact divisor accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero divisor DOMAIN")
    var tiny = Budget(0)
    try:
        _ = divide_exact(x, one, tiny)
        c.check(False, "exact division ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "division budget exhaustion")

def main() raises:
    var work = Budget()
    var c = Checker("division")
    check_ring[ZZ](NoContext(), c)
    check_ring[QQ](NoContext(), c)
    check_ring[ModInt](Modulus(101, work), c)
    var count = 0
    for record in read_cases():
        var op = record.field("op ")
        if op != "divmod" and op != "pseudo_divmod":
            continue
        var order = MonomialOrder()
        var ftext = record.field("in f ")
        var gtext = record.field("in g ")
        var qtext = record.field("out quotient ")
        var rtext = record.field("out remainder ")
        var f = parse_polynomial(ftext, order, work)
        var g = parse_polynomial(gtext, order, work)
        var q = parse_polynomial(qtext, order, work)
        var r = parse_polynomial(rtext, order, work)
        if op == "divmod":
            var got = divmod(DensePoly[QQ].from_sparse(f, work), DensePoly[QQ].from_sparse(g, work), work)
            c.check(got.quotient.to_sparse(order, work) == q and got.remainder.to_sparse(order, work) == r, record.name)
        else:
            var fz = as_zz(f, work)
            var gz = as_zz(g, work)
            var qz = as_zz(q, work)
            var rz = as_zz(r, work)
            var got = pseudo_divmod(DensePoly[ZZ].from_sparse(fz, work), DensePoly[ZZ].from_sparse(gz, work), work)
            c.check(got.quotient.to_sparse(order, work) == qz and got.remainder.to_sparse(order, work) == rz, record.name)
        count += 1
    c.check(count == 8, "all eight division corpus cases")
    c.finish()
