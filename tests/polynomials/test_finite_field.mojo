"""Corpus factors and exhaustive small-field irreducibles across explicit seeds."""
from polynomials import Budget, PolyErrorKind, ModInt, Modulus, DensePoly, MonomialOrder, factor_mod_p, divmod
from support.corpus import read_cases, parse_polynomial, as_mod
from support.check import Checker


def exhaustive(prime: UInt64, mut c: Checker) raises:
    var work = Budget()
    var ctx = Modulus(prime, work)
    var irreducibles = List[DensePoly[ModInt]]()
    var extent = Int(prime)
    for degree in range(1, 4):
        for code in range(extent):
            var values = List[ModInt]()
            var digits = code
            for _ in range(degree):
                values.append(ModInt(UInt64(digits % Int(prime)), ctx))
                digits //= Int(prime)
            values.append(ModInt.one(ctx))
            var candidate = DensePoly[ModInt].from_coefficients(0, ctx, values^, work)
            var irreducible = True
            # Every reducible quadratic or cubic has a linear factor.
            if degree > 1:
                for root in range(Int(prime)):
                    var value = UInt64(0)
                    for i in range(degree, -1, -1):
                        value = (value * UInt64(root) + candidate.coefficients[i].value) % prime
                    if value == 0:
                        irreducible = False
            if irreducible:
                irreducibles.append(candidate^)
        extent *= Int(prime)
    for i in range(len(irreducibles)):
        var j = (i + 1) % len(irreducibles)
        var product = irreducibles[i].mul(irreducibles[j], work)
        for seed in range(3):
            var got = factor_mod_p(product, work, UInt64(seed))
            c.check(len(got) == 2, "two planted irreducibles")
            c.check((got[0] == irreducibles[i] and got[1] == irreducibles[j]) or (got[1] == irreducibles[i] and got[0] == irreducibles[j]), "exhaustively verified irreducibles recovered")
    var one = DensePoly[ModInt].from_coefficients(0, ctx, [ModInt.one(ctx)], work)
    c.check(len(factor_mod_p(one, work, 0)) == 0, "one has no factors")
    var repeated = irreducibles[0].pow(2, work)
    try:
        _ = factor_mod_p(repeated, work, 0)
        c.check(False, "repeated factor accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "repeated factor DOMAIN")
    var zero = DensePoly[ModInt](0, ctx)
    try:
        _ = factor_mod_p(zero, work, 0)
        c.check(False, "zero accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero DOMAIN")
    if prime > 2:
        var nonmonic = one.scale(ModInt(2, ctx), work)
        try:
            _ = factor_mod_p(nonmonic, work, 0)
            c.check(False, "nonmonic accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "nonmonic DOMAIN")
    var before = repeated.copy()
    var tiny = Budget(5)
    try:
        _ = factor_mod_p(repeated, tiny, 0)
        c.check(False, "budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and repeated == before, "finite field budget preserves input")


def main() raises:
    var work = Budget()
    var c = Checker("finite_field")
    var count = 0
    for record in read_cases():
        if record.field("op ") != "factor_mod_p":
            continue
        var ring = record.field("ring ")
        var ctx = Modulus(UInt64(Int(ring[byte=3:])), work)
        var input_text = record.field("in f ")
        var parsed = parse_polynomial(input_text, MonomialOrder(), work)
        var f = DensePoly[ModInt].from_sparse(as_mod(parsed, ctx, work), work)
        var expected = List[DensePoly[ModInt]]()
        for line in record.lines:
            if line.startswith("out factor "):
                var text = String(line[byte=11:])
                var p = parse_polynomial(text, MonomialOrder(), work)
                expected.append(DensePoly[ModInt].from_sparse(as_mod(p, ctx, work), work))
        for seed in range(4):
            var got = factor_mod_p(f, work, UInt64(seed))
            c.check(got == expected, record.name)
            var product = DensePoly[ModInt].from_coefficients(0, ctx, [ModInt.one(ctx)], work)
            for g in got:
                product = product.mul(g, work)
                g.check_invariants(work)
            c.check(product == f, "factor product")
        count += 1
    c.check(count == 11, "all eleven finite-field cases")
    exhaustive(2, c)
    exhaustive(3, c)
    exhaustive(5, c)
    c.finish()
