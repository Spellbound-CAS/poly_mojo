"""Ring laws, exact arithmetic, mixed moduli, and budget failures."""
from apn_mojo import Integer, Rational
from polynomials import Budget, PolyErrorKind, CoefficientRing, NoContext, Modulus, ZZ, QQ, ModInt
from support.check import Checker
from support.rng import SplitMix64

def laws[C: CoefficientRing](a: C, b: C, d: C, mut c: Checker) raises:
    var work = Budget()
    var z = C.zero(a.context())
    var one = C.one(a.context())
    c.check(z.context() == a.context() and one.context() == a.context(), "constant contexts")
    c.check(a.add(z, work) == a and a.mul(one, work) == a, "identities")
    c.check(a.add(a.neg(work), work) == z and a.sub(a, work) == z, "additive inverse")
    c.check(a.add(b, work) == b.add(a, work), "sum commutes")
    c.check(a.mul(b, work) == b.mul(a, work), "product commutes")
    c.check(a.add(b, work).add(d, work) == a.add(b.add(d, work), work), "sum associates")
    c.check(a.mul(b, work).mul(d, work) == a.mul(b.mul(d, work), work), "product associates")
    c.check(a.mul(b.add(d, work), work) == a.mul(b, work).add(a.mul(d, work), work), "distributivity")

def main() raises:
    var c = Checker("coefficients")
    var rng = SplitMix64(42)
    var work = Budget()
    var ctx = Modulus(4611686018427387847, work)
    var small = Modulus(101, work)
    var primes: List[UInt64] = [2, 3, 101, 4611686018427387847]
    for prime in primes:
        var boundary_ctx = Modulus(prime, work)
        var values: List[UInt64] = [0, 1, 2, prime - 1, prime, prime + 1, UInt64.MAX]
        for raw_a in values:
            for raw_b in values:
                var x = ModInt(raw_a, boundary_ctx)
                var y = ModInt(raw_b, boundary_ctx)
                var a = Integer(raw_a) % Integer(prime)
                var b = Integer(raw_b) % Integer(prime)
                var sum = x.add(y, work)
                var difference = x.sub(y, work)
                var product = x.mul(y, work)
                var negative = x.neg(work)
                c.check(Integer(sum.value) == (a + b) % Integer(prime), "modular boundary sum")
                c.check(Integer(difference.value) == (a + Integer(prime) - b) % Integer(prime), "modular boundary subtraction")
                c.check(Integer(product.value) == (a * b) % Integer(prime), "modular boundary product")
                c.check(Integer(negative.value) == (Integer(prime) - a) % Integer(prime), "modular boundary negation")
                c.check(sum.value < prime and difference.value < prime and product.value < prime and negative.value < prime, "canonical boundary residues")
    for i in range(80):
        var a = rng.below(2001) - 1000
        var b = rng.below(2001) - 1000
        var d = rng.below(2001) - 1000
        laws(ZZ(Integer(a)), ZZ(Integer(b)), ZZ(Integer(d)), c)
        laws(QQ(Rational(a, 17)), QQ(Rational(b, 23)), QQ(Rational(d, 31)), c)
        var x = ModInt(rng.next(), ctx)
        var y = ModInt(rng.next(), ctx)
        laws(x, y, ModInt(rng.next(), ctx), c)
        var product = x.mul(y, work)
        var expected = (Integer(x.value) * Integer(y.value)) % Integer(ctx.value)
        c.check(Integer(product.value) == expected, "128-bit modular product")
        if not x.is_zero():
            c.check(x.mul(x.inverse(work), work) == ModInt.one(ctx), "field inverse")
        if a != 0:
            var q = QQ(Rational(a, 17))
            c.check(q.mul(q.inverse(work), work) == QQ.one(NoContext()), "rational inverse")
        var other = ModInt(UInt64(i), small)
        c.check(other.add(ModInt.one(small), work).modulus == 101, "interleaved fields")
    var composites: List[UInt64] = [0, 1, 4, 561, 3215031751, 3825123056546413051, UInt64(1) << 62]
    for p in composites:
        try:
            var invalid = Modulus(p, work)
            c.check(False, String("accepted nonprime ", invalid.value))
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "reject invalid modulus")
    var a = ModInt.one(ctx)
    var b = ModInt.one(small)
    for op in range(3):
        try:
            if op == 0:
                _ = a.add(b, work)
            elif op == 1:
                _ = a.sub(b, work)
            else:
                _ = a.mul(b, work)
            c.check(False, "accepted mixed moduli")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "mixed field DOMAIN")
    try:
        _ = QQ.zero(NoContext()).inverse(work)
        c.check(False, "zero rational inverse")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero rational DOMAIN")
    try:
        _ = ModInt.zero(ctx).inverse(work)
        c.check(False, "zero modular inverse")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero modular DOMAIN")
    var empty = Budget(0)
    var large = ZZ(Integer(1) << 256)
    var original = large.copy()
    try:
        _ = large.mul(large, empty)
        c.check(False, "uncharged multiply")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "precharged coefficient work")
    c.check(large == original and empty.remaining() == 0, "failure preserves inputs")
    c.finish()
