"""Evaluation homomorphisms, substitution, and the formal product rule."""
from polynomials import Budget, PolyErrorKind, CoefficientRing, ZZ, QQ, ModInt, Modulus, NoContext, MonomialOrder, SparsePoly, DensePoly, evaluate, evaluate_dense, derivative, substitute
from support.check import Checker
from support.gen import polynomial, coefficient
from support.rng import SplitMix64
from apn_mojo import Integer

def laws[C: CoefficientRing](ctx: C.Context, mut c: Checker) raises:
    var work = Budget()
    var rng = SplitMix64(411)
    var order = MonomialOrder()
    var x = SparsePoly[C].variable(0, order, ctx, work)
    var y = SparsePoly[C].variable(1, order, ctx, work)
    for _ in range(30):
        var a = polynomial[C](ctx, order, rng, work, variables=2)
        var b = polynomial[C](ctx, order, rng, work, variables=2)
        var point = List[Tuple[Int64, C]]()
        point.append((Int64(0), coefficient[C](ctx, rng.below(9) - 4, work)))
        point.append((Int64(1), coefficient[C](ctx, rng.below(9) - 4, work)))
        var av = evaluate(a, point, work)
        var bv = evaluate(b, point, work)
        c.check(evaluate(a.add(b, work), point, work) == av.add(bv, work), "sum homomorphism")
        c.check(evaluate(a.mul(b, work), point, work) == av.mul(bv, work), "product homomorphism")
        var da = derivative(a, 0, work)
        var db = derivative(b, 0, work)
        da.check_invariants(work)
        c.check(derivative(a.mul(b, work), 0, work) == da.mul(b, work).add(a.mul(db, work), work), "product rule")
        c.check(derivative(a.add(b, work), 0, work) == da.add(db, work), "derivative linearity")
        var replacement = x.add(y, work)
        c.check(substitute(a.add(b, work), 0, replacement, work) == substitute(a, 0, replacement, work).add(substitute(b, 0, replacement, work), work), "substitution sum")
        c.check(substitute(a, 0, x, work) == a, "identity substitution")
        c.check(substitute(a, 9, replacement, work) == a, "absent variable substitution")
    var f = x.pow(3, work).add(x, work)
    var dense = DensePoly[C].from_sparse(f, work)
    var three = coefficient[C](ctx, 3, work)
    var point = List[Tuple[Int64, C]]()
    point.append((Int64(0), three.copy()))
    c.check(evaluate(f, point, work) == coefficient[C](ctx, 30, work), "hand-computed cubic")
    c.check(evaluate_dense(dense, three, work) == evaluate(f, point, work), "dense Horner agrees")
    c.check(derivative(x.pow(3, work), 0, work) == x.pow(2, work).scale(three, work), "power rule")
    c.check(substitute(f, 0, SparsePoly[C].constant(three, order, work), work) == SparsePoly[C].constant(coefficient[C](ctx, 30, work), order, work), "constant substitution")
    try:
        _ = evaluate(y, point, work)
        c.check(False, "missing point accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "missing point DOMAIN")
    point.append((Int64(0), three.copy()))
    try:
        _ = evaluate(f, point, work)
        c.check(False, "duplicate point accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "duplicate point DOMAIN")
    var zero_budget = Budget(0)
    try:
        _ = derivative(f, 0, zero_budget)
        c.check(False, "uncharged derivative")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "derivative budget")

def word_blocks(mut c: Checker) raises:
    """Integer points: the word-level block Horner under the balanced halves agrees with plain Horner at sign and size extremes."""
    var work = Budget()
    var rng = SplitMix64(977)
    var top = Integer(Int64.MAX)
    var points = [Integer(0), Integer(1), Integer(-1), Integer(3), Integer(-7), Integer(2147483647), Integer(-2147483000), Integer(1) << Integer(62), -(Integer(1) << Integer(62)), top, -top, top + Integer(1)]
    for length in [31, 32, 33, 63, 64, 65, 100, 257]:
        var coefficients = List[ZZ]()
        for i in range(length):
            var r = rng.next()
            var v = Integer(Int64(r >> 1)) if r % 5 else (top if i % 2 else -top)
            coefficients.append(ZZ(-v if r % 3 == 0 else v))
        # Only the leading coefficient must be nonzero; zeros elsewhere exercise empty blocks.
        for i in range(min(length - 1, 40)):
            coefficients[i] = ZZ(Integer(0)) if i % 9 == 4 else coefficients[i].copy()
        var p = DensePoly[ZZ].from_coefficients(0, NoContext(), coefficients.copy(), work)
        for x in points:
            var expected = Integer(0)
            for i in range(len(p.coefficients) - 1, -1, -1):
                expected = expected * x + p.coefficients[i].value
            c.check(evaluate_dense(p, ZZ(x.copy()), work).value == expected, String("integer evaluation at ", x, " of ", length, " coefficients"))


def main() raises:
    var c = Checker("evaluation")
    var work = Budget()
    laws[ZZ](NoContext(), c)
    laws[QQ](NoContext(), c)
    laws[ModInt](Modulus(3, work), c)
    word_blocks(c)
    c.finish()
