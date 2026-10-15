"""Exact regression identities for shared reduction, grouping and CAD work."""
from apn_mojo import Rational
from polynomials import *
from polynomials._sparse_modular import _coefficients, _at, _support_interpolate
from polynomials._views import _split, _join
from polynomials._sparse_ops import _subtract_multiple
from polynomials._grouping import _coefficient_groups
from polynomials._powers import _powers
from polynomials._sampling import _SamplingPlan
from polynomials._sparse_gcd import _interpolate
from polynomials.groebner import _Pair, _heapify_pairs, _select
from support.check import Checker
from support.corpus import parse_polynomial, as_mod, as_zz


def division[C: CoefficientField](f: SparsePoly[C], g: SparsePoly[C], quotient: SparsePoly[C], remainder: SparsePoly[C], mut c: Checker, mut work: Budget) raises:
    var divisors: List[SparsePoly[C]] = [g.copy()]
    var result = divide_by_list(f, divisors, work)
    c.check(result.quotients[0] == quotient and result.remainder == remainder, "hand-derived ordered division")
    var exact = f.sub(remainder, work)
    var q = divide_exact(exact, g, work)
    c.check(Bool(q) and q.value() == quotient, "hand-derived exact quotient")
    c.check(not divide_exact(f, g, work), "nonzero remainder rejects exact division")
    var saved = f.copy()
    var saved_divisors = divisors.copy()
    for limit in [UInt64(0), UInt64(50), UInt64(500)]:
        var tiny = Budget(limit)
        try:
            var limited = divide_by_list(f, divisors, tiny)
            c.check(limited.quotients[0] == quotient and limited.remainder == remainder, "limited reduction completes exactly")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "limited reduction fails closed")
        c.check(f == saved and divisors == saved_divisors, "reduction preserves borrowed inputs")


def arithmetic(mut c: Checker, mut work: Budget) raises:
    for kind in range(1, 4):
        var order = MonomialOrder(UInt8(kind), [Int64(0), Int64(1), Int64(2)], work)
        var f = parse_polynomial("x0^2 + x0*x1 - x0 + 1", order, work)
        var g = parse_polynomial("x0 - 1", order, work)
        var q = parse_polynomial("x0 + x1", order, work)
        var r = parse_polynomial("x1 + 1", order, work)
        division(f, g, q, r, c, work)
        var ctx = Modulus(7, work)
        division(as_mod(f, ctx, work), as_mod(g, ctx, work), as_mod(q, ctx, work), as_mod(r, ctx, work), c, work)
        var zq = divide_exact(as_zz(f.sub(r, work), work), as_zz(g, work), work)
        c.check(Bool(zq) and zq.value() == as_zz(q, work), "integer exact division")
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(1), Int64(0), Int64(2)], work)
        var f = parse_polynomial("x0*x1^6 + 2*x0*x1^3 + x0 + x1^3", order, work)
        var q = parse_polynomial("x2 + 1", order, work)
        var expected = parse_polynomial("x0*x2^6 + 6*x0*x2^5 + 15*x0*x2^4 + 22*x0*x2^3 + 21*x0*x2^2 + 12*x0*x2 + 4*x0 + x2^3 + 3*x2^2 + 3*x2 + 1", order, work)
        c.check(substitute(f, 1, q, work) == expected, "sparse grouped binomial substitution")
        c.check(substitute(f, 2, q, work) == f, "absent substituted variable")
        var zero = SparsePoly[QQ].zero(order, NoContext(), work)
        c.check(substitute(f, 1, zero, work) == parse_polynomial("x0", order, work), "zero substitution retains constant coefficient")
        var parts = _split(f, 1, work)
        c.check(len(parts) == 7 and parts[1].is_zero() and parts[5].is_zero(), "dense coefficient adapter retains gaps")
        c.check(_join(parts, 1, f, work) == f, "split-join identity in every order")
        var grouped = _coefficients(f, 1, False, work)
        c.check(len(grouped) == 3, "occupied coefficient groups omit gaps")
        var saved = f.copy()
        for limit in [UInt64(0), UInt64(100), UInt64(1000)]:
            var tiny = Budget(limit)
            try:
                c.check(substitute(f, 1, q, tiny) == expected, "limited substitution completes exactly")
            except e:
                c.check(e.kind == PolyErrorKind.LIMIT, "substitution fails closed")
            c.check(f == saved, "substitution preserves input")


def prepared_helpers(mut c: Checker, mut work: Budget) raises:
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(0), Int64(1)], work)
        var f = parse_polynomial("x0^2 + x1 + 1", order, work)
        var g = parse_polynomial("x0 + x1", order, work)
        var term = Term[QQ](Monomial([(Int64(0), UInt32(1))], work), QQ(Rational(2)))
        var expected = parse_polynomial("-x0^2 - 2*x0*x1 + x1 + 1", order, work)
        c.check(_subtract_multiple(f, g, term, work) == expected, "fused hand-derived sparse subtraction")
        var zero_term = Term[QQ](term.monomial.copy(), QQ.zero(NoContext()))
        c.check(_subtract_multiple(f, g, zero_term, work) == f, "zero multiple preserves the polynomial")
        var ctx = Modulus(5, work)
        var modular_term = Term[ModInt](term.monomial.copy(), ModInt(2, ctx))
        c.check(_subtract_multiple(as_mod(f, ctx, work), as_mod(g, ctx, work), modular_term, work) == as_mod(expected, ctx, work), "fused modular subtraction")
        var foreign = Term[ModInt](term.monomial.copy(), ModInt(2, Modulus(7, work)))
        var modular_f = as_mod(f, ctx, work)
        var modular_g = as_mod(g, ctx, work)
        try:
            _ = _subtract_multiple(modular_f, modular_g, foreign, work)
            c.check(False, "foreign multiplier context accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "foreign multiplier is DOMAIN")
        var saved = f.copy()
        for limit in [UInt64(0), UInt64(15), UInt64(60)]:
            var tiny = Budget(limit)
            try:
                c.check(_subtract_multiple(f, g, term, tiny) == expected, "limited fused merge completes")
            except e:
                c.check(e.kind == PolyErrorKind.LIMIT, "fused merge fails closed")
            c.check(f == saved, "fused merge preserves input")
        var high = SparsePoly[QQ].from_terms([Term[QQ](Monomial([(Int64(0), UInt32.MAX)], work), QQ.one(NoContext()))], order, NoContext(), work)
        try:
            _ = _subtract_multiple(f, high, term, work)
            c.check(False, "fused exponent overflow accepted")
        except e:
            c.check(e.kind == PolyErrorKind.CAPACITY, "fused exponent overflow is CAPACITY")
        var sparse = high.add(SparsePoly[QQ].constant(QQ.one(NoContext()), order, work), work)
        c.check(len(_coefficient_groups(sparse, 0, False, work)[0]) == 2, "large sparse exponent allocates only occupied groups")
        var tiny = Budget(10000)
        c.check(substitute(sparse, 0, SparsePoly[QQ].constant(QQ.one(NoContext()), order, work), tiny) == SparsePoly[QQ].constant(QQ(Rational(2)), order, work), "large exponent substitution avoids a degree-sized array")
    var powers = _powers(QQ(Rational(-2, 3)), work)
    c.check(powers.power(0, work) == QQ.one(NoContext()), "zeroth cached power")
    c.check(powers.power(7, work) == QQ(Rational(-128, 2187)), "exact rational seventh power")
    c.check(powers.power(2, work) == QQ(Rational(4, 9)), "out-of-order exponent cache insertion")
    var first = Budget()
    _ = powers.power(19, first)
    var reused = Budget()
    _ = powers.power(19, reused)
    c.check(reused.remaining() > first.remaining(), "repeated exponent reuses arithmetic")
    var exhausted = Budget(0)
    try:
        _ = powers.power(31, exhausted)
        c.check(False, "unbudgeted cached power accepted")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and powers.power(7, work) == QQ(Rational(-128, 2187)), "cache remains valid after exhaustion")
    var expected_power = QQ(Rational(4294967296, 1853020188851841))
    for limit in [UInt64(1), UInt64(4), UInt64(32), UInt64(128), UInt64(512)]:
        var interrupted = _powers(QQ(Rational(-2, 3)), work)
        _ = interrupted.power(7, work)
        var limited = Budget(limit)
        try:
            c.check(interrupted.power(32, limited) == expected_power, "limited square-cache extension completes exactly")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "square-cache extension fails closed")
        c.check(interrupted.power(7, work) == QQ(Rational(-128, 2187)) and interrupted.power(32, work) == expected_power, "cache remains usable after partial square growth")
    var mod5 = _powers(ModInt(3, Modulus(5, work)), work)
    var mod7 = _powers(ModInt(3, Modulus(7, work)), work)
    c.check(mod5.power(4, work) == ModInt(1, Modulus(5, work)) and mod7.power(4, work) == ModInt(4, Modulus(7, work)), "power caches isolate coefficient contexts")
    for kind in range(1, 4):
        var order = MonomialOrder(UInt8(kind), [Int64(0)], work)
        var pairs = List[_Pair]()
        for i in range(32, -1, -1):
            var monomial = Monomial([(Int64(0), UInt32(i // 8))], work) if i >= 8 else Monomial()
            pairs.append(_Pair((i % 8) // 2, 100 + i, monomial^))
        _heapify_pairs(pairs, order, work)
        for i in range(33):
            var pair = _select(pairs, order, work)
            c.check(pair.first == (i % 8) // 2 and pair.second == 100 + i and pair.lcm.exponent(0, work) == UInt32(i // 8), "heap preserves monomial and both index tie breaks")
        c.check(len(pairs) == 0, "pair heap drains exactly once")


def interpolation(mut c: Checker, mut work: Budget) raises:
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0), Int64(1)], work)
    var ctx = Modulus(101, work)
    var f = as_mod(parse_polynomial("x0^3*x1^2 + 2*x0^3 + 3*x0*x1 + 4*x0 + 5*x1^2 + 6", order, work), ctx, work)
    var support = List[Monomial]()
    for term in f.terms:
        support.append(term.monomial.copy())
    var points: List[ModInt] = [ModInt(0, ctx), ModInt(1, ctx), ModInt(2, ctx)]
    var images = List[SparsePoly[ModInt]]()
    for point in points:
        images.append(_at(f, 1, point, work))
    var result = _support_interpolate(support, 1, points, images, work)
    c.check(Bool(result) and result.value() == f, "interpolation identity with overlapping exponents")
    var repeated: List[ModInt] = [ModInt(0, ctx), ModInt(0, ctx), ModInt(0, ctx)]
    var repeated_images: List[SparsePoly[ModInt]] = [images[0].copy(), images[0].copy(), images[0].copy()]
    c.check(not _support_interpolate(support, 1, repeated, repeated_images, work), "singular sampling is inconclusive")
    try:
        var missing = _SamplingPlan(f.terms, 0, [Int64(0)], work)
        c.check(False, "incomplete sampling coordinates accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "missing sampling coordinate is DOMAIN")
    var univariate = as_mod(parse_polynomial("x0 + 1", order, work), ctx, work)
    c.check(not _interpolate(univariate, univariate, support, 0, 1, work), "support with a vanished variable requests rediscovery")
    c.check(not _interpolate(univariate, univariate, List[Monomial](), 0, 1, work), "empty interpolation support is inconclusive")


def proofs(mut c: Checker, mut work: Budget) raises:
    var order = MonomialOrder()
    var f = parse_polynomial("x1^4 + x0*x1 + 2", order, work)
    var g = parse_polynomial("x1^2 - 1", order, work)
    var expected = parse_polynomial("9 - x0^2", order, work)
    var proof = projection_resultant(f, g, 1, work)
    c.check(proof.resultant == expected and check_projection_resultant(f, g, proof, work), "resultant equals f(1)*f(-1)")
    var saved = proof.copy()
    c.check(check_projection_resultant(f, g, proof, work) and proof.resultant == saved.resultant, "determinant checker preserves its input")
    for i in range(len(proof.steps)):
        c.check(proof.steps[i].trailing == saved.steps[i].trailing and proof.steps[i].pivot_row == saved.steps[i].pivot_row, "checker leaves every transcript matrix unchanged")
    var bad = proof.copy()
    bad.steps[0].pivot_row = len(proof.steps) + 1
    c.check(not check_projection_resultant(f, g, bad, work), "invalid pivot rejected")
    bad = proof.copy()
    bad.steps[0].pivot_row = -2
    c.check(not check_projection_resultant(f, g, bad, work), "negative pivot other than singular marker rejected")
    bad = proof.copy()
    _ = bad.steps[0].trailing.pop()
    c.check(not check_projection_resultant(f, g, bad, work), "missing trailing row rejected")
    bad = proof.copy()
    _ = bad.steps[0].trailing[0].pop()
    c.check(not check_projection_resultant(f, g, bad, work), "ragged trailing matrix rejected")
    bad = proof.copy()
    _ = bad.steps.pop()
    c.check(not check_projection_resultant(f, g, bad, work), "unfinished determinant rejected")
    bad = proof.copy()
    bad.steps.append(proof.steps[len(proof.steps) - 1].copy())
    c.check(not check_projection_resultant(f, g, bad, work), "step beyond terminal determinant rejected")
    bad = proof.copy()
    bad.steps[0].pivot_row = -1
    bad.steps[0].trailing = List[List[SparsePoly[QQ]]]()
    c.check(not check_projection_resultant(f, g, bad, work), "false singular first column rejected")
    var constant = parse_polynomial("2", order, work)
    var other = parse_polynomial("3", order, work)
    var empty = projection_resultant(constant, other, 1, work)
    c.check(len(empty.steps) == 0 and empty.resultant == parse_polynomial("1", order, work) and check_projection_resultant(constant, other, empty, work), "empty Sylvester determinant is one")
    var linear = parse_polynomial("x1 + x0", order, work)
    var single = projection_resultant(constant, linear, 1, work)
    c.check(len(single.steps) == 0 and single.resultant == constant and check_projection_resultant(constant, linear, single, work), "single-entry Sylvester determinant")
    var limited = Budget(0)
    try:
        _ = check_projection_resultant(f, g, proof, limited)
        c.check(False, "unbudgeted determinant checking accepted")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and proof.steps[0].trailing == saved.steps[0].trailing, "determinant exhaustion preserves borrowed transcript")
    var repeated = parse_polynomial("x0^2 + x0 + 1", order, work)
    var inputs: List[SparsePoly[QQ]] = [repeated.copy(), repeated.copy(), repeated.copy()]
    var projection = cad_projection(inputs, [Int64(0)], work)
    c.check(len(projection.factorizations) == 3 and check_cad_projection(inputs, projection, work), "repeated projection retains every scheduled witness")
    var factor = projection.factorizations[0].factors[0][0].copy()
    projection.factorizations[1].factors[0][0].terms[0].coefficient = QQ(Rational(2))
    c.check(projection.factorizations[0].factors[0][0] == factor and projection.factorizations[2].factors[0][0] == factor, "cached factor witnesses own independent polynomials")
    c.check(not check_cad_projection(inputs, projection, work), "corrupted copied factor rejected")
    projection.factorizations[1].factors[0][0] = factor^
    projection.factorizations[1].unit = QQ(Rational(2))
    c.check(not check_cad_projection(inputs, projection, work), "corrupted repeated witness rejected")


def main() raises:
    var c = Checker("prepared operations")
    var work = Budget()
    arithmetic(c, work)
    prepared_helpers(c, work)
    interpolation(c, work)
    proofs(c, work)
    c.finish()
