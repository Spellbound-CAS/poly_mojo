"""Final-stack verification: exact geometry, adversarial data, and SIMD fallbacks."""
from apn_mojo import Integer, Rational
from polynomials import Budget, PolyErrorKind, DensePoly, QQ, NoContext
from polynomials._cad_final import _construct_stack, _verified_stack, _check_final_stack, _FinalStack
from polynomials._cad_final_signs import _signs_at, _SignStats
from polynomials._root_ops import _sign_at
from support.check import Checker


def polynomial(values: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var coefficients = List[QQ]()
    for value in values:
        coefficients.append(QQ(Rational(value)))
    return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work)


def rejected(inputs: List[DensePoly[QQ]], bad: _FinalStack, mut work: Budget) raises -> Bool:
    var scalar = _SignStats()
    var vector = _SignStats()
    return not _check_final_stack(inputs, bad, scalar, work, vectorized=False) and not _check_final_stack(inputs, bad, vector, work)


def signs(mut c: Checker, mut work: Budget) raises:
    # Every tail length, mixed denominators, and an independent exact oracle.
    var p = polynomial([6, -5, -2, 1], work)
    for count in range(1, 26):
        var points = List[QQ]()
        for i in range(count):
            points.append(QQ(Rational(i - 12, i % 5 + 1)))
        var a = _SignStats()
        var b = _SignStats()
        var scalar = _signs_at(p, points, a, work, vectorized=False)
        var vector = _signs_at(p, points, b, work)
        c.check(b.word_points == count and b.big_points == 0, "bounded batches use word arithmetic")
        for i in range(count):
            var expected = _sign_at(p, points[i], work)
            c.check(scalar[i] == expected and vector[i] == expected, "scalar and SIMD agree with arbitrary-precision oracle")

    # Integer products that would wrap to zero must instead take the fallback.
    var large = Integer(1) << 90
    var huge = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(large)), QQ(Rational(-1)), QQ(Rational(1))], work)
    var points: List[QQ] = [QQ(Rational(0)), QQ(Rational(1)), QQ(Rational(-1)), QQ(Rational(large)), QQ(Rational(1, large))]
    var stats = _SignStats()
    var result = _signs_at(huge, points, stats, work)
    c.check(stats.big_points == len(points), "wide coefficients use exact fallback")
    for i in range(len(points)):
        c.check(result[i] == _sign_at(huge, points[i], work), "wide homogeneous Horner remains exact")

    var square = polynomial([0, 0, 1], work)
    points = [QQ(Rational(Int64.MAX)), QQ(Rational(Int64.MIN)), QQ(Rational(1, Integer(1) << 63))]
    stats = _SignStats()
    result = _signs_at(square, points, stats, work)
    c.check(stats.big_points == 3, "product, denominator and signed-minimum overflow fall back")
    for sign in result:
        c.check(sign == 1, "overflow cannot produce a zero or negative square")

    var rational = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1, 3)), QQ(Rational(-1, 7))], work)
    points = [QQ(Rational(7, 3)), QQ(Rational(0)), QQ(Rational(3))]
    stats = _SignStats()
    result = _signs_at(rational, points, stats, work)
    c.check(result == [0, 1, -1], "positive denominator clearing preserves exact zero and signs")
    result = _signs_at(polynomial([], work), points, stats, work)
    c.check(result == [0, 0, 0], "zero polynomial signs")

    # Adjacent points on opposite sides of the checked Int64 product bound.
    for x in range(3037000498, 3037000502):
        points = [QQ(Rational(x)), QQ(Rational(-x))]
        stats = _SignStats()
        result = _signs_at(square, points, stats, work)
        c.check(result == [1, 1], "exact signs across the word-product boundary")
        c.check(stats.word_points == (2 if x <= 3037000499 else 0), "checked bound selects native or arbitrary precision")


def stacks(mut c: Checker, mut work: Budget) raises:
    var f = polynomial([-1, 0, 1], work)
    var repeated = polynomial([1, -2, 1], work)
    var irrational = polynomial([-2, 0, 1], work)
    var inputs: List[DensePoly[QQ]] = [f.copy(), repeated.copy(), irrational.copy(), f.mul(irrational, work), polynomial([-3], work), polynomial([], work)]
    var stack = _verified_stack(inputs, work)
    c.check(len(stack.roots) == 4 and len(stack.sectors) == 5 and len(stack.signs) == 9, "four distinct boundaries and all nine cells")
    var stats = _SignStats()
    c.check(_check_final_stack(inputs, stack, stats, work, vectorized=False), "scalar final verifier accepts independently")
    for j in range(len(inputs)):
        for i in range(len(stack.sectors)):
            c.check(stack.signs[2 * i][j] == _sign_at(inputs[j], stack.sectors[i], work), "sector signs equal direct input evaluation")
    for k in range(4):
        c.check(stack.signs[2 * k + 1][0] == (0 if k == 1 or k == 2 else 1), "x squared minus one section signs")
        c.check(stack.signs[2 * k + 1][1] == (0 if k == 2 else 1), "even multiplicity does not flip sign")
        c.check(stack.signs[2 * k + 1][2] == (0 if k == 0 or k == 3 else -1), "irrational section signs")
        c.check(stack.signs[2 * k + 1][3] == 0, "product vanishes at every section")

    var bad = stack.copy()
    bad.roots.clear()
    bad.sectors = [QQ(Rational(0))]
    bad.signs = [stack.signs[0].copy()]
    c.check(rejected(inputs, bad, work), "omitting all roots fails completeness even with consistent shapes")
    bad = stack.copy()
    bad.roots[0].basis = len(stack.basis)
    c.check(rejected(inputs, bad, work), "out-of-range root owner")
    bad = stack.copy()
    bad.roots[0].basis = -1
    c.check(rejected(inputs, bad, work), "negative root owner")
    bad = stack.copy()
    bad.roots[0] = stack.roots[1].copy()
    c.check(rejected(inputs, bad, work), "duplicate or reordered roots")
    bad = stack.copy()
    bad.roots[0].low = QQ(Rational(100))
    bad.roots[0].high = QQ(Rational(101))
    c.check(rejected(inputs, bad, work), "forged root bracket")
    bad = stack.copy()
    bad.sectors[0] = QQ(Rational(100))
    c.check(rejected(inputs, bad, work), "sector sample outside its claimed gap")
    bad = stack.copy()
    bad.signs[0][0] = -bad.signs[0][0]
    c.check(rejected(inputs, bad, work), "wrong sector sign")
    bad = stack.copy()
    bad.signs[3][0] = 1
    c.check(rejected(inputs, bad, work), "wrong section sign")
    bad = stack.copy()
    bad.signs[0][0] = 2
    c.check(rejected(inputs, bad, work), "invalid sign encoding")
    bad = stack.copy()
    bad.leading[0] = QQ(Rational(-1))
    c.check(rejected(inputs, bad, work), "forged leading sign")
    bad = stack.copy()
    bad.multiplicities[0][0] = Int.MAX
    c.check(rejected(inputs, bad, work), "huge multiplicity rejected before exponentiation")
    bad = stack.copy()
    bad.multiplicities[0][0] = -1
    c.check(rejected(inputs, bad, work), "negative multiplicity")
    bad = stack.copy()
    bad.basis[0] = bad.basis[0].mul(bad.basis[0], work)
    c.check(rejected(inputs, bad, work), "repeated-root basis rejected")
    bad = stack.copy()
    bad.basis[0].variable = 1
    c.check(rejected(inputs, bad, work), "different-variable basis rejected")
    bad = stack.copy()
    _ = bad.signs.pop()
    c.check(rejected(inputs, bad, work), "missing unbounded sector")

    # Exact rational sections are also valid final data, although the initial
    # producer isolates rational roots in open brackets.
    inputs = [polynomial([0, 1], work)]
    stack = _construct_stack(inputs, work)
    stack.roots[0].low = QQ(Rational(0))
    stack.roots[0].high = QQ(Rational(0))
    stats = _SignStats()
    c.check(_check_final_stack(inputs, stack, stats, work), "exact rational section checks")
    bad = stack.copy()
    bad.roots[0].low = QQ(Rational(1))
    bad.roots[0].high = QQ(Rational(1))
    c.check(rejected(inputs, bad, work), "forged exact rational root")
    bad = stack.copy()
    bad.sectors[0] = QQ(Rational(0))
    c.check(rejected(inputs, bad, work), "sector may not equal an exact section")

    inputs = [polynomial([5], work), polynomial([], work)]
    stack = _verified_stack(inputs, work)
    c.check(len(stack.roots) == 0 and stack.signs[0] == [1, 0], "constant and identically-zero inputs")
    stack = _verified_stack(List[DensePoly[QQ]](), work)
    c.check(len(stack.signs) == 1 and len(stack.signs[0]) == 0, "empty input covers the line")
    inputs = [polynomial([1, 0, 1], work)]
    stack = _verified_stack(inputs, work)
    c.check(len(stack.roots) == 0 and stack.signs[0] == [1], "nonconstant positive polynomial has no real sections")

    # One polynomial, with no extra inputs to split its multiplicity blocks.
    inputs = [polynomial([-2, 5, -4, 1], work)]  # (x-1)^2 (x-2)
    stack = _verified_stack(inputs, work)
    c.check(len(stack.roots) == 2 and len(stack.basis) == 2, "different multiplicities split a single input")
    c.check(stack.signs == [[-1], [0], [-1], [0], [1]], "even and odd multiplicities in a single input")

    inputs = [DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-(Integer(1) << 90))), QQ(Rational(0)), QQ(Rational(1))], work)]
    stack = _verified_stack(inputs, work)
    stats = _SignStats()
    c.check(_check_final_stack(inputs, stack, stats, work) and stats.big_points > 0, "whole-stack arbitrary-precision root validation")
    c.check(stack.signs == [[1], [0], [-1], [0], [1]], "wide quadratic final signs")

    var limited = Budget(0)
    var saw_limit = False
    try:
        _ = _check_final_stack(inputs, stack, stats, limited)
    except e:
        saw_limit = e.kind == PolyErrorKind.LIMIT
    c.check(saw_limit, "verification exhaustion propagates as LIMIT")

    # Exhaustion late in verification must not produce an accepted result.
    for mode in range(2):
        var measured = Budget()
        stats = _SignStats()
        c.check(_check_final_stack(inputs, stack, stats, measured, vectorized=Bool(mode)), "measure exact verification budget")
        var needed = UInt64.MAX - measured.remaining()
        limited = Budget(needed - 1)
        saw_limit = False
        stats = _SignStats()
        try:
            _ = _check_final_stack(inputs, stack, stats, limited, vectorized=Bool(mode))
        except e:
            saw_limit = e.kind == PolyErrorKind.LIMIT
        c.check(saw_limit, "one step short cannot return an unchecked result")


def main() raises:
    var c = Checker("final CAD stack")
    var work = Budget()
    signs(c, work)
    stacks(c, work)
    c.finish()
