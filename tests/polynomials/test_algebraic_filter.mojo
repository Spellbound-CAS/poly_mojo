"""Fast field routines of the algebraic filters agree with their exact definitions."""
from apn_mojo import Integer, Rational
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, SparsePoly, RealAlgebraic, real_root, real_rational
from polynomials._algebraic import _Frame, _NormalizationCache
from polynomials._root_ops import _compare, _sign_at
from polynomials.resultant import resultant
from polynomials.certificates import _dense_equal
from polynomials._algebraic_filter import _norm_image, _gcd_modular, _gcd_euclid, _gcd_over_field, _reduce, _frame_value, _check_frame_direct
from support.check import Checker


def q(n: Int, d: Int = 1) raises -> QQ:
    return QQ(Rational(Integer(n), Integer(d)))


def poly(values: List[QQ], mut work: Budget) raises -> DensePoly[QQ]:
    return DensePoly[QQ].from_coefficients(0, NoContext(), values.copy(), work)


def constant(value: QQ, mut work: Budget) raises -> DensePoly[QQ]:
    return poly([value.copy()], work)


def embed(p: DensePoly[QQ], variable: Int64, mut work: Budget) raises -> SparsePoly[QQ]:
    var moved = p.copy()
    moved.variable = variable
    return moved.to_sparse(MonomialOrder(), work)


def same(a: List[DensePoly[QQ]], b: List[DensePoly[QQ]], mut work: Budget) raises -> Bool:
    if len(a) != len(b):
        return False
    for i in range(len(a)):
        if not _dense_equal(a[i], b[i], work):
            return False
    return True


def product(a: List[DensePoly[QQ]], b: List[DensePoly[QQ]], m: DensePoly[QQ], mut work: Budget) raises -> List[DensePoly[QQ]]:
    var result = List[DensePoly[QQ]]()
    for _ in range(len(a) + len(b) - 1):
        result.append(DensePoly[QQ](0, NoContext()))
    for i in range(len(a)):
        for j in range(len(b)):
            result[i + j] = _reduce(result[i + j].add(a[i].mul(b[j], work), work), m, work)
    return result^


def main() raises:
    var work = Budget()
    var c = Checker("algebraic filter fields")
    # m is the minimal polynomial of b = sqrt(2) + sqrt(3); a non-monic cubic as well.
    var m = poly([q(1), q(0), q(-10), q(0), q(1)], work)
    var cubic = poly([q(-5), q(0), q(0), q(2)], work)
    var order = MonomialOrder()
    var images: List[DensePoly[QQ]] = [poly([q(1), q(0), q(1)], work), poly([q(1, 2), q(-1), q(0), q(3)], work), poly([q(0), q(7, 3)], work)]
    for field in [m.copy(), cubic.copy()]:
        for r in images:
            var reduced = _reduce(r, field, work)
            if reduced.degree() < 1:
                continue
            var shifted = SparsePoly[QQ].variable(0, order, NoContext(), work).sub(embed(reduced, 1, work), work)
            var expected = DensePoly[QQ].from_sparse(resultant(embed(field, 1, work), shifted, 1, work).resultant, work, variable=Int64(0))
            c.check(_dense_equal(_norm_image(reduced, field, work), expected, work), "trace norm equals the resultant")
    # sqrt(2) = (b^3 - 9b)/2 in Q(b).
    var s = poly([q(0), q(-9, 2), q(0), q(1, 2)], work)
    var minus_s = DensePoly[QQ](0, NoContext()).sub(s, work)
    var two: List[DensePoly[QQ]] = [constant(q(-2), work), constant(q(0), work), constant(q(1), work)]
    var linear: List[DensePoly[QQ]] = [minus_s.copy(), constant(q(1), work)]
    var other: List[DensePoly[QQ]] = [constant(q(-7, 999999999999999989), work), constant(q(1), work)]
    var big: List[DensePoly[QQ]] = [poly([q(3), q(1, 1000000007)], work), constant(q(5), work), constant(q(1), work)]
    var cases = List[Tuple[List[DensePoly[QQ]], List[DensePoly[QQ]], Int]]()
    cases.append((two.copy(), product(linear, other, m, work), 1))
    cases.append((product(two, other, m, work), product(two, big, m, work), 2))
    cases.append((product(linear, big, m, work), product(linear, other, m, work), 1))
    cases.append((big.copy(), other.copy(), 0))
    for item in cases:
        var exact = _gcd_euclid(item[0], item[1], m, work)
        var fast = _gcd_modular(item[0], item[1], m, work)
        c.check(len(exact) == item[2] + 1, "exact gcd degree")
        c.check(Bool(fast) and same(fast.value(), exact, work), "modular gcd equals the exact gcd")
        c.check(same(_gcd_over_field(item[0], item[1], m, work), exact, work), "dispatch equals the exact gcd")
    # A point (sqrt 2, sqrt 3, 1/2) framed by theta = sqrt 2 + sqrt 3.
    var root2 = real_root(poly([q(-2), q(0), q(1)], work), q(1), q(2), work)
    var root3 = real_root(poly([q(-3), q(0), q(1)], work), q(1), q(2), work)
    var theta = real_root(m, q(3), q(4), work)
    var point: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), root2.copy()), (Int64(1), root3.copy()), (Int64(2), real_rational(q(1, 2), work))]
    var frame = _Frame()
    frame.generator = theta.copy()
    frame.variables = [Int64(0), Int64(1), Int64(2)]
    frame.representations = [s.copy(), poly([q(0), q(11, 2), q(0), q(-1, 2)], work), constant(q(1, 2), work)]
    var x0 = SparsePoly[QQ].variable(0, order, NoContext(), work)
    var x1 = SparsePoly[QQ].variable(1, order, NoContext(), work)
    var x2 = SparsePoly[QQ].variable(2, order, NoContext(), work)
    var one = SparsePoly[QQ].constant(q(1), order, work)
    var cache = _NormalizationCache()
    # x0 x1 - 2 x2 = sqrt 6 - 1, and x0^2 + x1^2 - 10 x2 = 0.
    var f = x0.mul(x1, work).sub(x2.add(x2, work), work)
    var framed = _frame_value(f, point, frame, cache, work)
    c.check(Bool(framed), "framed evaluation found")
    var parts = framed.take()
    c.check(not parts[0].rational and parts[0].polynomial.degree() == 2 and _compare(parts[0].low, q(1), work) >= 0 and _compare(parts[0].high, q(2), work) <= 0 and _sign_at(poly([q(-5), q(2), q(1)], work), parts[0].low, work) != _sign_at(poly([q(-5), q(2), q(1)], work), parts[0].high, work), "framed value is sqrt 6 - 1")
    c.check(_check_frame_direct(f, point, parts[0], [parts[1].copy()], parts[2], work), "framed evaluation accepted")
    var swapped = parts[2].copy()
    swapped[0] = parts[2][1].copy()
    swapped[1] = parts[2][0].copy()
    c.check(not _check_frame_direct(f, point, parts[0], [parts[1].copy()], swapped, work), "swapped representations rejected")
    var wrong = real_rational(q(1), work)
    c.check(not _check_frame_direct(f, point, wrong, [parts[1].copy()], parts[2], work), "wrong framed value rejected")
    var g = x0.mul(x0, work).add(x1.mul(x1, work), work).sub(x2.mul(SparsePoly[QQ].constant(q(10), order, work), work), work)
    var zero = _frame_value(g, point, frame, cache, work)
    c.check(Bool(zero) and zero.value()[0].rational and zero.value()[0].low.is_zero(), "framed zero is exact")
    _ = one
    c.finish()
