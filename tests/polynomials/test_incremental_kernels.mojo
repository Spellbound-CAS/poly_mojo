"""Algebraic identities for bounded polynomial intermediates."""
from apn_mojo import Rational, Integer, ExactComplex
from polynomials import *
from polynomials._bounded_products import _truncated_product, _PolynomialSum
from polynomials._grouping import _CoefficientGroups, _coefficient_groups
from polynomials._product_stream import _subtract_column
from support.check import Checker
from support.products import geometric


def domains[C: CoefficientRing](one: C, mut c: Checker, mut work: Budget) raises:
    var order = MonomialOrder(MonomialOrder.GREVLEX, [Int64(0), Int64(1)], work)
    var a = geometric(17, 0, False, one, order, work)
    var b = geometric(17, 1, False, one, order, work)
    var unit = SparsePoly[C].constant(one, order, work)
    var zero = SparsePoly[C].zero(order, one.context(), work)
    var quotients: List[SparsePoly[C]] = [a.copy(), a.copy(), unit.copy(), zero.copy()]
    var rows: List[List[SparsePoly[C]]] = [[b.copy()], [b.neg(work)], [b.copy()], [b.copy()]]
    c.check(_subtract_column(b, quotients, [0, 1, 2, 3], rows, 0, work).is_zero(), "domain-generic signed products with an empty group")
    quotients.append(unit.copy())
    rows.append([zero.copy()])
    c.check(_subtract_column(b, quotients, [0, 1, 2, 3, 4], rows, 0, work).is_zero(), "outer-heap signed products with empty columns")
    var trimmed = _truncated_product(a, b, 0, 0, work)
    c.check(trimmed == a.scale(one, work), "only parameter-constant pairs survive")
    var small = geometric(16, 1, False, one, order, work)
    var leading = SparsePoly[C].variable(1, order, one.context(), work).pow(30, work)
    c.check(_truncated_product(small, small, 0, 29, work) == small.mul(small, work).sub(leading, work), "tiny dense cut removes exactly the monic leading term")
    var limited = Budget(800)
    try:
        _ = _truncated_product(small, small, 0, 29, limited)
        c.check(False, "tiny dense cut budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "tiny dense cut fails closed")
    c.check(small == geometric(16, 1, False, one, order, work), "tiny dense cut preserves borrowed operands")
    var parameter = SparsePoly[C].variable(1, order, one.context(), work)
    var mixed = geometric(15, 0, False, one, order, work).add(parameter, work)
    c.check(_truncated_product(mixed, mixed, 0, 1, work) == mixed.mul(mixed, work).sub(parameter.pow(2, work), work), "owned compaction removes an interior term without reordering later terms")
    for allowance in [0, 32, 500]:
        var tiny = Budget(UInt64(allowance))
        try:
            _ = _subtract_column(b, quotients, [0, 1, 2, 3, 4], rows, 0, tiny)
            c.check(False, "signed product budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "signed product fails closed")
        tiny = Budget(UInt64(allowance))
        try:
            _ = _truncated_product(a, b, 0, 0, tiny)
            c.check(False, "truncated product budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "filtered multiplication fails closed")
    c.check(a == geometric(17, 0, False, one, order, work) and b == geometric(17, 1, False, one, order, work), "limited operations preserve borrowed inputs")


def main() raises:
    var c = Checker("incremental kernels")
    var work = Budget()
    domains(ZZ(Integer(1)), c, work)
    domains(QQ(Rational(1)), c, work)
    domains(QQi(ExactComplex(Rational(1), Rational(0))), c, work)
    domains(ModInt(1, Modulus(2, work)), c, work)
    for kind in range(4):
        var order = MonomialOrder() if kind == 0 else MonomialOrder(UInt8(kind), [Int64(0), Int64(1)], work)
        var x = SparsePoly[QQ].variable(0, order, NoContext(), work)
        var y = SparsePoly[QQ].variable(1, order, NoContext(), work)
        var one = SparsePoly[QQ].constant(QQ(Rational(1)), order, work)
        var a = geometric(40, 1, False, QQ(Rational(1)), order, work).mul(x.add(one, work), work)
        var b = geometric(40, 1, True, QQ(Rational(1)), order, work)
        for degree in [0, 1, 7, 39, 80]:
            var got = _truncated_product(a, b, 0, degree, work)
            # (1+y+...)(1-y+...) has even coefficients 1 below degree 40.
            var expected = SparsePoly[QQ].zero(order, NoContext(), work)
            var full = a.mul(b, work)
            for t in full.terms:
                if Int(t.monomial.exponent(1, work)) <= degree:
                    expected.terms.append(t.copy())
            c.check(got == expected, "filtered multiplication agrees with full exact identity")
        var groups = _CoefficientGroups(a, 0, work)
        var eager = _coefficient_groups(a, 0, False, work, first=False)
        for i in range(len(groups.starts) - 1):
            c.check(groups.coefficient(a, i, work) == eager[1][i], "indexed group matches occupied coefficient")
        var accumulator = _PolynomialSum[QQ]()
        for i in range(65):
            accumulator.add(a.neg(work) if i % 2 else a._copy(work), work)
        c.check(accumulator.finish(a, work) == a, "balanced runs preserve cancellation")
        var assignments: List[Tuple[Int64, SparsePoly[QQ]]] = [(Int64(0), y.copy()), (Int64(1), x.copy())]
        var swapped = substitute_many(a, assignments, work)
        c.check(substitute_many(swapped, assignments, work) == a, "simultaneous substitution involution")
        if kind:
            var left: List[SparsePoly[QQ]] = [a.copy(), a.copy(), b.copy()]
            var rows: List[List[SparsePoly[QQ]]] = [[b.copy()], [b.neg(work)], [one.copy()]]
            c.check(_subtract_column(b, left, [0, 1, 2], rows, 0, work).is_zero(), "signed product columns cancel independently")
    for n in [0, 1, 2, 3, 17, 65]:
        var values = List[QQ]()
        for i in range(n):
            values.append(QQ(Rational(i + 1, i + 2)))
        var p = DensePoly[QQ].from_coefficients(0, NoContext(), values^, work)
        c.check(BasisPolynomial[QQ].chebyshev(p, work).to_polynomial(work) == p, "linear-storage Chebyshev round trip")
    c.finish()
