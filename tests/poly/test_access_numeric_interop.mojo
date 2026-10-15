"""Numerical variable maps preserve exact stored formats, bounds, and ownership."""
from poly import dense, sparse
from apn_mojo import Float, Ball, ComplexBall, Rational, ArithmeticContext, FloatFormat, BallContext
from polynomials import Budget, PolyErrorKind, NumericSparsePoly, FloatDomain, MonomialOrder
from std.testing import assert_equal, assert_true


def main() raises:
    var fc = ArithmeticContext(format=FloatFormat(precision=37))
    var p = dense.asarray([Rational(1, 3), Rational(2)], context=fc)
    var core = dense.to_core(p, variable=-7)
    assert_equal(core.variable, Int64(-7))
    assert_equal(core.coefficients[0].to_json(), p[0].to_json())
    var restored = dense.from_core(core)
    assert_equal(dense.to_core(restored).variable, Int64(0))
    core.coefficients[0] = Float(9, context=fc)
    assert_equal(restored[0].to_json(), p[0].to_json())
    var sp = sparse.asarray([((2, 0, 0), 3), ((0, 1, 0), 5), ((1, 1, 0), 7)], context=fc)
    var ids: List[Int64] = [19, -4, 200]
    var mapped = sparse.to_core(sp, variables=ids)
    assert_equal(mapped.order.variables, ids)
    var work = Budget()
    var grlex = MonomialOrder(MonomialOrder.GRLEX, List[Int64]([Int64(200), Int64(-4), Int64(19)]), work)
    var reordered = NumericSparsePoly[FloatDomain].from_terms(mapped.terms.copy(), grlex, fc, work)
    var imported = sparse.from_core[ndim=3](reordered, variables=ids)
    assert_equal(imported.to_list(), sp.to_list())
    assert_equal(sparse.to_core(imported).order.kind, MonomialOrder.LEX)
    var kernel = NumericSparsePoly[FloatDomain].from_terms(mapped.terms.copy(), MonomialOrder(), fc, work)
    assert_equal(sparse.from_core[ndim=3](kernel, variables=ids).to_list(), sp.to_list())
    var permuted = sparse.from_core[ndim=3](mapped, variables=List[Int64]([Int64(-4), Int64(19), Int64(200)]))
    assert_equal(sparse.evaluate(permuted, (2, 3, 4), context=fc).to_rational_exact(), Rational(79))
    var bc = BallContext(precision=53)
    var uncertain = sparse.asarray([((1, 0), ComplexBall(Ball(Rational(1), Rational(1, 8), precision=53), Ball(0, precision=29)))], context=bc)
    var balls = sparse.to_core(uncertain, variables=List[Int64]([Int64(5), Int64(8)]))
    var again = sparse.from_core[ndim=2](balls, variables=List[Int64]([Int64(5), Int64(8)]))
    assert_equal(again[0][1].to_json(), uncertain[0][1].to_json())
    var empty = sparse.asarray(List[Tuple[Tuple[Int, Int], Int]](), context=fc)
    var empty_core = sparse.to_core(empty, variables=List[Int64]([Int64(7), Int64(9)]))
    assert_equal(len(sparse.from_core[ndim=2](empty_core, variables=List[Int64]([Int64(7), Int64(9)]))), 0)
    var scalar_zero = Float(0, context=fc)
    for mode in range(8):
        var rejected = False
        try:
            if mode == 0: _ = sparse.from_core[ndim=3](mapped)
            elif mode == 1: _ = sparse.to_core(sp, variables=List[Int64]([Int64(1), Int64(1), Int64(3)]))
            elif mode == 2: _ = sparse.to_core(sp, variables=List[Int64]([Int64(1)]))
            elif mode == 3:
                var bad = sparse.to_core(sp)
                bad.terms.swap_elements(0, 1)
                _ = sparse.from_core[ndim=3](bad)
            elif mode == 4:
                var bad = sparse.to_core(empty)
                bad.order.kind = 255
                _ = sparse.from_core[ndim=2](bad)
            elif mode == 5:
                var bad = dense.to_core(p)
                bad.coefficients.append(scalar_zero.copy())
                _ = dense.from_core(bad)
            elif mode == 6:
                var bad = sparse.to_core(sp)
                bad.terms[0].monomial.pairs[0] = (Int64(1), UInt32(0))
                _ = sparse.from_core[ndim=3](bad)
            else:
                var bad = sparse.to_core(sp)
                bad.terms[1] = bad.terms[0].copy()
                _ = sparse.from_core[ndim=3](bad)
        except e:
            rejected = e.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    # Exercise reordering at every small heap size, including single-child nodes.
    for n in range(33):
        var terms = List[Tuple[Tuple[Int, Int], Int]]()
        var expected = List[Tuple[Tuple[Int, Int], Int]]()
        for i in range(n):
            terms.append(((i, n - i), i + 1))
            expected.append(((n - i, i), i + 1))
        var input = sparse.asarray(terms, context=fc)
        var exported = sparse.to_core(input, variables=List[Int64]([Int64(-9), Int64(22)]))
        var swapped = sparse.from_core[ndim=2](exported, variables=List[Int64]([Int64(22), Int64(-9)]))
        assert_equal(swapped.to_list(), sparse.asarray(expected, context=fc).to_list())
    var meter = Budget()
    _ = sparse.from_core[ndim=3](mapped, meter, variables=ids)
    var used = UInt64.MAX - meter.remaining()
    var exact = Budget(used)
    _ = sparse.from_core[ndim=3](mapped, exact, variables=ids)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(used - 1)
    var limited = False
    try: _ = sparse.from_core[ndim=3](mapped, short, variables=ids)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(sparse.from_core[ndim=3](mapped, variables=ids).to_list(), sp.to_list())
    print("access numerical interoperability checks passed")
