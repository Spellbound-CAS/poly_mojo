"""Sparse algorithm selection, order-sensitive division, and shared budgets."""
from std.testing import assert_equal, assert_true, assert_false
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, PolyError, PolyErrorKind, MonomialOrder
from poly import dense, sparse, exact


def workflow(mut budget: Budget) raises PolyError:
    _ = dense.sign_on_interval([-2, 0, 1], 2, 3, budget)
    var resultant = dense.resultant_certified([-2, 1], [-3, 1], budget)
    _ = exact.check(resultant, budget)
    _ = exact.value(resultant, budget)
    var common = sparse.gcd_checked([(0, -1), (2, 1)], [(0, -1), (1, 1)], budget, method="modular")
    _ = exact.check(common, budget)
    var division = sparse.divide_checked([((1, 0), 1)], [((1, 0), 1), ((0, 2), 1)], budget, order="grlex")
    _ = exact.check(division, budget)
    _ = exact.remainder(division, budget)


def main() raises:
    var h = sparse.asarray([((1, 0), 1), ((0, 1), 1)])
    var a = sparse.multiply(h, [((1, 0), 6), ((0, 0), 6)])
    var b = sparse.multiply(h, [((0, 1), 9), ((0, 0), 9)])
    var expected = sparse.asarray([((1, 0), 3), ((0, 1), 3)])
    var common = sparse.gcd_checked(a, b, method="modular")
    assert_true(exact.check(common))
    assert_equal(exact.value(common).to_list(), expected.to_list())
    assert_equal(sparse.gcd(a, b, method="modular").to_list(), expected.to_list())
    var cofactors = exact.cofactors(common)
    assert_equal(sparse.multiply(exact.value(common), cofactors[0]).to_list(), a.to_list())
    assert_equal(sparse.multiply(exact.value(common), cofactors[1]).to_list(), b.to_list())
    var rational_b = sparse.multiply(h, [((0, 1), Rational(1, 2)), ((0, 0), Rational(1, 2))])
    var rational_gcd = sparse.gcd_checked(a, rational_b, method="modular")
    assert_true(exact.check(rational_gcd))
    assert_equal(exact.value(rational_gcd).to_list(), sparse.asarray([((1, 0), Rational(1)), ((0, 1), Rational(1))]).to_list())
    assert_true(exact.check(sparse.gcd_checked([(0, 0)], [(0, 0)], method="modular")))
    assert_equal(sparse.gcd([(0, 0)], [(0, -2), (1, -2)], method="modular").to_list(), sparse.asarray([(0, 2), (1, 2)]).to_list())
    assert_equal(sparse.gcd([(0, -6)], [(0, 9)], method="modular")[0][1], 3)
    var x = sparse.asarray([((1, 0), 1)])
    var divisor = sparse.asarray([((1, 0), 1), ((0, 2), 1)])
    var negative_y2 = sparse.asarray([((0, 2), Rational(-1))])
    var rational_x = sparse.asarray([((1, 0), Rational(1))])
    var lex = sparse.divmod(x, divisor)
    assert_equal(lex[0].to_list(), sparse.asarray([((0, 0), Rational(1))]).to_list())
    assert_equal(lex[1].to_list(), negative_y2.to_list())
    for order in [String("grlex"), String("grevlex")]:
        var divided = sparse.divmod(x, divisor, order=order)
        assert_equal(len(divided[0]), 0)
        assert_equal(divided[1].to_list(), rational_x.to_list())
        var checked = sparse.divide_checked(x, divisor, order=order)
        assert_true(exact.check(checked))
        assert_equal(exact.remainder(checked).to_list(), rational_x.to_list())
        var many = sparse.divide_many(x, [divisor.copy()], order=order)
        assert_equal(many[1].to_list(), rational_x.to_list())
        var retained = sparse.divide_many_checked(x, [divisor.copy()], order=order)
        assert_true(exact.check(retained))
        assert_equal(exact.remainder(retained).to_list(), rational_x.to_list())
    var reversed = sparse.divide_checked(x, divisor, axes=[1, 0])
    assert_true(exact.check(reversed))
    assert_equal(exact.remainder(reversed).to_list(), rational_x.to_list())
    assert_equal(reversed.evidence().remainder.order.kind, MonomialOrder.LEX)
    var budget = Budget()
    var multi = sparse.divide_many_checked(x, [divisor.copy(), x.copy()], budget, axes=[1, 0])
    assert_true(exact.check(multi))
    assert_equal(len(exact.remainder(multi)), 0)
    var ordinary = sparse.divide_many(x, [divisor.copy(), x.copy()], budget, axes=[1, 0])
    assert_equal(ordinary[0][1].to_list(), sparse.asarray([((0, 0), Rational(1))]).to_list())
    # Output arrays retain canonical term order and remain usable by other operations.
    var canonical = sparse.divmod(sparse.add(x, [((0, 1), 1)]), divisor, budget, axes=[1, 0])[1].copy()
    assert_equal(canonical.to_list(), sparse.asarray([((0, 1), Rational(1)), ((1, 0), Rational(1))]).to_list())
    assert_equal(sparse.add(canonical, canonical).to_list(), sparse.asarray([((0, 1), Rational(2)), ((1, 0), Rational(2))]).to_list())
    var empty = sparse.divide_many_checked(x, List[type_of(x)](), order="grevlex", axes=[1, 0])
    assert_true(exact.check(empty))
    assert_equal(exact.remainder(empty).to_list(), rational_x.to_list())
    var finite_x = sparse.asarray([((1, 0), 1)], modulus=5)
    var finite_d = sparse.asarray([((1, 0), 1), ((0, 2), 1)], modulus=5)
    assert_true(exact.check(sparse.divide_checked(finite_x, finite_d, order="grlex")))
    var complex_x = sparse.asarray([((1, 0), ExactComplex(0, 1))])
    assert_true(exact.check(sparse.divide_checked(complex_x, complex_x, axes=[1, 0])))
    for scenario in range(11):
        var failed = False
        try:
            if scenario == 0:
                _ = sparse.divmod(x, divisor, order="bad")
            elif scenario == 1:
                _ = sparse.divide_many(x, [divisor.copy()], axes=[0, 0])
            elif scenario == 2:
                _ = sparse.divide_checked(x, divisor, axes=[0])
            elif scenario == 3:
                _ = sparse.divide_many_checked(x, [divisor.copy()], axes=[0, 2])
            elif scenario == 4:
                _ = sparse.divmod(x, divisor, axes=[-1, 0])
            elif scenario == 5:
                _ = sparse.gcd(a, b, method="bad")
            elif scenario == 6:
                _ = sparse.divide_many(x, List[type_of(x)](), order="bad")
            elif scenario == 7:
                _ = sparse.divide_checked(x, sparse.subtract(x, x), order="grlex")
            elif scenario == 8:
                _ = sparse.gcd(finite_x, finite_d, method="modular")
            elif scenario == 9:
                _ = sparse.gcd_checked(complex_x, complex_x, method="modular")
            else:
                _ = sparse.divide_many_checked(x, [divisor.copy()], max_steps=0)
        except e:
            failed = e.kind == (PolyErrorKind.LIMIT if scenario == 10 else PolyErrorKind.UNSUPPORTED if scenario >= 8 else PolyErrorKind.DOMAIN)
        assert_true(failed)
    # Warm deterministic prime caches before measuring exact shared-budget costs.
    var warm = Budget()
    workflow(warm)
    var meter = Budget()
    workflow(meter)
    var cost = UInt64.MAX - meter.remaining()
    var enough = Budget(cost)
    workflow(enough)
    assert_equal(enough.remaining(), UInt64(0))
    var short = Budget(cost - 1)
    var limited = False
    try:
        workflow(short)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(a.to_list(), sparse.multiply(h, [((1, 0), 6), ((0, 0), 6)]).to_list())
    common._data[].proof.cofactors[0].terms.clear()
    assert_false(exact.check(common))
    print("array sparse policy checks passed")
