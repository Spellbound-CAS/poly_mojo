"""Domain boundaries, zero data, shared work limits, and borrowing after failure."""
import poly
from apn_mojo import Integer, Rational, ExactComplex, Batch
from polynomials import Budget, PolyError, PolyErrorKind, ModInt
from std.testing import assert_equal, assert_true


def workflow(p: poly.dense.Array[0], mut budget: Budget) raises PolyError:
    var core = poly.dense.to_core(p, budget, variable=12)
    var copied = poly.dense.from_core(core, budget)
    var common = poly.dense.gcd(copied, [-1, 1], budget)
    var blocks = poly.dense.squarefree(copied, budget)
    var plan = poly.dense.prepare_points([0, 1, 2], budget, interpolation=True)
    var values = plan.evaluate(common, budget)
    _ = plan.interpolate(values, budget)
    _ = blocks[1][0][0].to_list(budget)


def main() raises:
    var p = poly.dense.asarray([1, -2, 1])
    var meter = Budget()
    workflow(p, meter)
    var cost = UInt64.MAX - meter.remaining()
    assert_true(cost > 0)
    var exact = Budget(cost)
    workflow(p, exact)
    assert_equal(exact.remaining(), UInt64(0))
    var failed = False
    var short = Budget(cost - 1)
    try:
        workflow(p, short)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    assert_equal(p.to_list(), poly.dense.asarray([1, -2, 1]).to_list())
    assert_equal(len(poly.dense.gcd(List[Int](), List[Int]())), 0)
    assert_equal(poly.dense.content(List[Int]()), Integer(0))
    assert_equal(len(poly.dense.primitive_part(List[Int]())), 0)
    var zero_factors = poly.dense.factor(List[Int]())
    assert_equal(zero_factors[0], Integer(0))
    assert_equal(len(zero_factors[1]), 0)
    assert_equal(poly.dense.factor([-7])[0], Integer(-7))
    assert_equal(len(poly.dense.interpolate(List[Int](), List[Int]())), 0)
    var finite = poly.dense.asarray([-1, 0, 1], modulus=5)
    var common = poly.dense.gcd(finite, [-1, 1])
    assert_equal(common.to_list(), poly.dense.asarray([-1, 1], modulus=5).to_list())
    assert_equal(poly.dense.resultant(finite, [-2, 1]).value, UInt64(3))
    var sf = poly.sparse.factor(poly.to_sparse(finite))
    assert_equal(len(sf[1]), 2)
    var complex_p = poly.coefficients((ExactComplex(0, 1), 1))
    var complex_sq = poly.dense.squarefree(poly.dense.power(complex_p, 2))
    assert_equal(complex_sq[1][0][1], 2)
    var complex_bezout = poly.dense.extended_gcd(complex_p, [1])
    assert_equal(complex_bezout[0].domain(), "QQi")
    var empty_mod = poly.dense.subtract(finite, finite)
    var empty_plan = poly.dense.prepare_points(empty_mod)
    assert_equal(empty_plan.domain(), "GF(5)")
    assert_equal(len(empty_plan.interpolate(List[ModInt]())), 0)
    assert_equal(len(poly.dense.interpolate(empty_mod, List[ModInt]())), 0)
    var other_mod = poly.dense.asarray([0, 1], modulus=7)
    var plan = poly.dense.prepare_points(poly.coefficients((poly.modint(0, 5), 1)))
    var grid = Batch[Integer]([1, 2, 3, 4]).reshape([2, 2])
    for scenario in range(11):
        failed = False
        try:
            if scenario == 0:
                _ = poly.dense.resultant([1], [2])
            elif scenario == 1:
                _ = poly.dense.pseudo_divmod(p, List[Int]())
            elif scenario == 2:
                _ = poly.dense.interpolate([0, 1], [1])
            elif scenario == 3:
                _ = poly.dense.hermite_interpolate([0], [List[Int]()])
            elif scenario == 4:
                _ = poly.dense.prepare_points([0, 0], interpolation=True)
            elif scenario == 5:
                _ = poly.dense.jet_evaluate(p, 0, -1)
            elif scenario == 6:
                _ = plan.evaluate(other_mod)
            elif scenario == 7:
                _ = plan.interpolate(other_mod)
            elif scenario == 8:
                _ = poly.dense.prepare_points(grid)
            elif scenario == 9:
                _ = poly.sparse.extended_gcd([((1, 0), 1)], [((0, 0), 1)])
            elif scenario == 10:
                _ = poly.sparse.content([((1, 0), 1)], axis=2)
        except e:
            failed = e.kind == PolyErrorKind.DOMAIN
        assert_true(failed)
    print("array exact boundary checks passed")
