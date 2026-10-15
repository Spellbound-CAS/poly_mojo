"""Algebraic identities, promotion, iteration, and resource limits through poly."""
import poly
from apn_mojo import Integer, Rational, ExactComplex
from polynomials import Budget, PolyError, PolyErrorKind, ModInt
from std.testing import assert_equal, assert_true


def test_dense() raises:
    var p = poly.dense.asarray([1, 2, 3])
    assert_equal(poly.dense.subtract(p, [1, 1]).to_list(), poly.dense.asarray([0, 1, 3]).to_list())
    assert_equal(poly.dense.add(p, poly.dense.negative(p)).to_list(), List[Integer]())
    assert_equal(poly.dense.power([1, 1], 3).to_list(), poly.dense.asarray([1, 3, 3, 1]).to_list())
    assert_equal(poly.dense.power(List[Int](), 0).to_list(), poly.dense.asarray([1]).to_list())
    assert_equal(poly.dense.power(List[Int](), 4).to_list(), List[Integer]())
    assert_equal(poly.dense.degree([1, 0, 0]), 0)
    assert_equal(poly.dense.degree(List[Int]()), -1)
    var integral = poly.dense.integral([1, 1])
    assert_equal(integral.domain(), "QQ")
    assert_equal(integral.to_list(), poly.coefficients((0, 1, Rational(1, 2))).to_list())
    assert_equal(poly.dense.derivative(integral).to_list(), poly.coefficients((Rational(1), 1)).to_list())
    var division = poly.dense.divmod([1, 0, 1], [1, 2])
    assert_equal(division[0].to_list(), poly.coefficients((Rational(-1, 4), Rational(1, 2))).to_list())
    assert_equal(division[1].to_list(), poly.coefficients((Rational(5, 4),)).to_list())
    assert_equal(poly.dense.add(poly.dense.multiply(division[0], [1, 2]), division[1]).to_list(), poly.coefficients((Rational(1), 0, 1)).to_list())
    assert_equal(poly.dense.compose([1, 2, 3], [1, 1]).to_list(), poly.dense.asarray([6, 8, 3]).to_list())
    assert_equal(poly.dense.from_roots([0, 1, 1]).to_list(), poly.dense.asarray([0, 1, -2, 1]).to_list())
    assert_equal(poly.dense.from_roots(List[Rational]()).to_list(), poly.coefficients((Rational(1),)).to_list())
    var output = poly.dense.evaluate_many([0, 1], [2, 0, 0])
    assert_equal(len(output), 3)
    assert_equal(output.to_list(), poly.coefficients((2, 0, 0)).to_list())
    assert_equal(poly.dense.evaluate_many(p, List[Int]()).to_list(), List[Integer]())
    var points = poly.coefficients((1, Rational(1, 2), 0))
    var evaluated = poly.dense.evaluate_many(p, points)
    assert_equal(evaluated.domain(), "QQ")
    assert_equal(evaluated.to_list(), poly.coefficients((Rational(6), Rational(11, 4), 1)).to_list())
    assert_equal(p.to_list(), poly.dense.asarray([1, 2, 3]).to_list())


def test_sparse() raises:
    var p = poly.sparse.asarray([((1, 0), 2), ((0, 2), 3)])
    var q = poly.sparse.asarray([((1, 0), 1), ((0, 0), 1)])
    var difference = poly.sparse.subtract(p, q)
    assert_equal(poly.sparse.evaluate(difference, (2, 3)), Integer(28))
    assert_equal(len(poly.sparse.add(p, poly.sparse.negative(p))), 0)
    var squared = poly.sparse.power(p, 2)
    assert_equal(poly.sparse.evaluate(squared, (2, 3)), Integer(961))
    assert_equal(poly.sparse.degree(squared), 4)
    assert_equal(poly.sparse.degree(squared, axis=0), 2)
    assert_equal(poly.sparse.degree(squared, axis=1), 4)
    var zero = poly.sparse.subtract(p, p)
    assert_equal(poly.sparse.degree(zero), -1)
    assert_equal(poly.sparse.evaluate(poly.sparse.power(zero, 0), (0, 0)), Integer(1))
    var integrated = poly.sparse.integral(p, axis=1)
    assert_equal(integrated.domain(), "QQ")
    assert_equal(poly.sparse.evaluate(poly.sparse.derivative(integrated, axis=1), (2, 3)), Rational(31))
    # x^2*y + y divided by x*y+1 yields x and y-x under lex order.
    var dividend = [((2, 1), 1), ((0, 1), 1)]
    var divisor = [((1, 1), 1), ((0, 0), 1)]
    var division = poly.sparse.divmod(dividend, divisor)
    assert_equal(poly.sparse.evaluate(division[0], (2, 3)), Rational(2))
    assert_equal(poly.sparse.evaluate(division[1], (2, 3)), Rational(1))
    var restored = poly.sparse.add(poly.sparse.multiply(division[0], divisor), division[1])
    assert_equal(poly.sparse.evaluate(restored, (2, 3)), Rational(15))
    assert_equal(poly.sparse.evaluate_many(p, [(0, 0), (2, 3), (0, 0)]).to_list(), poly.coefficients((0, 31, 0)).to_list())
    assert_equal(poly.sparse.evaluate_many([(1, 1)], [0, 2, 0]).to_list(), poly.coefficients((0, 2, 0)).to_list())
    var rational_points = [(Rational(1, 2), Rational(1, 2))]
    assert_equal(poly.sparse.evaluate_many(p, rational_points).to_list(), poly.coefficients((Rational(7, 4),)).to_list())


def test_complex_and_modular() raises:
    var z = poly.coefficients((1, ExactComplex(0, 1)))
    var division = poly.dense.divmod(poly.dense.power(z, 2), z)
    assert_equal(division[0].domain(), "QQi")
    assert_equal(division[0].to_list(), z.to_list())
    assert_equal(len(division[1]), 0)
    assert_equal(poly.dense.derivative(poly.dense.integral(z)).to_list(), z.to_list())
    var sparse_z = poly.to_sparse(z)
    assert_equal(poly.sparse.derivative(poly.sparse.integral(sparse_z)).to_list(), sparse_z.to_list())
    var sparse_division = poly.sparse.divmod(poly.sparse.power(sparse_z, 2), sparse_z)
    assert_equal(sparse_division[0].to_list(), sparse_z.to_list())
    var field = poly.dense.asarray([1, 2], modulus=7)
    var mod_division = poly.dense.divmod(poly.dense.power(field, 2), field)
    assert_equal(mod_division[0].to_list(), field.to_list())
    assert_equal(mod_division[1].modulus(), UInt64(7))
    assert_equal(poly.dense.derivative(poly.dense.integral(field)).to_list(), field.to_list())
    var mod_sparse = poly.to_sparse(field)
    var mod_sparse_division = poly.sparse.divmod(poly.sparse.power(mod_sparse, 2), mod_sparse)
    assert_equal(mod_sparse_division[0].to_list(), mod_sparse.to_list())
    assert_equal(mod_sparse_division[1].modulus(), UInt64(7))
    assert_equal(poly.sparse.derivative(poly.sparse.integral(mod_sparse)).to_list(), mod_sparse.to_list())
    var modular_points = poly.coefficients((poly.modint(2, 7), 0, 0))
    var results = poly.dense.evaluate_many([0, 1], modular_points)
    assert_equal(results.modulus(), UInt64(7))
    assert_equal(results.to_list(), modular_points.to_list())
    assert_equal(poly.sparse.evaluate_many([(1, 1)], modular_points).to_list(), modular_points.to_list())
    var roots = poly.dense.from_roots(modular_points)
    assert_equal(poly.dense.evaluate(roots, 2).value, UInt64(0))
    assert_equal(poly.dense.evaluate(roots, 0).value, UInt64(0))
    var empty_field = poly.dense.asarray(List[Int](), modulus=7)
    assert_equal(poly.dense.from_roots(empty_field).modulus(), UInt64(7))
    assert_equal(poly.dense.evaluate_many(empty_field, List[Int]()).modulus(), UInt64(7))


def test_conversions() raises:
    var input = poly.dense.asarray([3, 0, 0, -4, 0])
    var sparse = poly.to_sparse(input)
    assert_equal(sparse.to_list(), poly.sparse.asarray([(0, 3), (3, -4)]).to_list())
    assert_equal(poly.to_dense(sparse).to_list(), input.to_list())
    assert_equal(poly.to_dense([((2,), 1)]).to_list(), poly.dense.asarray([0, 0, 1]).to_list())
    var zero = poly.dense.asarray(List[Int](), modulus=7)
    assert_equal(poly.to_dense(poly.to_sparse(zero)).modulus(), UInt64(7))
    # Construct the endpoint at runtime so literal narrowing cannot turn a
    # capacity-limit check into a negative-power rejection on Mojo 1.1.
    var maximum_power = Int("4294967295")
    assert_equal(UInt64(maximum_power), UInt64(4294967295))
    var failed = False
    try:
        _ = poly.to_dense([(maximum_power, 1)], max_steps=10000)
    except error:
        assert_equal(error.kind, PolyErrorKind.LIMIT)
        failed = True
    assert_true(failed)


def workflow(mut budget: Budget) raises PolyError:
    var p = poly.dense.power([1, 1], 4, budget)
    var division = poly.dense.divmod(p, [2, 1], budget)
    var sparse = poly.to_sparse(division[0], budget)
    var integral = poly.sparse.integral(sparse, budget)
    var restored = poly.to_dense(poly.sparse.derivative(integral, budget), budget)
    _ = poly.dense.evaluate_many(restored, [0, 1, 2, 0], budget)
    _ = poly.sparse.evaluate_many(sparse, [0, 1, 2, 0], budget)
    _ = poly.dense.from_roots([0, 1, 1], budget)


def test_budget_and_errors() raises:
    var measured = Budget()
    workflow(measured)
    var used = UInt64.MAX - measured.remaining()
    var exact = Budget(used)
    workflow(exact)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(used - 1)
    var failed = False
    try:
        workflow(short)
    except error:
        assert_equal(error.kind, PolyErrorKind.LIMIT)
        failed = True
    assert_true(failed)
    # Each API reports typed domain errors without changing prepared inputs.
    var p = poly.dense.asarray([1, 2])
    for scenario in range(7):
        var caught = False
        try:
            if scenario == 0:
                _ = poly.dense.power(p, -1)
            elif scenario == 1:
                _ = poly.dense.divmod(p, List[Int]())
            elif scenario == 2:
                _ = poly.sparse.divmod([(1, 1)], List[Tuple[Int, Int]]())
            elif scenario == 3:
                _ = poly.sparse.integral([((0, 1), 1)], axis=2)
            elif scenario == 4:
                _ = poly.to_dense([((0, 0), 1)])
            elif scenario == 5:
                _ = poly.dense.integral(poly.dense.asarray([0, 1], modulus=2))
            else:
                _ = poly.sparse.integral(poly.sparse.asarray([(1, 1)], modulus=2))
        except error:
            assert_equal(error.kind, PolyErrorKind.DOMAIN)
            caught = True
        assert_true(caught)
    assert_equal(p.to_list(), poly.dense.asarray([1, 2]).to_list())


def main() raises:
    test_dense()
    test_sparse()
    test_complex_and_modular()
    test_conversions()
    test_budget_and_errors()
    print("array operation checks passed")
