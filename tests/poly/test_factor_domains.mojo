"""Prime-field inseparability, Gaussian splitting, and exact factor evidence."""
import poly
from apn_mojo import ExactComplex, Rational
from polynomials import Budget, PolyErrorKind, ModInt, Modulus, DensePoly, factor_finite_field, factor_gaussian, check_factorization, check_squarefree_product, squarefree, divmod, derivative_dense, extended_gcd, NoContext, QQi
from std.testing import assert_equal, assert_true


def finite_checks() raises:
    for prime in [2, 3, 5]:
        var a = poly.dense.asarray([1, 1], modulus=UInt64(prime))
        var b = poly.dense.asarray([0, 1], modulus=UInt64(prime))
        var f = poly.dense.multiply(poly.dense.power(a, prime * prime), poly.dense.power(b, prime + 1))
        var checked = poly.dense.factor_checked(f)
        assert_true(poly.exact.check(checked))
        var factors = poly.exact.factors(checked)
        assert_equal(len(factors), 2)
        assert_equal(factors[0][0].to_list(), b.to_list())
        assert_equal(factors[0][1], prime + 1)
        assert_equal(factors[1][0].to_list(), a.to_list())
        assert_equal(factors[1][1], prime * prime)
        var blocks = poly.dense.squarefree_checked(f)
        assert_true(poly.exact.check(blocks))
        var grouped = poly.exact.factors(blocks)
        assert_equal(grouped[0][1], prime + 1)
        assert_equal(grouped[1][1], prime * prime)
        var repeated = poly.dense.factor(f, seed=77)
        assert_equal(repeated[1][0][0].to_list(), factors[0][0].to_list())
        assert_equal(repeated[1][1][0].to_list(), factors[1][0].to_list())
        for c in [0, 1, prime - 1]:
            var constant = poly.dense.asarray([c], modulus=UInt64(prime))
            var factorized = poly.dense.factor(constant)
            assert_equal(factorized[0].value, UInt64(c))
            assert_equal(len(factorized[1]), 0)
            assert_true(poly.exact.check(poly.dense.squarefree_checked(constant)))
    # Derivative-free residual plus an ordinary repeated irreducible block.
    var linear = poly.dense.asarray([1, 1], modulus=2)
    var quadratic = poly.dense.asarray([1, 1, 1], modulus=2)
    var f = poly.dense.multiply(poly.dense.power(linear, 6), poly.dense.power(quadratic, 3))
    var result = poly.dense.factor(f)
    assert_equal(len(result[1]), 2)
    assert_equal(result[1][0][1], 6)
    assert_equal(result[1][1][0].to_list(), quadratic.to_list())
    assert_equal(result[1][1][1], 3)
    var grouped = poly.dense.squarefree(f)
    assert_equal(grouped[1][0][1], 3)
    assert_equal(grouped[1][1][1], 6)
    var large_prime = poly.dense.asarray([1, 2, 1], modulus=2147483647)
    assert_equal(poly.dense.factor(large_prime)[1][0][1], 2)


def gaussian_checks() raises:
    # x^2+1 splits over QQ(i), while x^2+2 does not.
    var f = poly.dense.asarray([ExactComplex(1), ExactComplex(0), ExactComplex(1)])
    var result = poly.dense.factor(f)
    assert_equal(result[0], ExactComplex(1))
    assert_equal(len(result[1]), 2)
    assert_equal(result[1][0][0].to_list(), [ExactComplex(0, -1), ExactComplex(1)])
    assert_equal(result[1][1][0].to_list(), [ExactComplex(0, 1), ExactComplex(1)])
    var irreducible = poly.dense.asarray([ExactComplex(2), ExactComplex(0), ExactComplex(1)])
    var untouched = poly.dense.factor(irreducible)
    assert_equal(len(untouched[1]), 1)
    assert_equal(untouched[1][0][0].to_list(), irreducible.to_list())
    # Nonreal coefficients, a nonreal unit, and a denominator in a root.
    var a = poly.dense.asarray([ExactComplex(Rational(1, 2), Rational(-1)), ExactComplex(1)])
    var b = poly.dense.asarray([ExactComplex(0, 2), ExactComplex(1)])
    var product = poly.dense.multiply(poly.dense.power(a, 2), b)
    product = poly.dense.multiply(product, [ExactComplex(2, -3)])
    var proof = poly.dense.factor_checked(product)
    assert_true(poly.exact.check(proof))
    assert_equal(poly.exact.unit(proof), ExactComplex(2, -3))
    var factors = poly.exact.factors(proof)
    assert_equal(len(factors), 2)
    assert_equal(factors[0][0].to_list(), b.to_list())
    assert_equal(factors[0][1], 1)
    assert_equal(factors[1][0].to_list(), a.to_list())
    assert_equal(factors[1][1], 2)
    for c in [ExactComplex(0), ExactComplex(2, -3)]:
        var constant = poly.dense.factor([c.copy()])
        assert_equal(constant[0], c)
        assert_equal(len(constant[1]), 0)
    var core = poly.dense.to_core(product)
    var budget = Budget()
    var bad = factor_gaussian(core, budget)
    bad.factors[0][1] += 1
    assert_true(not check_factorization(core, bad, budget))


def exhaustive_finite() raises:
    # Every small monic input: verify products, squarefree blocks, and
    # irreducibility independently by testing all possible monic divisors.
    var budget = Budget()
    for p in [2, 3, 5]:
        var ctx = Modulus(UInt64(p), budget)
        var limit = 4 if p == 2 else (3 if p == 3 else 2)
        var combinations = 1
        for degree in range(1, limit + 1):
            combinations *= p
            for code in range(combinations):
                var digits = code
                var values = List[ModInt]()
                for _ in range(degree):
                    values.append(ModInt(UInt64(digits % p), ctx))
                    digits //= p
                values.append(ModInt.one(ctx))
                var f = DensePoly[ModInt].from_coefficients(0, ctx, values^, budget)
                var decomposition = factor_finite_field(f, budget)
                assert_true(check_factorization(f, decomposition, budget))
                var blocks = squarefree(f, budget)
                assert_true(check_squarefree_product(f, blocks, budget))
                for block in blocks.factors:
                    var derivative = derivative_dense(block[0], budget)
                    var gcd = extended_gcd(block[0], derivative, budget)
                    assert_equal(gcd.gcd.total_degree(budget), 0)
                for i in range(len(blocks.factors)):
                    for j in range(i):
                        var gcd = extended_gcd(blocks.factors[i][0], blocks.factors[j][0], budget)
                        assert_equal(gcd.gcd.total_degree(budget), 0)
                for entry in decomposition.factors:
                    var possibilities = 1
                    for d in range(1, entry[0].degree() // 2 + 1):
                        possibilities *= p
                        for n in range(possibilities):
                            var lower = n
                            var trial = List[ModInt]()
                            for _ in range(d):
                                trial.append(ModInt(UInt64(lower % p), ctx))
                                lower //= p
                            trial.append(ModInt.one(ctx))
                            var divisor = DensePoly[ModInt].from_coefficients(0, ctx, trial^, budget)
                            assert_true(not divmod(entry[0], divisor, budget).remainder.is_zero())
                var bad = decomposition.copy()
                bad.factors[0][1] += 1
                assert_true(not check_factorization(f, bad, budget))


def budgets() raises:
    var finite = poly.dense.asarray([1, 0, 1], modulus=2)
    var gaussian = poly.dense.asarray([ExactComplex(1), ExactComplex(0), ExactComplex(1)])
    for steps in [UInt64(0), UInt64(20), UInt64(100)]:
        for operation in range(4):
            var failed = False
            try:
                if operation == 0:
                    _ = poly.dense.factor(finite, max_steps=steps)
                elif operation == 1:
                    _ = poly.dense.factor_checked(gaussian, max_steps=steps)
                elif operation == 2:
                    _ = poly.dense.squarefree(finite, max_steps=steps)
                else:
                    _ = poly.dense.squarefree_checked(finite, max_steps=steps)
            except e:
                failed = e.kind == PolyErrorKind.LIMIT
            assert_true(failed)
    assert_equal(gaussian.to_list(), [ExactComplex(1), ExactComplex(0), ExactComplex(1)])
    assert_equal(poly.dense.degree(finite), 2)


def main() raises:
    finite_checks()
    gaussian_checks()
    exhaustive_finite()
    budgets()
    print("array factor domain checks passed")
