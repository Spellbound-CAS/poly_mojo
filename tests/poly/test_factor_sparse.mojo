"""Multivariable field factors, sparse multiplicity blocks and resource boundaries."""
import poly
from apn_mojo import ExactComplex, Rational
from polynomials import Budget, PolyErrorKind, check_multivariate_factorization, check_sparse_squarefree_product, squarefree_multivariate, factor_multivariate_finite_field, factor_multivariate_gaussian, MonomialOrder
from std.testing import assert_equal, assert_true


def finite_checks() raises:
    for p in [2, 3, 5]:
        var a = poly.sparse.asarray([((1, 0), 1), ((0, 1), 1)], modulus=UInt64(p))
        var b = poly.sparse.asarray([((0, 1), 1), ((0, 0), 1)], modulus=UInt64(p))
        var f = poly.sparse.multiply(poly.sparse.power(a, p), b)
        var proof = poly.sparse.factor_checked(f)
        assert_true(poly.exact.check(proof))
        var factors = poly.exact.factors(proof)
        assert_equal(len(factors), 2)
        assert_equal(factors[0][0].to_list(), b.to_list())
        assert_equal(factors[0][1], 1)
        assert_equal(factors[1][0].to_list(), a.to_list())
        assert_equal(factors[1][1], p)
        var grouped = poly.sparse.squarefree_checked(f)
        assert_true(poly.exact.check(grouped))
        var blocks = poly.exact.factors(grouped)
        assert_equal(blocks[0][1], 1)
        assert_equal(blocks[1][1], p)
        assert_equal(blocks[1][0].to_list(), a.to_list())
        var raw = poly.sparse.factor(f, seed=3)
        assert_equal(raw[1][0][0].to_list(), b.to_list())
        assert_equal(raw[1][1][0].to_list(), a.to_list())
        var empty = poly.sparse.asarray(List[Tuple[Tuple[Int, Int], Int]](), modulus=UInt64(p))
        assert_equal(poly.sparse.factor(empty)[0].value, UInt64(0))
        assert_true(poly.exact.check(poly.sparse.squarefree_checked(empty)))
    # Repeated inverse Frobenius stages keep a high sparse multiplicity small.
    var high = poly.sparse.asarray([((4096, 0), 1), ((0, 4096), 1)], modulus=2)
    var high_result = poly.sparse.factor_checked(high, max_steps=1000000)
    var high_factors = poly.exact.factors(high_result)
    assert_equal(len(high_factors), 1)
    assert_equal(high_factors[0][1], 4096)
    assert_equal(poly.sparse.degree(high_factors[0][0]), 1)
    assert_true(poly.exact.check(high_result))
    # This irreducible maps to t^2+t, a split univariate image in GF(2).
    var irreducible = poly.sparse.asarray([((1, 0), 1), ((0, 1), 1)], modulus=2)
    var result = poly.sparse.factor(irreducible)
    assert_equal(len(result[1]), 1)
    assert_equal(result[1][0][0].to_list(), irreducible.to_list())
    # Several factors can map to the same univariate factor; retain multiplicity.
    var monomial = poly.sparse.asarray([((2, 3), 2)], modulus=5)
    var monomial_factors = poly.sparse.factor(monomial)
    assert_equal(monomial_factors[0].value, UInt64(2))
    assert_equal(len(monomial_factors[1]), 2)
    assert_equal(monomial_factors[1][0][1], 3)
    assert_equal(monomial_factors[1][1][1], 2)
    var constant = poly.sparse.asarray([((0, 0), 3)], modulus=5)
    assert_equal(poly.sparse.squarefree(constant)[0].value, UInt64(3))
    assert_equal(len(poly.sparse.factor(constant)[1]), 0)


def complex_checks() raises:
    var a = poly.sparse.asarray([((1, 0), ExactComplex(1)), ((0, 1), ExactComplex(0, 1))])
    var b = poly.sparse.asarray([((1, 0), ExactComplex(1)), ((0, 1), ExactComplex(0, -1))])
    var f = poly.sparse.multiply(a, b)
    var proof = poly.sparse.factor_checked(f)
    assert_true(poly.exact.check(proof))
    var factors = poly.exact.factors(proof)
    assert_equal(len(factors), 2)
    assert_equal(factors[0][0].to_list(), b.to_list())
    assert_equal(factors[1][0].to_list(), a.to_list())
    var repeated = poly.sparse.multiply(poly.sparse.power(a, 2), b)
    assert_true(poly.exact.check(poly.sparse.factor_checked(repeated)))
    var blocks = poly.sparse.squarefree(repeated)
    assert_equal(blocks[1][0][0].to_list(), b.to_list())
    assert_equal(blocks[1][0][1], 1)
    assert_equal(blocks[1][1][0].to_list(), a.to_list())
    assert_equal(blocks[1][1][1], 2)
    # The one-axis evidence format and Int / tuple power shapes remain usable.
    var uni = poly.sparse.asarray([(0, ExactComplex(1)), (2, ExactComplex(1))])
    assert_equal(len(poly.sparse.factor(uni)[1]), 2)
    assert_true(poly.exact.check(poly.sparse.squarefree_checked(uni)))


def rational_blocks() raises:
    var x = poly.sparse.asarray([((1, 0, 0), 1)])
    var y = poly.sparse.asarray([((0, 1, 0), 1), ((0, 0, 0), 1)])
    var z = poly.sparse.asarray([((0, 0, 1), 1), ((0, 0, 0), -1)])
    var f = poly.sparse.multiply(poly.sparse.power(x, 2), poly.sparse.power(poly.sparse.multiply(y, z), 3))
    var result = poly.sparse.squarefree_checked(f)
    assert_equal(poly.exact.unit(result), Rational(1))
    assert_true(poly.exact.check(result))
    var blocks = poly.exact.factors(result)
    assert_equal(len(blocks), 2)
    assert_equal(blocks[0][1], 2)
    assert_equal(blocks[1][1], 3)
    assert_equal(poly.sparse.degree(blocks[1][0]), 2)
    var work = Budget()
    var input = poly.sparse.to_core(poly.sparse.asarray([((1, 0), Rational(1)), ((0, 1), Rational(1))]))
    var square = input.pow(2, work)
    var evidence = squarefree_multivariate(square, work)
    evidence.factors[0][1] = 1
    assert_true(not check_sparse_squarefree_product(square, evidence, work))
    evidence.factors[0][1] = 0
    assert_true(not check_sparse_squarefree_product(square, evidence, work))
    evidence.factors[0][1] = 2
    evidence.unit = evidence.unit.neg(work)
    assert_true(not check_sparse_squarefree_product(square, evidence, work))



def boundaries() raises:
    var f = poly.sparse.asarray([((1, 0), 1), ((0, 1), 1)], modulus=3)
    var before = f.to_list()
    for steps in [UInt64(0), UInt64(10), UInt64(100)]:
        for operation in range(4):
            var failed = False
            try:
                if operation == 0:
                    _ = poly.sparse.factor(f, max_steps=steps)
                elif operation == 1:
                    _ = poly.sparse.factor_checked(f, max_steps=steps)
                elif operation == 2:
                    _ = poly.sparse.squarefree(f, max_steps=steps)
                else:
                    _ = poly.sparse.squarefree_checked(f, max_steps=steps)
            except e:
                failed = e.kind == PolyErrorKind.LIMIT
            assert_true(failed)
    assert_equal(f.to_list(), before)
    # Reject the dense image before allocation when its degree box is too large.
    var huge = poly.sparse.asarray([((1000000, 1000000), 1), ((0, 0), 1)], modulus=2)
    var failed = False
    try:
        _ = poly.sparse.factor(huge, max_steps=10000)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var maximum_power = Int("4294967295")
    var overflow = poly.sparse.asarray([((maximum_power, maximum_power), 1)], modulus=2)
    failed = False
    try:
        _ = poly.sparse.factor(overflow, max_steps=10000)
    except e:
        failed = e.kind == PolyErrorKind.CAPACITY or e.kind == PolyErrorKind.LIMIT
    assert_true(failed)


def core_orders() raises:
    var work = Budget()
    var source = poly.sparse.to_core(poly.sparse.asarray([((2, 0), 1), ((0, 2), 1)], modulus=5), variables=[Int64(7), Int64(2)])
    var gaussian = poly.sparse.to_core(poly.sparse.asarray([((2, 0), ExactComplex(1)), ((0, 2), ExactComplex(1))]), variables=[Int64(7), Int64(2)])
    for kind in [MonomialOrder.LEX, MonomialOrder.GRLEX, MonomialOrder.GREVLEX]:
        var order = MonomialOrder(kind, [Int64(7), Int64(2)], work)
        var f = source.reorder(order, work)
        var result = factor_multivariate_finite_field(f, work, 42)
        assert_equal(len(result.factors), 2)
        assert_true(check_multivariate_factorization(f, result, work))
        for entry in result.factors:
            assert_true(entry[0].order == order)
            assert_equal(entry[0].terms[len(entry[0].terms) - 1].coefficient.value, UInt64(1))
        var bad = result.copy()
        bad.factors[0][1] += 1
        assert_true(not check_multivariate_factorization(f, bad, work))
        var square = f.pow(2, work)
        var blocks = squarefree_multivariate(square, work)
        assert_true(check_sparse_squarefree_product(square, blocks, work))
        assert_equal(blocks.factors[0][1], 2)
        assert_true(blocks.factors[0][0].order == order)
        var g = gaussian.reorder(order, work)
        var complex_factors = factor_multivariate_gaussian(g, work)
        assert_true(check_multivariate_factorization(g, complex_factors, work))
        assert_equal(len(complex_factors.factors), 2)
    # KERNEL storage order is not multiplicative: leading scalars must be
    # recovered from normalized factors, not assumed to multiply as in lex.
    var x = poly.sparse.to_core(poly.sparse.asarray([((1, 0), 1)], modulus=5))
    var mixed = poly.sparse.to_core(poly.sparse.asarray([((1, 0), 1), ((0, 1), 2)], modulus=5))
    var kernel = x.mul(mixed.pow(2, work), work).reorder(MonomialOrder(), work)
    var kernel_factors = factor_multivariate_finite_field(kernel, work)
    assert_true(check_multivariate_factorization(kernel, kernel_factors, work))
    assert_equal(kernel_factors.unit.value, UInt64(4))
    var kernel_blocks = squarefree_multivariate(kernel, work)
    assert_true(check_sparse_squarefree_product(kernel, kernel_blocks, work))
    assert_equal(kernel_blocks.unit.value, UInt64(4))
    for entry in kernel_blocks.factors:
        assert_equal(entry[0].terms[len(entry[0].terms) - 1].coefficient.value, UInt64(1))
    # A genuinely irreducible polynomial with a reducible substitution image.
    var irreducible = poly.sparse.asarray([((1, 1), 1), ((0, 0), 1)], modulus=2)
    var factors = poly.sparse.factor(irreducible)
    assert_equal(len(factors[1]), 1)
    assert_equal(factors[1][0][0].to_list(), irreducible.to_list())
    var constant = poly.sparse.asarray([((0, 0), ExactComplex(2, -1))])
    assert_equal(poly.sparse.factor(constant)[0], ExactComplex(2, -1))
    var complex_product = poly.sparse.multiply(poly.sparse.asarray([((1, 0), ExactComplex(1)), ((0, 1), ExactComplex(1))]), constant)
    var complex_result = poly.sparse.factor_checked(complex_product)
    assert_equal(poly.exact.unit(complex_result), ExactComplex(2, -1))
    assert_true(poly.exact.check(complex_result))
    var tuple_input = poly.sparse.asarray([((0,), 1), ((2,), 1)], modulus=2)
    var tuple_proof = poly.sparse.squarefree_checked(tuple_input)
    assert_true(poly.exact.check(tuple_proof))
    assert_equal(poly.exact.factors(tuple_proof)[0][1], 2)


def main() raises:
    finite_checks()
    complex_checks()
    rational_blocks()
    boundaries()
    core_orders()
    print("array sparse factor checks passed")
