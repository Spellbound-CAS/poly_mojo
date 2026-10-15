"""Checked exact identities across rational, Gaussian and prime-field domains."""
from std.testing import assert_true, assert_false, assert_equal
from apn_mojo import Rational, ExactComplex
from polynomials import PolyErrorKind
from poly import dense, sparse, exact, modint


def main() raises:
    var p = sparse.asarray([((0,), -1), ((2,), 1)], modulus=5)
    var f = sparse.factor_checked(p)
    assert_true(exact.check(f))
    assert_equal(exact.unit(f).value, UInt64(1))
    assert_equal(len(exact.factors(f)), 2)
    var zero = sparse.factor_checked(sparse.asarray(List[Tuple[Int, Int]](), modulus=5))
    assert_true(exact.check(zero))
    assert_equal(len(exact.factors(zero)), 0)
    var scalar = dense.factor_checked(dense.asarray([2], modulus=5))
    assert_true(exact.check(scalar))
    assert_equal(exact.unit(scalar).value, UInt64(2))
    var iq = dense.interpolate_checked([modint(0, 5), modint(1, 5)], [modint(2, 5), modint(4, 5)])
    assert_true(exact.check(iq))
    assert_equal(exact.value(iq)[0].value, UInt64(2))
    assert_equal(exact.value(iq)[1].value, UInt64(2))
    var z = [ExactComplex(1, 1), ExactComplex(1, 0)]
    assert_true(exact.check(dense.gcd_checked(z, z)))
    var c = dense.coprimality(z, [ExactComplex(1, 0)])
    assert_true(Bool(c) and exact.check(c.value()))
    var square = sparse.asarray([(0, ExactComplex(0, 2)), (1, ExactComplex(2, 2)), (2, ExactComplex(1, 0))])
    var sq = sparse.squarefree_checked(square)
    assert_true(exact.check(sq))
    assert_equal(exact.factors(sq)[0][1], 2)
    var rational = sparse.factor_checked([((1, 0), Rational(1, 2)), ((0, 0), Rational(-1, 2))])
    assert_true(exact.check(rational))
    assert_equal(exact.unit(rational), Rational(1, 2))
    var result = dense.resultant_checked([Rational(-2), Rational(1)], [Rational(-3), Rational(1)])
    assert_true(exact.check(result))
    var multipliers = exact.cofactors(result)
    var combined = dense.add(dense.multiply(multipliers[0], [Rational(-2), Rational(1)]), dense.multiply(multipliers[1], [Rational(-3), Rational(1)]))
    assert_equal(combined.to_list(), exact.value(result).to_list())
    var gaussian_factors = dense.factor_checked(z)
    assert_true(exact.check(gaussian_factors))
    assert_equal(exact.unit(gaussian_factors), ExactComplex(1))
    assert_equal(exact.factors(gaussian_factors)[0][0].to_list(), dense.asarray(z).to_list())
    var multivariate = sparse.asarray([((1, 0), 1)], modulus=5)
    var finite_factors = sparse.factor_checked(multivariate)
    assert_true(exact.check(finite_factors))
    assert_equal(exact.factors(finite_factors)[0][0].to_list(), multivariate.to_list())
    print("array parity domain checks passed")
