"""Defective PRS scaling over polynomial coefficient rings, checked against minors."""
from apn_mojo import Rational
from polynomials import Budget, QQ, ZZ, NoContext, SparsePoly, DensePoly, MonomialOrder, resultant, subresultants, check_resultant_identity, substitute
from support.corpus import as_zz
from support.determinant import subresultant_coefficients
from std.testing import assert_equal, assert_true


def cubic(z: SparsePoly[QQ], a: Int, b: Int, c: Int, d: Int, mut budget: Budget) raises -> SparsePoly[QQ]:
    var constant = SparsePoly[QQ].constant(QQ(Rational(d)), z.order, budget)
    return z.pow(3, budget).scale(QQ(Rational(a)), budget).add(z.pow(2, budget).scale(QQ(Rational(b)), budget), budget).add(z.scale(QQ(Rational(c)), budget), budget).add(constant, budget)


def main() raises:
    var budget = Budget()
    var order = MonomialOrder()
    var x = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var z = SparsePoly[QQ].variable(1, order, NoContext(), budget)
    var f = x.pow(6, budget).scale(QQ(Rational(64)), budget).sub(x.pow(4, budget).scale(QQ(Rational(96)), budget), budget).add(x.pow(2, budget).scale(QQ(Rational(36)), budget), budget).sub(SparsePoly[QQ].constant(QQ(Rational(31)), order, budget), budget)
    var g = z.pow(6, budget).scale(QQ(Rational(64)), budget).sub(z.pow(4, budget).mul(x.pow(2, budget), budget).scale(QQ(Rational(96)), budget), budget).add(z.pow(2, budget).mul(x.pow(4, budget), budget).scale(QQ(Rational(36)), budget), budget).sub(x.pow(6, budget).scale(QQ(Rational(31)), budget), budget)
    # The final PRS polynomial has degree 26 in z; the resultant has degree 36.
    # Their quotient has a nonconstant denominator, but the scaled cofactors
    # and every Sylvester minor are polynomials.
    var expected = cubic(z, 64, -48, 0, 31, budget).pow(4, budget).mul(cubic(z, 64, 48, 0, -31, budget).pow(4, budget), budget).mul(cubic(z, 64, -96, 36, -31, budget).pow(2, budget), budget).mul(cubic(z, 64, 96, 36, 31, budget).pow(2, budget), budget)
    var proof = resultant(f, g, 0, budget)
    assert_equal(proof.resultant, expected)
    assert_true(check_resultant_identity(f, g, proof, budget))
    var fz = as_zz(f, budget)
    var gz = as_zz(g, budget)
    var integer_proof = resultant(fz, gz, 0, budget)
    assert_equal(integer_proof.resultant, as_zz(expected, budget))
    assert_true(check_resultant_identity(fz, gz, integer_proof, budget))
    var family = subresultants(f, g, 0, budget)
    assert_equal(family[0], expected)
    for value in [-2, -1, 0, 1, 2]:
        var q = SparsePoly[QQ].constant(QQ(Rational(value)), order, budget)
        var fs = DensePoly[ZZ].from_sparse(as_zz(substitute(f, 1, q, budget), budget), budget, variable=Int64(0))
        var gs = DensePoly[ZZ].from_sparse(as_zz(substitute(g, 1, q, budget), budget), budget, variable=Int64(0))
        for index in range(6):
            var coefficients = subresultant_coefficients(fs, gs, index)
            var minor = DensePoly[ZZ].from_coefficients(0, NoContext(), coefficients^, budget)
            var actual = DensePoly[ZZ].from_sparse(as_zz(substitute(family[index], 1, q, budget), budget), budget, variable=Int64(0))
            assert_equal(actual, minor)
    print("defective resultant scaling checks passed")
