"""Independent nonlinear Gaussian factors and repeated norm reductions."""
from apn_mojo import ExactComplex
from polynomials import QQi, NoContext, DensePoly, Budget, PolyErrorKind, factor_gaussian, check_factorization
from std.testing import assert_equal, assert_true


def main() raises:
    var work = Budget()
    var f = DensePoly[QQi].from_coefficients(3, NoContext(), [QQi(ExactComplex(1)), QQi(ExactComplex(0)), QQi(ExactComplex(0)), QQi(ExactComplex(0)), QQi(ExactComplex(1))], work)
    # x^4+1 = (x^2-i)(x^2+i), with both quadratics irreducible over QQ(i).
    var result = factor_gaussian(f, work)
    assert_equal(len(result.factors), 2)
    assert_true(check_factorization(f, result, work))
    for i in range(2):
        ref factor = result.factors[i][0]
        assert_equal(factor.variable, Int64(3))
        assert_equal(factor.degree(), 2)
        assert_equal(factor.coefficients[0].value, ExactComplex(0, 2 * i - 1))
        assert_equal(factor.coefficients[1].value, ExactComplex(0))
        assert_equal(factor.coefficients[2].value, ExactComplex(1))
    var repeated = f.pow(2, work).scale(QQi(ExactComplex(1, 1)), work)
    var blocks = factor_gaussian(repeated, work)
    assert_equal(blocks.unit.value, ExactComplex(1, 1))
    assert_equal(len(blocks.factors), 2)
    assert_equal(blocks.factors[0][1], 2)
    assert_equal(blocks.factors[1][1], 2)
    assert_true(check_factorization(repeated, blocks, work))
    # x^4+4 has four Gaussian-rational roots: +/-1 +/- i.
    var split = f.copy()
    split.coefficients[0] = QQi(ExactComplex(4))
    var linear = factor_gaussian(split, work)
    assert_equal(len(linear.factors), 4)
    for i in range(4):
        assert_equal(linear.factors[i][0].degree(), 1)
        assert_equal(linear.factors[i][0].coefficients[0].value, ExactComplex(2 * (i // 2) - 1, 2 * (i % 2) - 1))
        assert_equal(linear.factors[i][1], 1)
    assert_true(check_factorization(split, linear, work))
    var malformed = f.copy()
    malformed.coefficients.append(QQi(ExactComplex(0)))
    var rejected = False
    try:
        _ = factor_gaussian(malformed, work)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    for steps in [UInt64(0), UInt64(100), UInt64(10000)]:
        var small = Budget(steps)
        var limited = False
        try:
            _ = factor_gaussian(f, small)
        except e:
            limited = e.kind == PolyErrorKind.LIMIT
        assert_true(limited)
        assert_equal(f.coefficients[0].value, ExactComplex(1))
    print("Gaussian factorization checks passed")
