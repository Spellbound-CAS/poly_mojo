"""Irreducibility, normalization and completeness checks reject false certificates."""
from poly import dense, sparse, algebraic, coefficients, exact
from polynomials import DensePoly, NoContext
from std.testing import assert_equal, assert_true, assert_false


def main() raises:
    var p = dense.asarray(coefficients((algebraic.asreal(-1), 0, 1)))
    var result = dense.factor_checked(p)
    var good = result.evidence().copy()
    # This unsplit quadratic has the right product but is reducible over AA.
    result._data[].proof.factors = [(dense.to_core(p), 1)]
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    result._data[].proof.factors[0][1] = 0
    assert_false(exact.check(result))
    result._data[].proof.factors[0][1] = Int.MAX
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    result._data[].proof.factors.swap_elements(0, 1)
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    _ = result._data[].proof.factors.pop()
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    result._data[].proof.unit = algebraic.asreal(2)
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    result._data[].proof.factors[0][0].coefficients[0] = algebraic.asreal(-3)
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    result._data[].proof.factors[0][0].variable = 1
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    # Rescaling a factor and the unit preserves the product, but not monicity.
    result._data[].proof.factors[0][0] = dense.to_core(coefficients((algebraic.asreal(-2), 2)))
    result._data[].proof.unit = algebraic.divide(1, 2)
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    result._data[].proof.factors[0][0] = dense.to_core(coefficients((algebraic.asreal(-2), 2)))
    result._data[].proof.factors[1][0] = dense.to_core([algebraic.divide(1, 2), algebraic.divide(1, 2)])
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    assert_true(exact.check(result))
    var repeated_p = dense.asarray(coefficients((algebraic.asreal(1), -2, 1)))
    var repeated = dense.factor_checked(repeated_p)
    var linear = dense.to_core(coefficients((algebraic.asreal(-1), 1)))
    repeated._data[].proof.factors = [(linear.copy(), 1), (linear.copy(), 1)]
    assert_false(exact.check(repeated))
    # Discriminant zero also fails the real-quadratic irreducibility check.
    repeated._data[].proof.factors = [(dense.to_core(repeated_p), 1)]
    assert_false(exact.check(repeated))
    var complex_p = dense.asarray(coefficients((algebraic.ascomplex(1), 0, 1)))
    var complex_result = dense.factor_checked(complex_p)
    complex_result._data[].proof.factors = [(dense.to_core(complex_p), 1)]
    assert_false(exact.check(complex_result))
    var constant = dense.factor_checked([algebraic.asreal(1)])
    constant._data[].proof.factors = [(dense.to_core([algebraic.asreal(1)]), 1)]
    assert_false(exact.check(constant))
    var zero = dense.factor_checked([algebraic.asreal(0)])
    zero._data[].proof.unit = algebraic.asreal(1)
    assert_false(exact.check(zero))
    var s = sparse.asarray((((0, 2), algebraic.asreal(1)), ((0, 0), -1)))
    var sr = sparse.factor_checked(s)
    var sg = sr.evidence().copy()
    sr._data[].proof.factors = [(sparse.to_core(s), 1)]
    assert_false(exact.check(sr))
    sr._data[].proof = sg.copy()
    sr._data[].proof.factors[0][0] = sparse.to_core(sparse.asarray((((1, 0), algebraic.asreal(1)), ((0, 0), -1))))
    assert_false(exact.check(sr))
    sr._data[].proof = sg.copy()
    sr._data[].proof.factors[0][0] = sparse.to_core(sparse.asarray((((1, 1), algebraic.asreal(1)),)))
    assert_false(exact.check(sr))
    sr._data[].proof = sg.copy()
    assert_true(exact.check(sr))
    sr._data[].proof.factors[0][0] = sparse.to_core(sparse.asarray((((0, 1000000000), algebraic.asreal(1)),)))
    assert_false(exact.check(sr, max_steps=1000000))
    sr._data[].proof = sg.copy()
    # A valid witness with a forged zero cache cannot establish a coefficient.
    var forged = algebraic.asreal(7)
    forged._canonical[].zero = True
    result._data[].proof.unit = forged^
    assert_false(exact.check(result))
    result._data[].proof = good.copy()
    # Output copies and copies of the evidence do not mutate the retained result.
    var outputs = exact.factors(result)
    outputs[0][0] = dense.asarray([algebraic.asreal(99)])
    var detached = result.evidence().copy()
    detached.factors.clear()
    assert_equal(len(exact.factors(result)), 2)
    assert_true(exact.check(result))
    print("array algebraic factor evidence checks passed")
