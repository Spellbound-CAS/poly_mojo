"""Exact sign queries, determinant proofs, and sparse algorithm choices."""
from std.testing import assert_equal, assert_true
from apn_mojo import Rational
from poly import dense, sparse, exact


def main() raises:
    # x^2 - 2 stays positive throughout [2, 3].
    var sign = dense.sign_on_interval([-2, 0, 1], 2, 3)
    assert_equal(sign.value(), 1)
    print("Positive throughout [2, 3]")

    # The determinant gives Res(x^2 + 1, x^2 - 2) = 9.
    var resultant = dense.resultant_certified([1, 0, 1], [-2, 0, 1])
    assert_true(exact.check(resultant))
    assert_equal(exact.value(resultant)[0], Rational(9))
    print("Resultant 9 verified from its determinant")

    var h = sparse.asarray([((1, 0), 1), ((0, 1), 1)])  # x + y
    var a = sparse.multiply(h, [((1, 0), 1), ((0, 0), 1)])
    var b = sparse.multiply(h, [((0, 1), 1), ((0, 0), 1)])
    var common = sparse.gcd_checked(a, b, method="modular")
    assert_true(exact.check(common))
    assert_equal(exact.value(common).to_list(), h.to_list())
    print("Sparse modular GCD: x + y")

    # Dividing x by x + y^2 illustrates how monomial order affects reduction.
    var x = sparse.asarray([((1, 0), 1)])
    var divisor = sparse.asarray([((1, 0), 1), ((0, 2), 1)])
    var lex = sparse.divide_checked(x, divisor)
    var graded = sparse.divide_checked(x, divisor, order="grlex")
    assert_true(exact.check(lex) and exact.check(graded))
    assert_equal(exact.remainder(lex).to_list(), sparse.asarray([((0, 2), Rational(-1))]).to_list())
    assert_equal(exact.remainder(graded).to_list(), sparse.asarray([((1, 0), Rational(1))]).to_list())
    print("Division remainders: -y^2 under lex, x under grlex")
