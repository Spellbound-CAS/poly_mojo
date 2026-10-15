"""Full determinant evidence across inferred exact coefficient fields."""
from poly import dense, sparse, exact, algebraic, coefficients, modint
from polynomials import Budget
from apn_mojo import ExactComplex, Rational
from std.testing import assert_equal, assert_true, assert_false

def main() raises:
    var gaussian = dense.resultant_certified([ExactComplex(0, 1)], [0, 1])
    assert_true(exact.check(gaussian))
    assert_equal(exact.value(gaussian)[0], ExactComplex(0, 1))
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray(coefficients((algebraic.negative(a), 1)))
    var proof = dense.resultant_certified(p, [0, 1])
    assert_equal(algebraic.compare(exact.value(proof)[0], a), 0)
    assert_true(exact.check(proof))
    proof._data[].proof.resultant = sparse.to_core(sparse.asarray([(0, algebraic.asreal(0))]))
    assert_false(exact.check(proof))
    var gf = dense.resultant_certified([modint(1, 5), modint(2, 5), modint(1, 5)], [modint(2, 5), modint(1, 5)])
    assert_true(exact.check(gf))
    assert_equal(exact.value(gf)[0].value, UInt64(1))
    var z = algebraic.ascomplex(a, 1)
    var complex = dense.resultant_certified(coefficients((algebraic.negative(z), 1)), [0, 1])
    assert_true(exact.check(complex))
    assert_true(exact.value(complex)[0] == z)
    var sq = sparse.asarray([((1, 0), a.copy()), ((0, 1), algebraic.asreal(1))])
    var right = sparse.asarray([((1, 0), algebraic.asreal(1))])
    var eliminated = sparse.resultant_certified(sq, right, axis=0)
    assert_true(exact.check(eliminated))
    assert_equal(algebraic.compare(sparse.evaluate(exact.value(eliminated), (0, 3)), -3), 0)
    print("completion determinant checks passed")
