"""Equivalent exact coordinate witnesses share a stream and keep their archive evidence."""
from poly import cad, dense, algebraic, io
from polynomials import Budget, QQ, sturm_sequence
from apn_mojo import Rational
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var budget = Budget()
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var witness = a.witness().copy()
    witness.polynomial = witness.polynomial.mul(dense.to_core([Rational(-3), Rational(1)]), budget)
    witness.sequence = sturm_sequence(witness.polynomial, budget)
    var same = algebraic.asreal(witness)
    assert_true(same == a)
    var p = [((1, 0), 1), ((0, 1), -1)]
    var checker = cad.checker([p.copy()], fixed=[(1, a.copy())])
    var receipt = cad.decompose_to([p.copy()], checker, budget, fixed=[(1, same.copy())])
    assert_true(checker.finish(receipt, budget))
    assert_equal(receipt.cells, 3)
    var reverse = cad.checker([p.copy()], fixed=[(1, witness.copy())])
    var reverse_receipt = cad.decompose_to([p.copy()], reverse, budget, fixed=[(1, a.copy())])
    assert_true(reverse.finish(reverse_receipt, budget))
    var result = cad.decompose([p.copy()], fixed=[(1, same.copy())])
    assert_true(cad.check(result))
    assert_equal(cad.cells(result)[1].signs, [0])
    assert_true(cad.cells(result)[1].sample[0] == a)
    assert_equal(cad.fixed_coordinates(result)[0][1].witness().polynomial.degree(), 3)
    var saved = io.dumps(result)
    var restored = io.loads[type_of(result)](saved)
    assert_true(cad.check(restored))
    assert_equal(io.dumps(restored), saved)
    assert_equal(cad.fixed_coordinates(restored)[0][1].witness().polynomial.coefficients, witness.polynomial.coefficients)
    # The requested bounds can be much narrower than the stack's root sample.
    var narrow = dense.root_between([-2, 0, 1], Rational(1414, 1000), Rational(1415, 1000))
    var narrow_checker = cad.checker([p.copy()], fixed=[(1, narrow.copy())])
    var narrow_receipt = cad.decompose_to([p.copy()], narrow_checker, budget, fixed=[(1, a.copy())])
    assert_true(narrow_checker.finish(narrow_receipt, budget))
    var narrow_result = cad.decompose([p.copy()], fixed=[(1, narrow.copy())])
    assert_true(cad.check(narrow_result))
    var wrong = cad.decompose([p.copy()], fixed=[(1, same.copy())])
    wrong._data[].fixed[0][1] = algebraic.negative(narrow)
    assert_false(cad.check(wrong))
    # Even an irrational-style witness can denote an exact rational value.
    witness.polynomial = dense.to_core([Rational(-1), Rational(1), Rational(-1), Rational(1)])
    witness.sequence = sturm_sequence(witness.polynomial, budget)
    witness.low = QQ(Rational(0))
    witness.high = QQ(Rational(2))
    var one = algebraic.asreal(witness)
    assert_true(one == algebraic.asreal(1))
    var q = [(0, -1), (1, 1)]
    var rational_checker = cad.checker([q.copy()], fixed=[(0, 1)])
    var rational_receipt = cad.decompose_to([q.copy()], rational_checker, budget, fixed=[(0, one.copy())])
    assert_true(rational_checker.finish(rational_receipt, budget))
    var solution = cad.solve([cad.zero(q)], fixed=[(0, one.copy())])
    assert_true(solution.satisfiable() and cad.check(solution))
    var full = cad.cover([cad.zero(q)], fixed=[(0, one.copy())])
    assert_true(full.satisfiable() and cad.check(full))
    var impossible = cad.solve([cad.nonzero(q)], fixed=[(0, one.copy())])
    assert_false(impossible.satisfiable())
    assert_true(cad.check(impossible))
    print("fixed CAD witness checks passed")
