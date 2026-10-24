"""Exact algebraic values and original witnesses survive dense/sparse text round trips."""
from poly import dense, sparse, algebraic, coefficients
from polynomials import Budget, QQ, PolyError, PolyErrorKind
from apn_mojo import Rational
from std.testing import assert_equal, assert_true


def main() raises:
    # Value tag 2 + rational flag, then the endpoints and the polynomial.
    var fixture = "POLYNOMIALS-ALGEBRAIC 1\nDENSE\nAA\n0\n1\n3\n2\n2\n2\n-2\n1\n"
    var hand = dense.loads[algebraic.Real](fixture)
    assert_true(hand[0] == algebraic.asreal(2))
    assert_equal(dense.dumps(hand), fixture)
    # Records written before Descartes witnesses carried Sturm sequences.
    var stale = False
    try:
        _ = dense.loads[algebraic.Real]("POLYNOMIALS-ALGEBRAIC 1\nDENSE\nAA\n0\n1\n1\n2\n2\n2\n-2\n1\n0\n0\n")
    except e:
        stale = e.kind == PolyErrorKind.ENCODING
    assert_true(stale)
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray(coefficients((a.copy(), 0, Rational(1, 3), 1)))
    var text = dense.dumps(p)
    var copy = dense.loads[algebraic.Real](text)
    assert_equal(copy.to_list(), p.to_list())
    assert_equal(dense.dumps(copy), text)
    assert_true(algebraic.check(copy[0]))
    assert_equal(copy[0].witness().polynomial.coefficients, a.witness().polynomial.coefficients)
    assert_true(copy[0].witness().low == a.witness().low and copy[0].witness().high == a.witness().high)
    var z = algebraic.ascomplex(a, dense.root_between([-3, 0, 1], 1, 2))
    var q = dense.asarray(coefficients((z.copy(), 2, Rational(1, 7))))
    var qcopy = dense.loads[algebraic.Complex](dense.dumps(q))
    assert_equal(q.to_list(), qcopy.to_list())
    assert_true(algebraic.check(qcopy[0]))
    assert_equal(dense.dumps(qcopy), dense.dumps(q))
    var sp = sparse.asarray((((0, 2, 0), z.copy()), ((1, 0, 0), 1)))
    var spcopy = sparse.loads[algebraic.Complex, ndim=3](sparse.dumps(sp))
    assert_equal(spcopy.to_list(), sp.to_list())
    assert_equal(sparse.dumps(spcopy), sparse.dumps(sp))
    var c = sparse.asarray([((), a.copy())])
    assert_equal(sparse.loads[algebraic.Real, ndim=0](sparse.dumps(c)).to_list(), c.to_list())
    var empty = dense.asarray(List[algebraic.Complex]())
    assert_equal(dense.loads[algebraic.Complex](dense.dumps(empty)).domain(), "QQbar")
    var sempty = sparse.asarray(List[Tuple[Tuple[Int, Int], algebraic.Real]]())
    assert_equal(sparse.loads[algebraic.Real, ndim=2](sparse.dumps(sempty)).domain(), "AA")
    # A valid nonminimal defining polynomial must be retained, not silently
    # replaced by the minimal polynomial used for scalar equality.
    var witness = a.witness().copy()
    var budget = Budget()
    witness.polynomial = witness.polynomial.mul(dense.to_core([Rational(-3), Rational(1)]), budget)
    var nonminimal = algebraic.asreal(witness)
    var restored = dense.loads[algebraic.Real](dense.dumps([nonminimal.copy()]))[0]
    assert_true(restored == a)
    assert_equal(restored.witness().polynomial.coefficients, witness.polynomial.coefficients)
    assert_equal(dense.dumps([restored.copy()]), dense.dumps([nonminimal.copy()]))
    # Loading restores value evidence; it does not invent arithmetic history.
    var arithmetic = algebraic.add(a, 1)
    var detached = dense.loads[algebraic.Real](dense.dumps([arithmetic.copy()]))[0]
    assert_true(detached == arithmetic)
    assert_true(algebraic.check(detached))
    assert_true(not detached.evidence())
    print("array algebraic serialization checks passed")
