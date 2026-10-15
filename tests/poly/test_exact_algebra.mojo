"""Exact algebra identities through coefficient and term sequences."""
import poly
from apn_mojo import Integer, Rational
from polynomials import PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    assert_equal(poly.dense.content([-6, 12, -6]), Integer(6))
    assert_equal(poly.dense.primitive_part([-6, 12, -6]).to_list(), poly.dense.asarray([-1, 2, -1]).to_list())
    assert_equal(poly.dense.gcd([-2, 0, 2], [-4, 4]).to_list(), poly.dense.asarray([-2, 2]).to_list())
    var bezout = poly.dense.extended_gcd([-1, 0, 1], [-1, 1])
    var reconstructed = poly.dense.add(poly.dense.multiply(bezout[1], [-1, 0, 1]), poly.dense.multiply(bezout[2], [-1, 1]))
    assert_equal(bezout[0].to_list(), reconstructed.to_list())
    var quotient = poly.dense.divide_exact([1, 2, 1], [1, 1])
    assert_true(Bool(quotient))
    assert_equal(quotient.value().to_list(), poly.dense.asarray([1, 1]).to_list())
    assert_true(not poly.dense.divide_exact([1, 1], [2]))
    assert_equal(poly.dense.resultant([-2, 1], [-3, 1]), Integer(-1))
    var sub = poly.dense.subresultants([-2, 1], [-3, 1])
    assert_equal(len(sub), 1)
    assert_equal(sub[0][0], Integer(-1))
    var pseudo = poly.dense.pseudo_divmod([1, 0, 1], [1, 2])
    assert_equal(pseudo[0], Integer(4))
    assert_equal(poly.dense.add(poly.dense.multiply(pseudo[1], [1, 2]), pseudo[2]).to_list(), poly.dense.asarray([4, 0, 4]).to_list())
    var p = poly.sparse.asarray([((2, 1), 2), ((1, 2), 2)])
    var content = poly.sparse.content(p, axis=0)
    assert_equal(content.to_list(), poly.sparse.asarray([((0, 1), 2)]).to_list())
    var primitive = poly.sparse.primitive_part(p)
    assert_equal(poly.sparse.multiply(primitive, content).to_list(), p.to_list())
    var gcd = poly.sparse.gcd(p, content)
    assert_equal(gcd.to_list(), content.to_list())
    var divided = poly.sparse.divide_exact(p, content)
    assert_true(Bool(divided))
    assert_equal(divided.value().to_list(), primitive.to_list())
    var x_y = [((1, 0), 1), ((0, 1), -1)]
    var x_1 = [((1, 0), 1), ((0, 0), -1)]
    assert_equal(poly.sparse.resultant(x_y, x_1).to_list(), poly.sparse.asarray([((0, 1), 1), ((0, 0), -1)]).to_list())
    var sequence = poly.sparse.subresultants(x_y, x_1)
    assert_equal(sequence[0].to_list(), poly.sparse.resultant(x_y, x_1).to_list())
    var limited = False
    try:
        _ = poly.dense.gcd([1, 2], [1], max_steps=0)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    print("array exact algebra checks passed")
