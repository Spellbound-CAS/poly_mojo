"""Field algorithms over complex algebraic coefficients and explicit unsupported calls."""
from poly import dense, sparse, algebraic, coefficients, exact, series, basis
from polynomials import PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var z = algebraic.ascomplex(a, 1)
    var p = dense.asarray(coefficients((algebraic.negative(z), 1)))
    var squared = dense.multiply(p, p)
    var divided = dense.divide_checked(squared, p)
    assert_true(exact.check(divided))
    assert_equal(exact.quotient(divided).to_list(), p.to_list())
    var common = dense.gcd_checked(squared, p)
    assert_true(exact.check(common))
    assert_equal(exact.value(common).to_list(), p.to_list())
    var sf = dense.squarefree_checked(squared)
    assert_true(exact.check(sf))
    assert_equal(exact.factors(sf)[0][1], 2)
    assert_true(exact.unit(sf) == algebraic.ascomplex(1))
    assert_equal(dense.derivative(dense.integral(p)).to_list(), p.to_list())
    var nodes = coefficients((algebraic.ascomplex(0), 1))
    var values = coefficients((z.copy(), algebraic.add(z, 1)))
    var interpolation = dense.interpolate_checked(nodes, values)
    assert_true(exact.check(interpolation))
    assert_equal(exact.value(interpolation).to_list(), dense.asarray(coefficients((z.copy(), 1))).to_list())
    var term = sparse.asarray((((1, 0), 1), ((0, 1), z.copy())))
    var sparse_sf = sparse.squarefree_checked(sparse.multiply(term, term))
    assert_true(exact.check(sparse_sf))
    assert_equal(exact.factors(sparse_sf)[0][1], 2)
    var core = sparse.to_core(term)
    assert_equal(sparse.from_core[ndim=2](core).to_list(), term.to_list())
    assert_equal(series.asarray(p, order=3).domain(), "QQbar")
    assert_equal(basis.asarray(p, basis="chebyshev").domain(), "QQbar")
    assert_equal(dense.factor(p)[1][0][0].to_list(), p.to_list())
    assert_true(exact.check(dense.factor_checked(p)))
    assert_equal(len(sparse.factor(term)[1]), 1)
    assert_true(exact.check(sparse.factor_checked(term)))
    assert_equal(dense.loads[algebraic.Complex](dense.dumps(p)).to_list(), p.to_list())
    var wrong_axis = False
    try:
        _ = algebraic.coefficients(term)
    except e:
        wrong_axis = e.kind == PolyErrorKind.DOMAIN
    assert_true(wrong_axis)
    var duplicate = False
    try:
        _ = sparse.specialize(term, [(1, z.copy()), (1, z.copy())])
    except e:
        duplicate = e.kind == PolyErrorKind.DOMAIN
    assert_true(duplicate)
    print("array complex algebraic operation checks passed")
