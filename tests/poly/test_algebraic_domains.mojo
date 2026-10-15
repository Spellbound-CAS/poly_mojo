"""General exact algebraic coefficients, promotion, cancellation, and field identities."""
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, NoContext, PolyErrorKind
from poly import dense, sparse, algebraic, exact, coefficients, to_dense, to_sparse
from std.testing import assert_equal, assert_true


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var another = dense.root_between([-6, 2, 3, -1], 1, 2)
    assert_true(a == another)
    assert_true(a != dense.root_between([-2, 0, 1], -2, 0))
    assert_true(algebraic.multiply(a, a) == algebraic.asreal(2))
    var positive = dense.root_between([-2, 0, 1], -1, 2)
    var negative = dense.root_between([-2, 0, 1], -2, 1)
    assert_true(positive != negative)
    assert_true(positive == a)
    var b = dense.root_between([-3, 0, 1], 1, 2)
    var c = dense.root_between([-6, 0, 1], 2, 3)
    assert_true(dense.multiply([a.copy()], [b.copy()])[0] == c)
    var p = dense.asarray(coefficients((algebraic.negative(a), 1)))
    assert_equal(p.domain(), "AA")
    assert_equal(dense.asarray(List[algebraic.Real]()).domain(), "AA")
    assert_equal(len(dense.subtract(p, p)), 0)
    assert_equal(dense.subtract(p, p).domain(), "AA")
    assert_true(dense.evaluate(p, a).is_zero())
    var plus = dense.asarray(coefficients((a.copy(), 1)))
    var product = dense.multiply(p, plus)
    assert_equal(product.to_list(), dense.asarray(coefficients((-2, algebraic.asreal(0), 1))).to_list())
    var division = dense.divide_checked(product, p)
    assert_true(exact.check(division))
    assert_equal(exact.quotient(division).to_list(), plus.to_list())
    division._data[].proof.quotient.coefficients[0] = algebraic.asreal(7)
    assert_true(not exact.check(division))
    assert_equal(dense.derivative(p)[0], algebraic.asreal(1))
    assert_equal(dense.derivative(dense.integral(p)).to_list(), p.to_list())
    assert_equal(dense.compose(p, plus).to_list(), dense.asarray(coefficients((algebraic.asreal(0), 1))).to_list())
    var core = dense.to_core(p)
    assert_equal(dense.from_core(core).to_list(), p.to_list())
    assert_equal(to_dense(to_sparse(p)).to_list(), p.to_list())
    var z = algebraic.ascomplex(a, 1)
    assert_true(algebraic.check(z))
    assert_true(algebraic.multiply(z, algebraic.conjugate(z)) == algebraic.ascomplex(3))
    assert_true(algebraic.multiply(z, algebraic.reciprocal(z)) == algebraic.ascomplex(1))
    var mixed = coefficients((a.copy(), ExactComplex(Rational(0), Rational(1))))
    assert_equal(mixed.domain(), "QQbar")
    var cp = dense.asarray(coefficients((algebraic.negative(z), 1)))
    assert_equal(cp.domain(), "QQbar")
    assert_true(dense.evaluate(cp, z).is_zero())
    assert_equal(dense.add(cp, p).domain(), "QQbar")
    assert_equal(dense.multiply(cp, [0]).domain(), "QQbar")
    assert_equal(len(dense.multiply(cp, [0])), 0)
    var terms = sparse.asarray((((1, 0), a.copy()), ((0, 1), 1), ((0, 0), Rational(1, 2))))
    assert_equal(terms.domain(), "AA")
    assert_true(sparse.evaluate(terms, (a.copy(), 0)) == algebraic.asreal(Rational(5, 2)))
    assert_equal(sparse.subtract(terms, terms).domain(), "AA")
    assert_equal(len(sparse.subtract(terms, terms)), 0)
    var complex_terms = sparse.asarray((((1, 0), z.copy()), ((0, 1), 1)))
    assert_equal(complex_terms.domain(), "QQbar")
    assert_true(sparse.evaluate(complex_terms, (algebraic.ascomplex(1), z.copy())) == algebraic.multiply(2, z))
    var limited = False
    try:
        _ = dense.multiply(p, plus, max_steps=0)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(p[0], algebraic.negative(a))
    var zero_inverse = False
    try:
        _ = algebraic.reciprocal(algebraic.ascomplex(0))
    except e:
        zero_inverse = e.kind == PolyErrorKind.DOMAIN
    assert_true(zero_inverse)
    print("array algebraic domain checks passed")
