"""Projection, query-point binding, and Lazard orders with AA coefficients."""
from poly import cad, dense, sparse, algebraic
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = sparse.asarray([(0, algebraic.negative(a)), (1, algebraic.asreal(1))])
    var projection = cad.project([p.copy()])
    assert_true(cad.check(projection))
    var projected = cad.projection_polynomials(projection)
    assert_true(len(projected) > 0)
    for polynomial in projected:
        assert_equal(polynomial.domain(), "AA")
    var query = cad.cell_at([p.copy()], (a.copy(),))
    assert_true(cad.check(query))
    assert_equal(cad.signs(query), [0])
    assert_equal(cad.sample(query), [a.copy()])
    var positive = cad.cell_at([p.copy()], [algebraic.asreal(2)])
    assert_true(cad.check(positive))
    assert_equal(cad.signs(positive), [1])
    # An equally valid opposite coefficient root must fail the facade binding.
    var wrong = cad.cell_at([p.copy()], [algebraic.asreal(2)])
    wrong._data[].point[0] = a.witness().copy()
    assert_false(cad.check(wrong))
    var x = sparse.asarray([((1, 0), algebraic.asreal(1)), ((0, 0), algebraic.negative(a))])
    var y = sparse.asarray([((0, 1), algebraic.asreal(1)), ((0, 0), algebraic.asreal(-1))])
    var product = sparse.multiply(sparse.power(x, 2), y)
    var lazard = cad.lazard_evaluate(product, [(0, a.copy())], axis=1)
    assert_true(cad.check(lazard))
    assert_equal(cad.vanishing_orders(lazard), [(0, 2)])
    var lifting = cad.lifting_coefficients(lazard)
    assert_equal(len(lifting), 2)
    assert_true(lifting[0] == algebraic.negative(lifting[1]))
    assert_equal(algebraic.compare(lifting[1], 0), 1)
    var real_rational = sparse.asarray([(0, algebraic.asreal(-1)), (1, algebraic.asreal(1))])
    var rational_projection = cad.project([real_rational.copy()], lazard=False)
    assert_true(cad.check(rational_projection))
    assert_equal(cad.projection_polynomials(rational_projection)[0].domain(), "AA")
    # Independent algebraic coefficient fields share one exact CAD problem.
    var b = dense.root_between([-3, 0, 1], 1, 2)
    var q = sparse.asarray([(0, algebraic.negative(b)), (1, algebraic.asreal(1))])
    var separate_fields = cad.cell_at([p.copy(), q.copy()], [algebraic.asreal(2)])
    assert_true(cad.check(separate_fields))
    assert_equal(cad.signs(separate_fields), [1, 1])
    assert_equal(len(cad.sample(separate_fields)), 1)
    print("array algebraic CAD specialized checks passed")
