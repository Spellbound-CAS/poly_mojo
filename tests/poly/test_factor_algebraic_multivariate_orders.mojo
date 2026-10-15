"""Preserve core term orders, variable labels, and normalization scalars."""
from poly import sparse, algebraic
from poly import _algebraic_factor as implementation
from polynomials import Budget, MonomialOrder, check_multivariate_factorization
from std.testing import assert_equal, assert_true, assert_false


def main() raises:
    var budget = Budget()
    var p = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 2), -1)))
    var source = sparse.to_core(p, variables=[Int64(7), Int64(2)])
    for kind in [MonomialOrder.LEX, MonomialOrder.GRLEX, MonomialOrder.GREVLEX]:
        var order = MonomialOrder(kind, [Int64(7), Int64(2)], budget)
        var input = source.reorder(order, budget)
        var output = implementation._factor_sparse_evidence(input, budget)
        assert_true(implementation._check_sparse_evidence(input, output[0], output[1], budget))
        assert_true(check_multivariate_factorization(input, output[0], budget))
        assert_equal(len(output[0].factors), 2)
        for entry in output[0].factors:
            assert_true(entry[0].order == order)
            assert_true(entry[0].terms[len(entry[0].terms) - 1].coefficient == algebraic.asreal(1))
        var corrupted = output[0].copy()
        corrupted.factors[0][0] = corrupted.factors[0][0].reorder(MonomialOrder(), budget)
        assert_false(implementation._check_sparse_evidence(input, corrupted, output[1], budget))
    # KERNEL order is not multiplicative: the scalar is NOT generally the
    # greatest stored coefficient of the original polynomial.
    var x = sparse.to_core(sparse.asarray((((1, 0), algebraic.asreal(1)),)))
    var line = sparse.to_core(sparse.asarray((((1, 0), algebraic.asreal(1)), ((0, 1), 2))))
    var kernel = x.mul(line.pow(2, budget), budget).reorder(MonomialOrder(), budget)
    var result = implementation._factor_sparse_evidence(kernel, budget)
    assert_true(implementation._check_sparse_evidence(kernel, result[0], result[1], budget))
    assert_true(check_multivariate_factorization(kernel, result[0], budget))
    assert_true(result[0].unit == algebraic.asreal(4))
    assert_equal(len(result[0].factors), 2)
    var constant = sparse.to_core(sparse.asarray((((0, 0), algebraic.asreal(-3)),)))
    var cr = implementation._factor_sparse_evidence(constant, budget)
    assert_true(implementation._check_sparse_evidence(constant, cr[0], cr[1], budget))
    assert_equal(len(cr[0].factors), 0)
    var zero = constant.sub(constant, budget)
    var zr = implementation._factor_sparse_evidence(zero, budget)
    assert_true(implementation._check_sparse_evidence(zero, zr[0], zr[1], budget))
    assert_true(zr[0].unit.is_zero())
    print("array multivariate algebraic order checks passed")
