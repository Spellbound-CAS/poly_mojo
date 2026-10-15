"""Incremental root and CAD checking uses the coefficient-array boundary."""
from poly import dense, sparse, cad
from polynomials import Budget, NumericOptions
from std.testing import assert_true


def main() raises:
    var budget = Budget()
    var options = NumericOptions(absolute_bits=2, max_subdivisions=256)
    var real_checker = dense.real_root_checker([-2, 0, 1], budget, options=options)
    var receipt = dense.real_root_cover_to([-2, 0, 1], real_checker, budget, options=options, complete_only=True)
    assert_true(real_checker.finish(receipt, budget))
    var complex_checker = sparse.complex_root_checker([(0, 1), (2, 1)], budget, options=options)
    var complex_receipt = sparse.complex_root_cover_to([(0, 1), (2, 1)], complex_checker, budget, options=options, complete_only=True)
    assert_true(complex_checker.finish(complex_receipt, budget))
    var real_sparse_checker = sparse.real_root_checker([(0, -2), (2, 1)], budget, options=options)
    var sparse_receipt = sparse.real_root_cover_to([(0, -2), (2, 1)], real_sparse_checker, budget, options=options)
    assert_true(real_sparse_checker.finish(sparse_receipt, budget))
    var direct_complex_checker = dense.complex_root_checker([1, 0, 1], budget, options=options)
    var direct_receipt = dense.complex_root_cover_to([1, 0, 1], direct_complex_checker, budget, options=options)
    assert_true(direct_complex_checker.finish(direct_receipt, budget))
    var inputs = [[(0, -1), (2, 1)]]
    var checker = cad.checker(inputs, budget)
    var cad_receipt = cad.decompose_to(inputs, checker, budget)
    assert_true(checker.finish(cad_receipt, budget))
    print("array solver stream checks passed")
