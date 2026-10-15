"""Use core certificates and variable IDs without parsing an expression."""
import poly
from polynomials import Budget, factor_integer, check_factorization
from std.testing import assert_equal, assert_true


def main() raises:
    var work = Budget(1000000)
    var p = poly.dense.asarray([-1, 0, 1], work)
    var core = poly.dense.to_core(p, work, variable=7)
    var result = factor_integer(core, work)
    assert_true(check_factorization(core, result, work))
    print("Core factor product checked")
    var first = poly.dense.from_core(result.factors[0][0], work)
    print("First factor:", first.to_list(work))
    var surface = poly.sparse.asarray([((1, 0), 2), ((0, 1), 3)], work)
    var mapped = poly.sparse.to_core(surface, work, variables=[7, 2])
    var restored = poly.sparse.from_core[ndim=2](mapped, work, variables=[7, 2])
    assert_equal(restored.to_list(work), surface.to_list(work))
    print("Variable map round trip checked")
