from poly import dense, sparse, roots, algebraic
from std.testing import assert_equal, assert_true


def main() raises:
    # (x^2 + 1)^2 has the distinct roots -i and i, each twice.
    var result = dense.complex_roots([1, 0, 2, 0, 1])
    assert_equal(len(result), 2)
    assert_true(result[0] == algebraic.ascomplex(0, -1))
    assert_true(result[1] == algebraic.ascomplex(0, 1))
    assert_equal(result.multiplicity(0), 2)
    assert_equal(result.multiplicity(1), 2)
    assert_true(roots.check(result))
    for root in result:
        assert_true(dense.evaluate([1, 0, 2, 0, 1], root).is_zero())
    var sparse_result = sparse.complex_roots([(0, 1), (2, 1)])
    assert_true(sparse_result[0] == result[0])
    assert_true(roots.check_complex([1, 0, 2, 0, 1], result.evidence()))
    print("Exact roots: -i and i, each with multiplicity 2")
    print("Complete root certificate verified")
