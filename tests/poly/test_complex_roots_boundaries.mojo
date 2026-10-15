"""Reject incomplete or tampered exact-root evidence and preserve input on limits."""
from poly import dense, sparse, roots, algebraic
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true, assert_false


def main() raises:
    var p = dense.asarray([1, 0, 1])
    var result = dense.complex_roots(p)
    var proof = result.evidence().copy()
    _ = proof.roots.pop()
    _ = proof.multiplicities.pop()
    assert_false(roots.check_complex(p, proof))
    proof = result.evidence().copy()
    proof.multiplicities[0] = 0
    assert_false(roots.check_complex(p, proof))
    proof.multiplicities[0] = Int.MAX
    assert_false(roots.check_complex(p, proof))
    proof = result.evidence().copy()
    proof.roots[1] = proof.roots[0].copy()
    assert_false(roots.check_complex(p, proof))
    proof = result.evidence().copy()
    proof.roots[0] = algebraic.ascomplex(0, -2)
    assert_false(roots.check_complex(p, proof))
    proof = result.evidence().copy()
    var saved = proof.roots[0].copy()
    proof.roots[0] = proof.roots[1].copy()
    proof.roots[1] = saved^
    assert_false(roots.check_complex(p, proof))
    assert_false(roots.check_complex([2, 0, 1], result.evidence()))
    assert_false(roots.check_complex(List[Int](), result.evidence()))
    assert_true(roots.check(result))
    assert_equal(len(result), 2)
    for steps in [0, 1, 100, 1000]:
        var budget = Budget(UInt64(steps))
        var limited = False
        try:
            _ = dense.complex_roots(p, budget)
        except e:
            limited = e.kind == PolyErrorKind.LIMIT
        assert_true(limited)
        assert_equal(p[0], 1)
    var measured = Budget()
    _ = dense.complex_roots(p, measured)
    var used = UInt64.MAX - measured.remaining()
    var exact_budget = Budget(used)
    _ = dense.complex_roots(p, exact_budget)
    assert_equal(exact_budget.remaining(), UInt64(0))
    var short_budget = Budget(used - 1)
    var short_rejected = False
    try:
        _ = dense.complex_roots(p, short_budget)
    except e:
        short_rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(short_rejected)
    var check_limited = False
    try:
        _ = roots.check(result, max_steps=0)
    except e:
        check_limited = e.kind == PolyErrorKind.LIMIT
    assert_true(check_limited)
    var zero_rejected = False
    try:
        _ = dense.complex_roots([0, 0])
    except e:
        zero_rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(zero_rejected)
    var shape_rejected = False
    try:
        _ = sparse.complex_roots([((1, 0), 1), ((0, 1), 1)])
    except e:
        shape_rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(shape_rejected)
    var axis_rejected = False
    try:
        _ = sparse.complex_roots([(2, 1), (0, 1)], axis=1)
    except e:
        axis_rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(axis_rejected)
    print("array complex root boundary checks passed")
