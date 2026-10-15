"""Reject false multivariate irreducibility and corrupted retained image proofs."""
from poly import sparse, dense, algebraic, coefficients, exact
from std.testing import assert_equal, assert_true, assert_false


def main() raises:
    var p = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 2), 1)))
    var r = sparse.factor_checked(p)
    assert_true(exact.check(r))
    var good = r.evidence().copy()
    var images = r._data[].algebraic_images.copy()
    assert_equal(len(images), 1)
    assert_true(Bool(images[0]))
    r._data[].algebraic_images.clear()
    assert_false(exact.check(r))
    r._data[].algebraic_images = images.copy()
    r._data[].algebraic_images[0] = None
    assert_false(exact.check(r))
    r._data[].algebraic_images = images.copy()
    r._data[].algebraic_images[0].value().unit = algebraic.asreal(2)
    assert_false(exact.check(r))
    r._data[].algebraic_images = images.copy()
    r._data[].algebraic_images[0].value().factors[0][1] = Int.MAX
    assert_false(exact.check(r, max_steps=1000000))
    r._data[].algebraic_images = images.copy()
    r._data[].proof.factors[0][1] = Int.MAX
    assert_false(exact.check(r, max_steps=1000000))
    r._data[].proof = good.copy()
    r._data[].proof.factors[0][0] = sparse.to_core(sparse.asarray((((1000000000, 0), algebraic.asreal(1)), ((0, 1), 1))))
    assert_false(exact.check(r, max_steps=1000000))
    r._data[].proof = good.copy()
    r._data[].proof.factors[0][0] = sparse.to_core(sparse.asarray((((2, 0), algebraic.asreal(2)), ((0, 2), 2))))
    r._data[].proof.unit = algebraic.divide(1, 2)
    assert_false(exact.check(r))
    r._data[].proof = good.copy()
    # The product and the image factorization are BOTH correct, but the claimed
    # multivariate factor is reducible. Only the independent search rejects it.
    var reducible = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 2), -1)))
    r._data[].input = sparse.to_core(reducible)
    r._data[].proof.factors = [(sparse.to_core(reducible), 1)]
    var image_result = dense.factor_checked(coefficients((algebraic.asreal(0), 0, 1, 0, 0, 0, -1)))
    r._data[].algebraic_images[0] = image_result.evidence().copy()
    assert_false(exact.check(r))
    r._data[].input = sparse.to_core(p)
    r._data[].proof = good.copy()
    r._data[].algebraic_images = images.copy()
    assert_true(exact.check(r))
    # A simple factor needs no image proof; arbitrary extra data is rejected.
    var linear = sparse.asarray((((1, 0), algebraic.asreal(1)), ((0, 1), 1)))
    var simple = sparse.factor_checked(linear)
    simple._data[].algebraic_images[0] = images[0].copy()
    assert_false(exact.check(simple))
    var repeated = sparse.factor_checked(sparse.power(linear, 2))
    repeated._data[].proof.factors = [(sparse.to_core(linear), 1), (sparse.to_core(linear), 1)]
    repeated._data[].algebraic_images = [None, None]
    assert_false(exact.check(repeated))
    var two = sparse.factor_checked(reducible)
    two._data[].proof.factors.swap_elements(0, 1)
    two._data[].algebraic_images.swap_elements(0, 1)
    assert_false(exact.check(two))
    var detached = r.evidence().copy()
    detached.factors.clear()
    assert_equal(len(exact.factors(r)), 1)
    assert_true(exact.check(r))
    print("array multivariate algebraic evidence checks passed")
