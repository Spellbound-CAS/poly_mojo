"""Prove an ideal membership using an independently checked identity."""
from std.testing import assert_equal, assert_true
from polynomials import Budget, QQ, NoContext, PolynomialRing, SparsePoly, groebner, membership, check_membership


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x", "y"], NoContext(), work)
    var generators: List[SparsePoly[QQ]] = [ring.parse("x - y", work), ring.parse("x^2 - 1", work)]
    var basis = groebner(generators, generators[0].order, work)
    var candidate = ring.parse("y^2 - 1", work)
    var witness = membership(candidate, generators, candidate.order, work)
    assert_true(Bool(witness))
    var verified = check_membership(candidate, generators, witness.value(), work)
    assert_equal(len(basis.basis), 2)
    assert_true(verified)
    print("Basis elements:", len(basis.basis))
    print("Membership checked:", verified)
