"""Reuse a basis and recover a witness against the original generators."""
from poly import sparse
from std.testing import assert_equal, assert_true


def main() raises:
    var first = [((1, 0), 1), ((0, 1), -1)]  # x - y
    var second = [((2, 0), 1), ((0, 0), -1)]  # x**2 - 1
    var basis = sparse.groebner([first.copy(), second.copy()], order="lex")
    assert_true(sparse.check_groebner(basis))
    var candidate = [((0, 2), 1), ((0, 0), -1)]  # y**2 - 1
    print("Basis polynomials:", len(basis))
    print("Member:", sparse.contains(candidate, basis))
    var witness = sparse.membership(candidate, basis)
    assert_true(Bool(witness))
    assert_true(sparse.check_membership(witness.value()))
    var recovered = sparse.add(sparse.multiply(witness.value().cofactors[0], first), sparse.multiply(witness.value().cofactors[1], second))
    assert_equal(len(sparse.subtract(recovered, candidate)), 0)
    print("Original-generator identity checked")
