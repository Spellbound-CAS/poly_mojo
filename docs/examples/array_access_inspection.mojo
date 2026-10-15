"""Inspect sparse support without changing the polynomial's storage order."""
from poly import sparse
from apn_mojo import Integer
from std.testing import assert_equal


def main() raises:
    var p = sparse.asarray([((2, 0, 0), 3), ((0, 3, 0), 5), ((0, 0, 0), 7)])
    assert_equal(sparse.leading_term(p)[0], (2, 0, 0))
    assert_equal(sparse.leading_term(p, order="grlex")[0], (0, 3, 0))
    assert_equal(sparse.leading_term(p, axes=List[Int]([1, 0, 2]))[1], Integer(5))
    assert_equal(sparse.variables(p), List[Int]([0, 1]))
    print("Lex selects 3*x^2; graded lex selects 5*y^3")
    print("Occurring axes: 0, 1")
