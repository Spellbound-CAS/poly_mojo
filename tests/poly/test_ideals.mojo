"""Reusable bases, independent checks, and original-generator identities."""
from poly import sparse, modint
from polynomials import Budget, PolyErrorKind, SparsePoly, QQ, NoContext
from apn_mojo import Rational
from std.memory import ArcPointer
from std.testing import assert_equal, assert_true


def main() raises:
    var x_minus_y = [((1, 0), 1), ((0, 1), -1)]
    var x_squared_minus_one = [((2, 0), 1), ((0, 0), -1)]
    var generators = [x_minus_y.copy(), x_squared_minus_one.copy()]
    var basis = sparse.groebner(generators, order="lex")
    assert_equal(len(basis), 2)
    assert_equal(basis[0].domain(), "QQ")
    assert_true(sparse.check_groebner(basis))
    var member = [((0, 2), 1), ((0, 0), -1)]
    assert_true(sparse.contains(member, basis))
    assert_equal(len(sparse.reduce(member, basis)), 0)
    assert_true(not sparse.contains([((0, 0), 1)], basis))
    var identity = sparse.membership(member, basis)
    assert_true(Bool(identity))
    assert_true(sparse.check_membership(identity.value()))
    assert_equal(len(identity.value().cofactors), 2)
    var recovered = sparse.add(sparse.multiply(identity.value().cofactors[0], generators[0]), sparse.multiply(identity.value().cofactors[1], generators[1]))
    assert_equal(len(sparse.subtract(recovered, member)), 0)
    assert_true(not sparse.membership([((0, 0), 1)], basis))
    var reversed_basis = sparse.groebner(generators, order="lex", axes=[1, 0])
    assert_true(sparse.check_groebner(reversed_basis))
    assert_true(sparse.contains(member, reversed_basis))
    for p in basis:
        assert_true(sparse.contains(p, reversed_basis))
    var lean = sparse.groebner(generators, witnesses=False)
    assert_true(sparse.contains(member, lean))
    var failed = False
    try:
        _ = sparse.membership(member, lean)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var finite_a = sparse.asarray(x_minus_y, modulus=5)
    var finite_b = sparse.asarray(x_squared_minus_one, modulus=5)
    var finite = sparse.groebner([finite_a.copy(), finite_b.copy()])
    assert_true(sparse.contains(member, finite))
    assert_true(sparse.check_groebner(finite))
    var empty = List[List[Tuple[Tuple[Int, Int], Int]]]()
    var zero = sparse.groebner(empty)
    assert_equal(len(zero), 0)
    assert_true(sparse.check_groebner(zero))
    assert_true(not sparse.contains(member, zero))
    # Original-generator identities and containment alone are insufficient:
    # x**2-y and x*y-1 are not a Groebner basis in lexicographic order.
    var altered = basis._data[].copy()
    var first_bad = sparse.to_core([((2, 0), Rational(1)), ((0, 1), Rational(-1))])
    var second_bad = sparse.to_core([((1, 1), Rational(1)), ((0, 0), Rational(-1))])
    altered.generators = [first_bad.copy(), second_bad.copy()]
    var first_visible = sparse.from_core[ndim=2](first_bad)
    var second_visible = sparse.from_core[ndim=2](second_bad)
    altered.basis.basis = [first_bad^, second_bad^]
    var work = Budget()
    var one = SparsePoly[QQ].constant(QQ(Rational(1)), altered.basis.order, work)
    var none = SparsePoly[QQ].zero(altered.basis.order, NoContext(), work)
    var rows: List[List[SparsePoly[QQ]]] = [[one.copy(), none.copy()], [none.copy(), one.copy()]]
    altered.basis.cofactors = rows^
    var visible = basis._values.copy()
    var entries: List[sparse.Array[Tuple[Int, Int], 1]] = [first_visible^, second_visible^]
    visible._values = ArcPointer(entries^)
    var bogus = sparse.Basis[Tuple[Int, Int], 1](ArcPointer(altered^), visible^)
    assert_true(not sparse.check_groebner(bogus))
    assert_true(sparse.check_groebner(basis))
    print("array ideal checks passed")
