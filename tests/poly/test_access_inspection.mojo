"""Sparse support and leading terms use chosen orders without changing stored arrays."""
from poly import sparse, coefficients, algebraic, modint
from apn_mojo import Integer, Ball, BallContext
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var p = sparse.asarray([((2, 0, 0), 3), ((0, 3, 0), 5), ((0, 0, 0), 7)])
    assert_equal(sparse.leading_term(p)[0], (2, 0, 0))
    assert_equal(sparse.leading_term(p, order="grlex")[0], (0, 3, 0))
    assert_equal(sparse.leading_term(p, order="lex", axes=List[Int]([1, 0, 2]))[1], Integer(5))
    assert_equal(sparse.variables(p), List[Int]([0, 1]))
    var q = sparse.asarray([((1, 0, 1), 2), ((0, 2, 0), 4)])
    assert_equal(sparse.leading_term(q, order="grlex")[0], (1, 0, 1))
    assert_equal(sparse.leading_term(q, order="grevlex")[0], (0, 2, 0))
    assert_equal(sparse.variables([((0, 0), 3)]), List[Int]())
    assert_equal(sparse.variables(List[Tuple[Int, Int]]()), List[Int]())
    assert_equal(sparse.variables([(2, 1), (2, -1)]), List[Int]())
    assert_equal(sparse.leading_term([(1, modint(3, 7))])[1].value, UInt64(3))
    var a = sparse.asarray([(1, algebraic.asreal(2))])
    assert_equal(algebraic.compare(sparse.leading_term(a)[1], 2), 0)
    var bc = BallContext(precision=50)
    var uncertain = sparse.asarray([((2, 0), Ball(0, 1, precision=50)), ((0, 1), Ball(2, precision=50))], context=bc)
    assert_equal(sparse.leading_term(uncertain)[0], (2, 0))
    assert_equal(sparse.variables(uncertain), List[Int]([0, 1]))
    for mode in range(4):
        var rejected = False
        try:
            if mode == 0: _ = sparse.leading_term(List[Tuple[Int, Int]]())
            elif mode == 1: _ = sparse.leading_term(p, order="bad")
            elif mode == 2: _ = sparse.leading_term(p, axes=List[Int]([0, 0, 2]))
            else: _ = sparse.leading_term(p, axes=List[Int]([1]))
        except e: rejected = e.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    var budget = Budget()
    _ = sparse.leading_term(p, budget, order="grevlex")
    var used = UInt64.MAX - budget.remaining()
    var exact = Budget(used)
    _ = sparse.leading_term(p, exact, order="grevlex")
    assert_equal(exact.remaining(), UInt64(0))
    var limited = False
    try: _ = sparse.leading_term(p, order="grevlex", max_steps=used - 1)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(p[0][0], (0, 0, 0))
    print("access sparse inspection checks passed")
