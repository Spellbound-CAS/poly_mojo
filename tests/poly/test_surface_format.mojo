"""Hand-written display expectations, name validation, and shared work limits."""
from poly import dense, sparse, modint
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var p = dense.asarray([1, -2, 3])
    assert_equal(dense.format(p, variable="t"), "(3)*t^2 + (-2)*t + (1)")
    assert_equal(dense.format(List[Int]()), "0")
    assert_equal(dense.format([Rational(1, 3)]), "(1/3)")
    assert_equal(dense.format([modint(9, 7), modint(1, 7)]), "(1)*x + (2)")
    assert_equal(dense.format([ExactComplex(Rational(1), Rational(-2))]), "(1 + (-2)*i)")
    var s = sparse.asarray([((2, 1, 0), 3), ((0, 0, 0), -1)])
    assert_equal(sparse.format(s, variables=["u", "v", "unused"]), "(3)*u^2*v + (-1)")
    assert_equal(sparse.format(s), "(3)*x0^2*x1 + (-1)")
    assert_equal(sparse.format([(2, 1)]), "(1)*x^2")
    for mode in range(6):
        var rejected = False
        try:
            if mode == 0: _ = dense.format(p, variable="2x")
            elif mode == 1: _ = dense.format(List[Int](), variable="")
            elif mode == 2: _ = sparse.format(s, variables=["x", "x", "z"])
            elif mode == 3: _ = sparse.format(s, variables=["x", "y"])
            elif mode == 4: _ = dense.format([ExactComplex(Rational(1), Rational(2))], variable="i")
            else: _ = dense.format(p, variable="x+y")
        except e: rejected = e.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    var meter = Budget()
    _ = dense.format(p, meter)
    var used = UInt64.MAX - meter.remaining()
    var exact = Budget(used)
    _ = dense.format(p, exact)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(used - 1)
    var limited = False
    try: _ = dense.format(p, short)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(p.to_list(), dense.asarray([1, -2, 3]).to_list())
    print("surface named formatting checks passed")
