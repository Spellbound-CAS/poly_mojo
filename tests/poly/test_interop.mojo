"""Core conversion ownership, variable maps, validation, and shared limits."""
import poly
from apn_mojo import Integer
from polynomials import Budget, PolyError, PolyErrorKind, ZZ, ModInt, Modulus, DensePoly, MonomialOrder
from std.testing import assert_equal, assert_true


def main() raises:
    var a = poly.dense.asarray([1, 2, 3])
    var core = poly.dense.to_core(a, variable=19)
    assert_equal(core.variable, Int64(19))
    var b = poly.dense.from_core(core)
    core.coefficients[0] = ZZ(Integer(99))
    assert_equal(a[0], Integer(1))
    assert_equal(b[0], Integer(1))
    var p = poly.sparse.asarray([((2, 0), 3), ((1, 1), 2), ((0, 1), 1)])
    var ids: List[Int64] = [7, 2]
    var s = poly.sparse.to_core(p, variables=ids)
    assert_equal(s.order.variables, ids)
    var q = poly.sparse.from_core[ndim=2](s, variables=ids)
    assert_equal(q.to_list(), p.to_list())
    s.terms[0].coefficient = ZZ(Integer(99))
    assert_equal(q.to_list(), p.to_list())
    var failed = False
    try:
        _ = poly.sparse.from_core[ndim=2](s)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    failed = False
    try:
        _ = poly.sparse.to_core(p, variables=[4, 4])
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    failed = False
    try:
        _ = poly.dense.to_core(a, max_steps=0)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var bad = DensePoly[ModInt](0, Modulus(_prime_value=9))
    failed = False
    try:
        _ = poly.dense.from_core(bad)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var zero = poly.sparse.to_core(List[Tuple[Tuple[Int, Int], Int]]())
    zero.order.kind = 99
    failed = False
    try:
        _ = poly.sparse.from_core[ndim=2](zero)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    print("array interoperability checks passed")
