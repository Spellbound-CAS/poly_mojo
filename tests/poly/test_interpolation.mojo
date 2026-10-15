"""Exact sample constraints, reusable plans, and Taylor transformations."""
import poly
from apn_mojo import Integer, Rational, ExactComplex
from polynomials import PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    var p = poly.dense.interpolate([0, 1, 2], [1, 4, 9])
    assert_equal(p.to_list(), poly.coefficients((Rational(1), 2, 1)).to_list())
    var zero = poly.dense.interpolate([0, 1, 2], [0, 0, 0])
    assert_equal(len(zero), 0)
    var plan = poly.dense.prepare_points([0, 1, 2], interpolation=True)
    assert_equal(plan.domain(), "QQ")
    assert_equal(plan.point_count(), 3)
    assert_equal(plan.evaluate(p).to_list(), poly.coefficients((Rational(1), 4, 9)).to_list())
    assert_equal(plan.interpolate([1, 4, 9]).to_list(), p.to_list())
    var copied = plan.copy()
    assert_equal(copied.evaluate([0, 1]).to_list(), poly.coefficients((Rational(0), 1, 2)).to_list())
    var hermite = poly.dense.hermite_interpolate([0, 1], [[1, 2], [4, 4]])
    assert_equal(hermite.to_list(), p.to_list())
    assert_equal(poly.dense.hermite_interpolate([0], [[1, 2, 1]], derivatives=False).to_list(), p.to_list())
    assert_equal(poly.dense.affine_transform([1, 2, 1], 2, 1).to_list(), poly.dense.asarray([4, 8, 4]).to_list())
    assert_equal(poly.dense.taylor_shift([1, 2, 1], 1).to_list(), poly.dense.asarray([4, 4, 1]).to_list())
    assert_equal(poly.dense.jet_evaluate([1, 2, 1], 1, 5).to_list(), poly.coefficients((4, 4, 1, 0, 0)).to_list())
    assert_equal(poly.dense.definite_integral([0, 1], 0, 1), Rational(1, 2))
    var finite = poly.coefficients((poly.modint(0, 5), 1, 2))
    var fp = poly.dense.prepare_points(finite, interpolation=True)
    var fi = fp.interpolate([1, 4, 9])
    assert_equal(fp.evaluate(fi).to_list(), poly.coefficients((poly.modint(1, 5), 4, 4)).to_list())
    var direct = poly.dense.interpolate(finite, [1, 4, 9])
    assert_equal(direct.to_list(), fi.to_list())
    var complex_nodes = poly.coefficients((0, ExactComplex(0, 1)))
    var cp = poly.dense.prepare_points(complex_nodes)
    assert_equal(cp.interpolate([0, 1]).domain(), "QQi")
    var repeated = poly.dense.prepare_points([1, 1])
    assert_equal(repeated.evaluate([0, 1]).to_list(), poly.coefficients((Rational(1), 1)).to_list())
    var failed = False
    try:
        _ = repeated.interpolate([1, 1])
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    failed = False
    try:
        _ = plan.evaluate(p, max_steps=0)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    print("array interpolation checks passed")
