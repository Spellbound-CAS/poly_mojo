"""Reuse exact sample nodes and change variables in a sparse polynomial."""
import poly
from std.testing import assert_equal


def main() raises:
    var plan = poly.dense.prepare_points([0, 1, 2], interpolation=True)
    var p = plan.interpolate([1, 4, 9])
    print("Recovered:", p.to_list())
    print("Samples:", plan.evaluate(p).to_list())
    var hermite = poly.dense.hermite_interpolate([0, 1], [[1, 2], [4, 4]])
    assert_equal(hermite.to_list(), p.to_list())
    print("Taylor coefficients:", poly.dense.jet_evaluate(p, 1, 4).to_list())
    print("Integral:", poly.dense.definite_integral([0, 1], 0, 1))
    var surface = poly.sparse.asarray([((2, 1), 1), ((0, 1), 1)])
    var x = poly.sparse.asarray([((1, 0), 1)])
    var y = poly.sparse.asarray([((0, 1), 1)])
    var swapped = poly.sparse.substitute(surface, [y.copy(), x.copy()])
    print("Swapped axes:", swapped.to_list())
    print("Integrated x:", poly.sparse.definite_integral(surface, 0, 1).to_list())
