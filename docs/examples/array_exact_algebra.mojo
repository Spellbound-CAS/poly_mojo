"""Recover a common factor and factor a polynomial from coefficient inputs."""
import poly
from std.testing import assert_equal


def main() raises:
    var p = poly.dense.multiply(poly.dense.power([-1, 1], 2), [2, 0, 2])
    var common = poly.dense.gcd(p, [-2, 2])
    print("GCD:", common.to_list())
    var bezout = poly.dense.extended_gcd(p, [-2, 2])
    var identity = poly.dense.add(poly.dense.multiply(bezout[1], p), poly.dense.multiply(bezout[2], [-2, 2]))
    assert_equal(identity.to_list(), bezout[0].to_list())
    print("Bezout identity checked")
    var factors = poly.dense.factor(p)
    print("Scalar:", factors[0], "Distinct factors:", len(factors[1]))
    var rebuilt = poly.dense.asarray([factors[0].copy()])
    for entry in factors[1]:
        rebuilt = poly.dense.multiply(rebuilt, poly.dense.power(entry[0], entry[1]))
    assert_equal(rebuilt.to_list(), p.to_list())
    print("Factor product checked")
    var division = poly.dense.pseudo_divmod([1, 0, 1], [1, 2])
    print("Scale:", division[0], "Quotient:", division[1].to_list(), "Remainder:", division[2].to_list())
