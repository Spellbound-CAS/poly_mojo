"""Multivariate field certificates and independently justified irreducibility cases."""
from poly import rational, dense, algebraic
from poly._rational_certificates import _Irreducibility, _check_irreducible
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, QQ, DensePoly, NoContext
from std.testing import assert_true, assert_false, assert_equal


def _integer_case(values: List[Int], expected: Bool, mut budget: Budget) raises:
    var coefficients = List[QQ]()
    for value in values:
        coefficients.append(QQ(Rational(value)))
    var p = DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)
    assert_equal(_check_irreducible(p, _Irreducibility(0, 0), budget), expected)


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var aa = rational.from_terms([((1, 0), a.copy())], [((0, 1), 1)])
    assert_true(rational.check(aa))
    var qi = rational.from_terms([((1, 0), ExactComplex(0, 1))], [((0, 1), 1)])
    assert_true(rational.check(qi))
    var z = algebraic.ascomplex(a, 1)
    var cc = rational.from_terms([((1, 0), z.copy())], [((0, 1), 1)])
    assert_true(rational.check(cc))
    var out = rational.to_core(cc)
    out[0].terms[0].coefficient = algebraic.ascomplex(100)
    assert_true(rational.check(cc))
    assert_true(rational.evaluate(cc, (1, 1)) == z)
    var constant = rational.from_core[ndim=0](rational.to_core(rational.asarray([2], [3])))
    assert_equal(constant.ndim(), 0)
    assert_true(rational.check(constant))

    var budget = Budget()
    _integer_case([1, 0, 0, 0, 1], True, budget)   # cyclotomic Phi_8
    _integer_case([3, 0, 0, 0, 1], True, budget)   # Eisenstein at 3
    _integer_case([4, 0, 5, 0, 1], False, budget)  # (x^2+1)(x^2+4)
    _integer_case([1, 0, 2, 0, 1], False, budget)  # (x^2+1)^2
    _integer_case([-2, 0, 1], True, budget)       # no rational square root of 2
    _integer_case([-6, 1, 1], False, budget)      # (x-2)(x+3)
    _integer_case([1], False, budget)            # a unit is not irreducible
    var fractional = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1, 4)), QQ(Rational(0)), QQ(Rational(1))], budget)
    assert_true(_check_irreducible(fractional, _Irreducibility(0, 0), budget))
    assert_false(_check_irreducible(fractional, _Irreducibility(2, 0), budget))
    assert_false(_check_irreducible(fractional, _Irreducibility(3, 1), budget))
    print("array rational certificate field checks passed")
