"""Whole-call certificate budgets, retries, and preserved borrowed inputs."""
from poly import rational, dense
from poly._rational_certificates import _Irreducibility, _check_irreducible
from apn_mojo import Rational, ExactComplex, Integer
from polynomials import Budget, PolyErrorKind, DensePoly, QQ, NoContext
from std.testing import assert_true, assert_equal


def _checker[P: Copyable & Deinitable, k: Int](value: rational.Function[P, k]) raises:
    var budget = Budget()
    assert_true(rational.check(value, budget))
    var used = UInt64.MAX - budget.remaining()
    assert_true(used > 0)
    assert_true(rational.check(value, max_steps=used))
    var limited = False
    try:
        _ = rational.check(value, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(rational.check(value))


def _checker[P: Copyable & Deinitable, k: Int](value: rational.PartialFractions[P, k]) raises:
    var budget = Budget()
    assert_true(rational.check(value, budget))
    var used = UInt64.MAX - budget.remaining()
    assert_true(used > 0)
    assert_true(rational.check(value, max_steps=used))
    var limited = False
    try:
        _ = rational.check(value, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(rational.check(value))


def _checker[P: Copyable & Deinitable, k: Int](value: rational.Expansion[P, k]) raises:
    var budget = Budget()
    assert_true(rational.check(value, budget))
    var used = UInt64.MAX - budget.remaining()
    assert_true(used > 0)
    assert_true(rational.check(value, max_steps=used))
    var limited = False
    try:
        _ = rational.check(value, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(rational.check(value))


def main() raises:
    _ = rational.asarray([0, 1], [0, 2, 2])
    var numerator = dense.asarray([0, 1])
    var denominator = dense.asarray([0, 2, 2])
    var budget = Budget()
    var f = rational.asarray(numerator, denominator, budget)
    var used = UInt64.MAX - budget.remaining()
    assert_true(rational.check(rational.asarray(numerator, denominator, max_steps=used)))
    var limited = False
    try:
        _ = rational.asarray(numerator, denominator, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(numerator.to_list(), [Integer(0), Integer(1)])
    assert_equal(denominator.to_list(), [Integer(0), Integer(2), Integer(2)])
    _checker(f)
    _checker(rational.from_terms([((1, 0), 1)], [((0, 1), 1)]))

    var source = rational.asarray([1], [1, 0, 0, 0, 1])
    _ = rational.partial_fractions(source)
    budget = Budget()
    var parts = rational.partial_fractions(source, budget)
    used = UInt64.MAX - budget.remaining()
    assert_true(rational.check(rational.partial_fractions(source, max_steps=used)))
    limited = False
    try:
        _ = rational.partial_fractions(source, max_steps=used - 1)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    _checker(parts)
    _checker(rational.partial_fractions(rational.asarray([1], [-2, 0, 0, 1])))
    _checker(rational.partial_fractions(rational.asarray([ExactComplex(1)], [ExactComplex(-2), ExactComplex(0), ExactComplex(1)])))
    _checker(rational.laurent(source, 0, terms=4))
    assert_equal(rational.evaluate(source, 1), Rational(1, 2))

    # The finite but potentially expensive fallback must stop at the supplied
    # budget without retaining a partial answer or mutating its polynomial.
    budget = Budget()
    var big = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational((Integer(1) << 512) + 1)), QQ(Rational(0)), QQ(Rational(0)), QQ(Rational(0)), QQ(Rational(1))], budget)
    var before = big.copy()
    var small = Budget(1000)
    limited = False
    try:
        _ = _check_irreducible(big, _Irreducibility(0, 0), small)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(big == before)
    print("array rational certificate budget checks passed")
