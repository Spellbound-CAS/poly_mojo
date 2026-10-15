"""Failure categories, whole-call budgets, and independent identity rejection."""
from poly import rational, dense, algebraic, coefficients, sparse
from apn_mojo import ExactComplex
from poly._common import Values
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_true, assert_false, assert_equal


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var input = coefficients((a.copy(), 1))
    var budget = Budget()
    var f = rational.asarray(input, [1], budget)
    var used = UInt64.MAX - budget.remaining()
    assert_true(rational.check(rational.asarray(input, [1], max_steps=used)))
    var rejected = False
    try:
        _ = rational.asarray(input, [1], max_steps=used - 1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    rejected = False
    try:
        _ = rational.asarray(input, List[Int]())
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    var zero = rational.asarray(List[algebraic.Real](), [1])
    assert_true(rational.check(zero))
    assert_equal(rational.denominator(zero).to_list(), [(0, algebraic.asreal(1))])
    rejected = False
    try:
        _ = rational.reciprocal(zero)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    var pole = rational.asarray([1], coefficients((algebraic.negative(a), 1)))
    rejected = False
    try:
        _ = rational.evaluate(pole, a)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    rejected = False
    try:
        _ = rational.laurent(pole, a, terms=-1)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    var local = rational.laurent(pole, a, terms=2)
    local.coefficients = coefficients((algebraic.asreal(2), 0))
    assert_false(rational.check(local))
    var parts = rational.partial_fractions(pole)
    parts.terms._values[][0][2] = Int.MAX
    assert_false(rational.check(parts))
    budget = Budget()
    assert_true(rational.check(f, budget))
    used = UInt64.MAX - budget.remaining()
    assert_true(rational.check(f, max_steps=used))
    rejected = False
    try:
        _ = rational.check(f, max_steps=used - 1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    # A changed output or cancellation cofactor must not retain a valid identity.
    var corrupt = rational.asarray(input, [1])
    corrupt._data[].value.numerator.terms[0].coefficient = algebraic.asreal(7)
    assert_false(rational.check(corrupt))
    var gaussian = rational.asarray([ExactComplex(1)], [ExactComplex(1), ExactComplex(0), ExactComplex(1)])
    assert_true(rational.check(rational.partial_fractions(gaussian)))
    assert_equal(len(rational.poles(gaussian)[1]), 2)
    assert_true(rational.check(rational.laurent(gaussian, ExactComplex(0, 1), terms=2)))
    print("array algebraic rational limit checks passed")
