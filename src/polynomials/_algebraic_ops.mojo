"""Polynomial transformations and rational bounds shared by algebraic proofs."""
# Source: Z3 2d08bfc, algebraic_numbers.cpp (sum/product transforms).
# Z3
# Copyright (c) Microsoft Corporation
# All rights reserved.
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .sparse import SparsePoly, Term
from .monomial import Monomial, MonomialOrder
from ._dense_ops import _copy as _copy_poly
from ._root_ops import _compare as _compare_q, _copy as _copy_q
from .algebraic_types import RealAlgebraic as _Real


def _embed(f: DensePoly[QQ], variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    var moved = _copy_poly(f, budget)
    moved.variable = variable
    return moved.to_sparse(MonomialOrder(), budget)


def _compose(f: DensePoly[QQ], value: SparsePoly[QQ], mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    var result = SparsePoly[QQ].zero(value.order, NoContext(), budget)
    for i in range(len(f.coefficients) - 1, -1, -1):
        budget.charge(1)
        result = result.mul(value, budget).add(SparsePoly[QQ].constant(f.coefficients[i], value.order, budget), budget)
    return result^


def _product_polynomial(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    # Homogenize f(z/t) as t^degree(f) * f(z/t).
    var n = f.degree()
    if UInt64(max(0, n)) > UInt64(UInt32.MAX):
        raise PolyError(PolyErrorKind.CAPACITY, "algebraic resultant degree exceeds UInt32")
    _allocation[Term[QQ]](len(f.coefficients), budget)
    var terms = List[Term[QQ]]()
    for i in range(len(f.coefficients)):
        budget.charge(1)
        if f.coefficients[i].is_zero():
            continue
        _allocation[Tuple[Int64, UInt32]](2, budget)
        var powers = List[Tuple[Int64, UInt32]]()
        if i:
            powers.append((Int64(0), UInt32(i)))
        if n - i:
            powers.append((Int64(1), UInt32(n - i)))
        terms.append(Term[QQ](Monomial(powers^, budget), _copy_q(f.coefficients[i], budget)))
    return SparsePoly[QQ].from_terms(terms^, MonomialOrder(), NoContext(), budget)


def _product_bounds(a: _Real, b: _Real, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    var low = a.low.mul(b.low, budget)
    var high = _copy_q(low, budget)
    _allocation[QQ](3, budget)
    var products = [a.low.mul(b.high, budget), a.high.mul(b.low, budget), a.high.mul(b.high, budget)]
    for product in products:
        budget.charge(1)
        if _compare_q(product, low, budget) < 0:
            low = _copy_q(product, budget)
        if _compare_q(product, high, budget) > 0:
            high = _copy_q(product, budget)
    return (low^, high^)


def _binary_definitions(a: _Real, b: _Real, operation: Int, mut budget: Budget) raises PolyError -> Tuple[SparsePoly[QQ], SparsePoly[QQ]]:
    var first: SparsePoly[QQ]
    if operation == 2:
        first = _product_polynomial(a.polynomial, budget)
    else:
        var z = SparsePoly[QQ].variable(0, MonomialOrder(), NoContext(), budget)
        var t = SparsePoly[QQ].variable(1, MonomialOrder(), NoContext(), budget)
        first = _compose(a.polynomial, z.sub(t, budget) if operation == 0 else z.add(t, budget), budget)
    return (first^, _embed(b.polynomial, 1, budget))


def _unary_definition(a: _Real, operation: Int, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    _allocation[QQ](len(a.polynomial.coefficients), budget)
    var coefficients = List[QQ]()
    for i in range(len(a.polynomial.coefficients)):
        budget.charge(1)
        if operation == 3:
            coefficients.append(a.polynomial.coefficients[i].neg(budget) if i % 2 else _copy_q(a.polynomial.coefficients[i], budget))
        else:
            coefficients.append(_copy_q(a.polynomial.coefficients[a.polynomial.degree() - i], budget))
    return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)


def _operation_bounds(a: _Real, b: _Real, operation: Int, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    if operation == 0:
        return (a.low.add(b.low, budget), a.high.add(b.high, budget))
    if operation == 1:
        return (a.low.sub(b.high, budget), a.high.sub(b.low, budget))
    if operation == 2:
        return _product_bounds(a, b, budget)
    if operation == 3:
        return (a.high.neg(budget), a.low.neg(budget))
    return (a.high.inverse(budget), a.low.inverse(budget))


def _rational_operation(a: QQ, b: QQ, operation: Int, mut budget: Budget) raises PolyError -> QQ:
    if operation == 0:
        return a.add(b, budget)
    if operation == 1:
        return a.sub(b, budget)
    if operation == 2:
        return a.mul(b, budget)
    if operation == 3:
        return a.neg(budget)
    return a.inverse(budget)


def _is_zero_real(a: _Real, mut budget: Budget) raises PolyError -> Bool:
    budget.charge(1)
    if a.rational:
        return a.low.is_zero()
    return a.polynomial.coefficients[0].is_zero() and a.low.value.sign() < 0 and a.high.value.sign() > 0


def _norm_expression(coefficients: List[_Real], mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    var count = len(coefficients)
    while count and _is_zero_real(coefficients[count - 1], budget):
        budget.charge(1)
        count -= 1
    var order = MonomialOrder()
    var z = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var result = SparsePoly[QQ].zero(order, NoContext(), budget)
    for i in range(count - 1, -1, -1):
        budget.charge(1)
        var coefficient = SparsePoly[QQ].constant(coefficients[i].low, order, budget) if coefficients[i].rational else SparsePoly[QQ].variable(Int64(_sum(i, 1)), order, NoContext(), budget)
        result = result.mul(z, budget).add(coefficient, budget)
    return result^


def _cauchy_bound(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> QQ:
    var bound = QQ.one(NoContext())
    if f.degree() > 0:
        var inverse = f.coefficients[f.degree()].inverse(budget)
        for i in range(f.degree()):
            budget.charge(1)
            var ratio = f.coefficients[i].mul(inverse, budget)
            if ratio.value.sign() < 0:
                ratio = ratio.neg(budget)
            bound = bound.add(ratio, budget)
    return bound^
