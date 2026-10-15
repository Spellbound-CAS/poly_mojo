"""Subresultant PRS, with the Brown scaling used by SymPy euclidtools.

The recurrence follows dup_inner_subresultants/dmp_inner_subresultants in
SymPy 1.13.3, sympy/polys/euclidtools.py; cofactor propagation is added here.

Copyright (c) 2006-2023 SymPy Development Team

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of SymPy nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
"""
from .errors import Budget, PolyError, _allocation
from .coefficients import CoefficientRing
from .sparse import SparsePoly
from ._tower import _one, _zero, _lc, _pseudo, _exact_poly


@fieldwise_init
struct _PrsEntry[C: CoefficientRing](Copyable, Movable):
    var polynomial: SparsePoly[Self.C]
    var principal: SparsePoly[Self.C]
    var index: Int
    var u: SparsePoly[Self.C]
    var v: SparsePoly[Self.C]


def _prs[C: CoefficientRing](left: SparsePoly[C], right: SparsePoly[C], variable: Int64, carry: Bool, mut budget: Budget, history: Bool = True) raises PolyError -> List[_PrsEntry[C]]:
    left._same(right, budget)
    var f = left._copy(budget)
    var g = right._copy(budget)
    var n = f.degree_in(variable, budget)
    var m = g.degree_in(variable, budget)
    var swapped = n < m
    if swapped:
        swap(f, g)
        swap(n, m)
    var result = List[_PrsEntry[C]]()
    if g.is_zero():
        return result^
    var u0 = _one(f, budget)
    var v0 = _zero(f, budget)
    var u1 = _zero(f, budget)
    var v1 = _one(f, budget)
    if swapped:
        swap(u0, v0)
        swap(u1, v1)
    var d = n - m
    var sign = _one(f, budget)
    if (d + 1) % 2:
        sign = sign.neg(budget)
    var lc = _lc(g, variable, budget)
    var c = lc.pow(UInt64(d), budget).neg(budget)
    var divisor = sign^
    while not g.is_zero() and m > 0:
        budget.charge(1)
        var index = m - 1
        var qr = _pseudo(f, g, variable, carry, budget)
        var h = _exact_poly(qr[1], divisor, budget)
        if h.is_zero():
            break
        var u = _zero(f, budget)
        var v = _zero(f, budget)
        if carry:
            u = _exact_poly(u0.mul(qr[2], budget).sub(qr[0].mul(u1, budget), budget), divisor, budget)
            v = _exact_poly(v0.mul(qr[2], budget).sub(qr[0].mul(v1, budget), budget), divisor, budget)
        var k = h.degree_in(variable, budget)
        d = m - k
        divisor = lc.neg(budget).mul(c.pow(UInt64(d), budget), budget)
        lc = _lc(h, variable, budget)
        if d > 1:
            c = _exact_poly(lc.neg(budget).pow(UInt64(d), budget), c.pow(UInt64(d - 1), budget), budget)
        else:
            c = lc.neg(budget)
        var entry = _PrsEntry[C](h._copy(budget), c.neg(budget), index, u._copy(budget), v._copy(budget))
        if not history and len(result):
            _ = result.pop()
        _allocation[_PrsEntry[C]](1, budget)
        result.append(entry^)
        f = g^
        g = h^
        u0 = u1^
        v0 = v1^
        u1 = u^
        v1 = v^
        m = k
    return result^
