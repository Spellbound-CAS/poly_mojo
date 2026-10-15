"""Shared incremental finite-field elimination for sparse interpolation."""
# Source-guided sparse interpolation: Z3 2d08bfc, polynomial.cpp.
# Copyright (c) Microsoft Corporation. All rights reserved.
# MIT License
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
from .errors import Budget, PolyError, _allocation, _sum
from .coefficients import ModInt, Modulus


def _insert_row(mut rows: List[List[ModInt]], var row: List[ModInt], mut budget: Budget) raises PolyError -> Int:
    # One normalized row per pivot; the last column is the right-hand side.
    var count = len(rows)
    for j in range(count):
        budget.charge(1)
        if row[j].is_zero():
            continue
        if len(rows[j]):
            var scale = row[j].copy()
            for k in range(j, _sum(count, 1)):
                budget.charge(1)
                row[k] = row[k].sub(scale.mul(rows[j][k], budget), budget)
        else:
            var inverse = row[j].inverse(budget)
            for k in range(j, _sum(count, 1)):
                budget.charge(1)
                row[k] = row[k].mul(inverse, budget)
            rows[j] = row^
            return 1
    return 0 if row[count].is_zero() else -1


def _solve_rows(rows: List[List[ModInt]], ctx: Modulus, mut budget: Budget) raises PolyError -> List[ModInt]:
    var count = len(rows)
    _allocation[ModInt](count, budget)
    var solution = List[ModInt](length=count, fill=ModInt.zero(ctx))
    for j in range(count - 1, -1, -1):
        budget.charge(1)
        solution[j] = rows[j][count].copy()
        for k in range(j + 1, count):
            budget.charge(1)
            solution[j] = solution[j].sub(rows[j][k].mul(solution[k], budget), budget)
    return solution^
