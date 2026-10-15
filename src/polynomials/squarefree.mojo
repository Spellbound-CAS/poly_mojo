"""Squarefree decomposition in characteristic zero and over prime fields.

Reference: SymPy sqfreetools.py (Yun recurrence).
Our loops use explicit worklists and a shared budget.

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
from .coefficients import CoefficientField, ModInt
from .dense import DensePoly
from ._coefficient_ops import _copy as _copy_coefficient
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from .certificates import SquarefreeDecomposition
from ._dense_ops import _monic, _derivative, _gcd, _exact


def squarefree[C: CoefficientField](f: DensePoly[C], mut budget: Budget) raises PolyError -> SquarefreeDecomposition[C]:
    """Decompose a polynomial into monic, pairwise coprime multiplicity blocks.

    Parameters:
        C: A characteristic-zero coefficient field or ModInt prime field.

    Args:
        f: The polynomial to decompose, including zero or constants.
        budget: The shared work allowance.

    Returns:
        The leading coefficient and nonconstant blocks in increasing multiplicity.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    comptime if C == ModInt:
        return rebind_var[SquarefreeDecomposition[C]](_finite_squarefree(rebind[DensePoly[ModInt]](f), budget))
    f.check_invariants(budget)
    var unit = C.zero(f.ring)
    if not f.is_zero():
        unit = _copy_coefficient(f.coefficients[f.degree()], budget)
    var factors = List[Tuple[DensePoly[C], Int]]()
    if f.degree() <= 0:
        return SquarefreeDecomposition[C](unit^, factors^)
    var monic = _monic(f, budget)
    var derivative = _derivative(monic, budget)
    var common = _gcd(monic, derivative, budget)
    var b = _exact(monic, common, budget)
    var c = _exact(derivative, common, budget)
    var multiplicity = 1
    while b.degree() > 0:
        budget.charge(1)
        var d = c.sub(_derivative(b, budget), budget)
        var a = _gcd(b, d, budget)
        b = _exact(b, a, budget)
        c = _exact(d, a, budget)
        if a.degree() > 0:
            _allocation[Tuple[DensePoly[C], Int]](1, budget)
            factors.append((a^, multiplicity))
        multiplicity = _sum(multiplicity, 1)
    return SquarefreeDecomposition[C](unit^, factors^)


def _finite_squarefree(f: DensePoly[ModInt], mut budget: Budget) raises PolyError -> SquarefreeDecomposition[ModInt]:
    f.check_invariants(budget)
    var unit = ModInt.zero(f.ring)
    if not f.is_zero():
        unit = _copy_coefficient(f.coefficients[f.degree()], budget)
    var factors = List[Tuple[DensePoly[ModInt], Int]]()
    if f.degree() <= 0:
        return SquarefreeDecomposition[ModInt](unit^, factors^)
    var remaining = _monic(f, budget)
    var scale = 1
    while remaining.degree() > 0:
        budget.charge(1)
        var common = _gcd(remaining, _derivative(remaining, budget), budget)
        var visible = _exact(remaining, common, budget)
        var multiplicity = 1
        while visible.degree() > 0:
            budget.charge(1)
            var shared = _gcd(visible, common, budget)
            var block = _exact(visible, shared, budget)
            if block.degree() > 0:
                _allocation[Tuple[DensePoly[ModInt], Int]](1, budget)
                factors.append((block^, _product(scale, multiplicity)))
            visible = shared^
            common = _exact(common, visible, budget)
            multiplicity = _sum(multiplicity, 1)
        if common.degree() <= 0:
            break
        # Every remaining exponent is divisible by p. In the prime field,
        # inverse Frobenius fixes coefficients, so only exponents change.
        var p = Int(f.ring.value)
        var count = common.degree() // p + 1
        _allocation[ModInt](count, budget)
        var values = List[ModInt](length=count, fill=ModInt.zero(f.ring))
        for i in range(len(common.coefficients)):
            budget.charge(1)
            if i % p == 0:
                values[i // p] = _copy_coefficient(common.coefficients[i], budget)
            elif not common.coefficients[i].is_zero():
                raise PolyError(PolyErrorKind.DOMAIN, "inseparable factor has a non-Frobenius exponent")
        remaining = DensePoly[ModInt].from_coefficients(f.variable, f.ring, values^, budget)
        scale = _product(scale, p)
    for i in range(1, len(factors)):
        var j = i
        while j > 0 and factors[j][1] < factors[j - 1][1]:
            budget.charge(1)
            var held = (DensePoly[ModInt](f.variable, f.ring.copy()), 0)
            swap(held, factors[j])
            swap(held, factors[j - 1])
            swap(held, factors[j])
            j -= 1
    return SquarefreeDecomposition[ModInt](unit^, factors^)
