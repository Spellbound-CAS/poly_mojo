"""Exact Gaussian-rational factorization by squarefree conjugate norms.

Trager's norm reduction specializes to N(a + i*b) = a*a + b*b over QQ(i).
Integer imaginary shifts separate conjugate factors before rational factoring.
"""
from apn_mojo import ExactComplex
from .coefficients import QQ, NoContext
from .exact_complex import QQi
from .dense import DensePoly
from .errors import Budget, PolyError, _allocation, _sum, _charge_size
from .certificates import Factorization
from .squarefree import squarefree
from .factor import factor
from .calculus import taylor_shift
from ._dense_ops import _copy, _gcd, _derivative


def _lift(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQi]:
    _allocation[QQi](len(f.coefficients), budget)
    var values = List[QQi](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        _charge_size(c.size(), budget)
        values.append(QQi(ExactComplex(c.value)))
    return DensePoly[QQi].from_coefficients(f.variable, NoContext(), values^, budget)


def _norm(f: DensePoly[QQi], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    _allocation[QQ](len(f.coefficients), budget)
    _allocation[QQ](len(f.coefficients), budget)
    var real = List[QQ](capacity=len(f.coefficients))
    var imag = List[QQ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        _charge_size(c.size(), budget)
        real.append(QQ(c.value.real()))
        imag.append(QQ(c.value.imag()))
    var a = DensePoly[QQ].from_coefficients(f.variable, NoContext(), real^, budget)
    var b = DensePoly[QQ].from_coefficients(f.variable, NoContext(), imag^, budget)
    return a.mul(a, budget).add(b.mul(b, budget), budget)


def _less(a: DensePoly[QQi], b: DensePoly[QQi], mut budget: Budget) raises PolyError -> Bool:
    budget.charge(1)
    if a.degree() != b.degree():
        return a.degree() < b.degree()
    for i in range(a.degree(), -1, -1):
        var delta = a.coefficients[i].sub(b.coefficients[i], budget)
        if delta.value.real().sign() != 0:
            return delta.value.real().sign() < 0
        if delta.value.imag().sign() != 0:
            return delta.value.imag().sign() < 0
    return False


def factor_gaussian(f: DensePoly[QQi], mut budget: Budget) raises PolyError -> Factorization[QQi]:
    """Factor over QQ(i) into monic irreducibles with exact multiplicities.

    Args:
        f: The Gaussian-rational polynomial, including zero and constants.
        budget: The shared allowance for shifts, norms, rational factors and GCDs.

    Returns:
        The leading coefficient and factors sorted by degree, then coefficients
        from leading to constant, comparing real before imaginary components.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    var blocks = squarefree(f, budget)
    var factors = List[Tuple[DensePoly[QQi], Int]]()
    for block in blocks.factors:
        budget.charge(1)
        if block[0].degree() == 1:
            _allocation[Tuple[DensePoly[QQi], Int]](1, budget)
            factors.append((_copy(block[0], budget), block[1]))
            continue
        var shift = 0
        while True:
            budget.charge(1)
            var point = QQi(ExactComplex(0, shift))
            var shifted = taylor_shift(block[0], point, budget)
            var norm = _norm(shifted, budget)
            if _gcd(norm, _derivative(norm, budget), budget).degree() > 0:
                shift = _sum(shift, 1)
                continue
            var rational = factor(norm, budget)
            for entry in rational.factors:
                budget.charge(1)
                var recovered = _gcd(shifted, _lift(entry[0], budget), budget)
                if recovered.degree() > 0:
                    _allocation[Tuple[DensePoly[QQi], Int]](1, budget)
                    factors.append((taylor_shift(recovered, point.neg(budget), budget), block[1]))
            break
    for i in range(1, len(factors)):
        var j = i
        while j > 0 and _less(factors[j][0], factors[j - 1][0], budget):
            budget.charge(1)
            var held = (DensePoly[QQi](f.variable, f.ring.copy()), 0)
            swap(held, factors[j])
            swap(held, factors[j - 1])
            swap(held, factors[j])
            j -= 1
    var unit = QQi.zero(f.ring)
    swap(unit, blocks.unit)
    return Factorization[QQi](unit^, factors^)
