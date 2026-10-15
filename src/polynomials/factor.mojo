"""Univariate integer and rational factorization by Hensel–Zassenhaus (Z-6)."""
from apn_mojo import Integer
from .coefficients import CoefficientRing, ZZ, QQ, ModInt, Modulus, NoContext, _pow
from .dense import DensePoly
from ._subsets import _Recombination
from ._coefficient_ops import _copy as _copy_coefficient
from .errors import Budget, PolyError, _allocation, _sum
from .certificates import Factorization
from .squarefree import squarefree
from .finite_field import factor_mod_p
from .primes import prime, _residue
from ._dense_ops import _copy, _constant, _monic, _derivative, _gcd
from ._images import _integer_dense as _integer, _rational_dense as _rational, _mod_image
from ._coefficient_ops import _quotient
from ._hensel import _HenselState, _hensel, _truncate, _trial_divide
from ._integer_ops import _content, _primitive, _divide_scalar
from ._lattice_factor import _van_hoeij_lifted


def _bound(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> ZZ:
    # ||f||_2 <= ||f||_1. This exact, conservative Landau–Mignotte bound
    # avoids a square root; recombination always reconstructs the side
    # whose degree is at most half the original degree.
    var norm = ZZ.zero(NoContext())
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        norm = norm.add(c.neg(budget) if c.value.sign() < 0 else _copy_coefficient(c, budget), budget)
    var power = _pow(ZZ(Integer(2)), UInt64(f.degree() // 2 + 1), budget)
    return norm.mul(power, budget).mul(f.coefficients[f.degree()], budget)


def _recombine(f: DensePoly[ZZ], lifted: List[DensePoly[ZZ]], modulus: ZZ, mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    var search = _Recombination(len(lifted), budget)
    var remaining = _copy(f, budget)
    var result = List[DensePoly[ZZ]]()
    while search.ready():
        budget.charge(1)
        var degree = 0
        for i in search.selection.indices:
            budget.charge(1)
            degree = _sum(degree, lifted[search.active[i]].degree())
        var complement = degree > remaining.degree() // 2
        _allocation[Bool](len(search.active), budget)
        var chosen = List[Bool](length=len(search.active), fill=complement)
        for i in search.selection.indices:
            budget.charge(1)
            chosen[i] = not complement
        var candidate = _constant(remaining, remaining.coefficients[remaining.degree()], budget)
        for i in range(len(search.active)):
            budget.charge(1)
            if chosen[i]:
                candidate = _truncate(candidate.mul(lifted[search.active[i]], budget), modulus, budget)
        candidate = _primitive(candidate, budget)
        var quotient = Optional[DensePoly[ZZ]](None)
        if candidate.degree() > 0 and candidate.degree() < remaining.degree():
            if candidate.coefficients[0].is_zero() or Bool(_quotient(remaining.coefficients[0], candidate.coefficients[0], budget)):
                quotient = _trial_divide(remaining, candidate, budget)
        if quotient:
            _allocation[DensePoly[ZZ]](1, budget)
            if complement:
                # The selected (smaller-cardinality) side is irreducible:
                # every proper subset was tested earlier. Reconstructing
                # its lower-degree complement only supplies the quotient.
                result.append(quotient.take())
                remaining = candidate^
            else:
                result.append(candidate^)
                remaining = quotient.take()
            search.accept(budget)
        else:
            search.advance(budget)
    if remaining.degree() > 0:
        _allocation[DensePoly[ZZ]](1, budget)
        result.append(remaining^)
    return result^


def _squarefree_factors(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    if f.degree() <= 1:
        _allocation[DensePoly[ZZ]](1, budget)
        return [_copy(f, budget)]
    var index = 0
    while True:
        budget.charge(1)
        var p = prime(index, budget)
        index = _sum(index, 1)
        if _residue(f.coefficients[f.degree()].value, p, budget) == 0:
            continue
        var ctx = Modulus(_prime_value=p)
        var image = _monic(_mod_image(f, ctx, budget), budget)
        if _gcd(image, _derivative(image, budget), budget).degree() > 0:
            continue
        var modular = factor_mod_p(image, budget, seed=0)
        if len(modular) == 1:
            _allocation[DensePoly[ZZ]](1, budget)
            return [_copy(f, budget)]
        # Measured crossover: exact LLL loses on 8/16-way splitting but
        # wins at 20. Peel singleton groups before paying for the lattice.
        if len(modular) >= 20:
            return _large_recombine(f, modular, budget)
        var lifting = _hensel(f, modular, _bound(f, budget), budget)
        return _recombine(f, lifting[0], lifting[1], budget)


def _large_recombine(f: DensePoly[ZZ], modular: List[DensePoly[ModInt]], mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    var lifting = _hensel(f, modular, _bound(f, budget), budget)
    var remaining = _copy(f, budget)
    var result = List[DensePoly[ZZ]]()
    var active_modular = List[DensePoly[ModInt]]()
    var active_lifted = List[DensePoly[ZZ]]()
    for i in range(len(modular)):
        budget.charge(1)
        var candidate = _primitive(_truncate(lifting[0][i].scale(remaining.coefficients[remaining.degree()], budget), lifting[1], budget), budget)
        var quotient = _trial_divide(remaining, candidate, budget)
        if candidate.degree() > 0 and quotient:
            remaining = quotient.take()
            _allocation[DensePoly[ZZ]](1, budget)
            result.append(candidate^)
        else:
            _allocation[DensePoly[ModInt]](1, budget)
            _allocation[DensePoly[ZZ]](1, budget)
            active_modular.append(_copy(modular[i], budget))
            active_lifted.append(_copy(lifting[0][i], budget))
    if remaining.degree() > 0:
        var rest: List[DensePoly[ZZ]]
        if len(active_modular) >= 20:
            var state = _HenselState(active_lifted^, _copy_coefficient(lifting[1], budget), List[DensePoly[ModInt]]())
            rest = _van_hoeij_lifted(remaining, active_modular, _bound(remaining, budget), state, budget)
        else:
            rest = _recombine(remaining, active_lifted, lifting[1], budget)
        for g in rest:
            budget.charge(1)
            _allocation[DensePoly[ZZ]](1, budget)
            result.append(_copy(g, budget))
    return result^


def _less[C: CoefficientRing](a: DensePoly[C], b: DensePoly[C], mut budget: Budget) raises PolyError -> Bool:
    budget.charge(1)
    if a.degree() != b.degree():
        return a.degree() < b.degree()
    for i in range(a.degree(), -1, -1):
        budget.charge(1)
        var delta = a.coefficients[i].sub(b.coefficients[i], budget)
        if not delta.is_zero():
            comptime if C == ZZ:
                return rebind[ZZ](delta).value.sign() < 0
            elif C == QQ:
                return rebind[QQ](delta).value.sign() < 0
    return False


def _sort[C: CoefficientRing](mut factors: List[Tuple[DensePoly[C], Int]], template: DensePoly[C], mut budget: Budget) raises PolyError:
    for i in range(1, len(factors)):
        budget.charge(1)
        var current = (DensePoly[C](template.variable, template.ring.copy()), 0)
        swap(current, factors[i])
        var j = i
        while j > 0 and _less(current[0], factors[j - 1][0], budget):
            var moved = (DensePoly[C](template.variable, template.ring.copy()), 0)
            swap(moved, factors[j - 1])
            swap(moved, factors[j])
            j -= 1
        swap(current, factors[j])


def factor_integer(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Factorization[ZZ]:
    """Factor an integer polynomial into primitive factors with positive leading coefficients.

    Args:
        f: The integer polynomial, including zero or constants.
        budget: The shared work allowance.

    Returns:
        Signed content and factors sorted by degree, then leading-to-constant coefficients.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    var unit = _content(f, budget)
    var factors = List[Tuple[DensePoly[ZZ], Int]]()
    if f.degree() <= 0:
        return Factorization[ZZ](unit^, factors^)
    var primitive = _divide_scalar(f, unit, budget)
    var blocks = squarefree(_rational(primitive, budget), budget)
    for entry in blocks.factors:
        budget.charge(1)
        var block = _primitive(_integer(entry[0], budget), budget)
        var irreducibles = _squarefree_factors(block, budget)
        for i in range(len(irreducibles)):
            budget.charge(1)
            _allocation[Tuple[DensePoly[ZZ], Int]](1, budget)
            var irreducible = DensePoly[ZZ](f.variable, NoContext())
            swap(irreducible, irreducibles[i])
            factors.append((irreducible^, entry[1]))
    _sort(factors, f, budget)
    return Factorization[ZZ](unit^, factors^)


def factor(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> Factorization[QQ]:
    """Factor a rational polynomial into monic irreducibles with multiplicities.

    Args:
        f: The rational polynomial, including zero or constants.
        budget: The shared work allowance.

    Returns:
        The leading coefficient and factors in the deterministic coefficient order.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    var unit = QQ.zero(f.ring)
    if not f.is_zero():
        unit = _copy_coefficient(f.coefficients[f.degree()], budget)
    var factors = List[Tuple[DensePoly[QQ], Int]]()
    if f.degree() > 0:
        var integer_factors = factor_integer(_integer(f, budget), budget)
        for entry in integer_factors.factors:
            budget.charge(1)
            _allocation[Tuple[DensePoly[QQ], Int]](1, budget)
            factors.append((_monic(_rational(entry[0], budget), budget), entry[1]))
        _sort(factors, f, budget)
    return Factorization[QQ](unit^, factors^)
