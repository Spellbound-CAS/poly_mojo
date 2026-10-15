"""Complete budgeted Kronecker recombination for QQ(i) and prime fields.

Mixed-radix substitution is injective on every possible divisor's degree box.
Univariate factors supply candidates; only exact sparse divisibility accepts
them. Exhausting subsets proves irreducibility, though the search can be large.
"""
from .coefficients import CoefficientField, ModInt, _pow
from .exact_complex import QQi
from .dense import DensePoly
from .sparse import SparsePoly, Term
from .monomial import Monomial
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from .certificates import Factorization, _equal
from .multivariate_factor import MultivariateFactorization, _sort, _unit
from .finite_field import factor_finite_field
from .gaussian_factor import factor_gaussian
from .sparse_squarefree import squarefree_multivariate
from ._subsets import _Recombination
from ._coefficient_ops import _copy as _copy_coefficient
from ._views import _is_unit
from .division import divide_exact


def _decode[C: CoefficientField](f: DensePoly[C], variables: List[Int64], radices: List[Int], template: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var terms = List[Term[C]]()
    for i in range(len(f.coefficients)):
        budget.charge(1)
        if f.coefficients[i].is_zero():
            continue
        _allocation[Tuple[Int64, UInt32]](len(variables), budget)
        var pairs = List[Tuple[Int64, UInt32]](capacity=len(variables))
        var power = i
        for j in range(len(variables)):
            budget.charge(1)
            var exponent = power % radices[j]
            power //= radices[j]
            if exponent:
                pairs.append((variables[j], UInt32(exponent)))
        if power:
            raise PolyError(PolyErrorKind.CAPACITY, "Kronecker exponent exceeds its degree box")
        _allocation[Term[C]](1, budget)
        terms.append(Term[C](Monomial(pairs^, budget), _copy_coefficient(f.coefficients[i], budget)))
    var result = SparsePoly[C].from_terms(terms^, template.order, template.ring, budget)
    return result.scale(result.terms[len(result.terms) - 1].coefficient.inverse(budget), budget)


def _factor_block[C: CoefficientField](f: SparsePoly[C], mut budget: Budget, seed: UInt64) raises PolyError -> MultivariateFactorization[C]:
    f.check_invariants(budget)
    var factors = List[Tuple[SparsePoly[C], Int]]()
    if f.is_zero() or _is_unit(f, budget):
        return MultivariateFactorization[C](_unit(f, budget), factors^)
    var variables = f.variables(budget)
    _allocation[Int](len(variables), budget)
    _allocation[Int](len(variables), budget)
    var radices = List[Int](capacity=len(variables))
    var strides = List[Int](capacity=len(variables))
    var size = 1
    for v in variables:
        budget.charge(1)
        var radix = _sum(f.degree_in(v, budget), 1)
        radices.append(radix)
        strides.append(size)
        size = _product(size, radix)
    # Charge the entire encoded degree box before allocating dense storage.
    _allocation[C](size, budget)
    var coefficients = List[C](length=size, fill=C.zero(f.ring))
    for term in f.terms:
        budget.charge(1)
        var index = 0
        for j in range(len(variables)):
            index = _sum(index, _product(Int(term.monomial.exponent(variables[j], budget)), strides[j]))
        coefficients[index] = _copy_coefficient(term.coefficient, budget)
    var image = DensePoly[C].from_coefficients(variables[0], f.ring, coefficients^, budget)
    var split: Factorization[C]
    comptime if C == ModInt:
        split = rebind_var[Factorization[C]](factor_finite_field(rebind[DensePoly[ModInt]](image), budget, seed))
    else:
        split = rebind_var[Factorization[C]](factor_gaussian(rebind[DensePoly[QQi]](image), budget))
    var count = 0
    for entry in split.factors:
        budget.charge(1)
        count = _sum(count, entry[1])
    _allocation[Int](count, budget)
    var occurrences = List[Int](capacity=count)
    for i in range(len(split.factors)):
        for _ in range(split.factors[i][1]):
            budget.charge(1)
            occurrences.append(i)
    var remaining = f._copy(budget)
    var search = _Recombination(count, budget)
    while search.ready():
        budget.charge(1)
        _allocation[C](1, budget)
        var candidate_image = DensePoly[C].from_coefficients(image.variable, f.ring, [C.one(f.ring)], budget)
        for i in search.selection.indices:
            budget.charge(1)
            candidate_image = candidate_image.mul(split.factors[occurrences[search.active[i]]][0], budget)
        var candidate = _decode(candidate_image, variables, radices, f, budget)
        var quotient = divide_exact(remaining, candidate, budget)
        if quotient:
            _allocation[Tuple[SparsePoly[C], Int]](1, budget)
            factors.append((candidate^, 1))
            remaining = quotient.take()
            search.accept(budget)
        else:
            search.advance(budget)
    if not _is_unit(remaining, budget):
        var leading = _copy_coefficient(remaining.terms[len(remaining.terms) - 1].coefficient, budget)
        _allocation[Tuple[SparsePoly[C], Int]](1, budget)
        factors.append((remaining.scale(leading.inverse(budget), budget), 1))
        remaining = SparsePoly[C].constant(leading, f.order, budget)
    _sort(factors, budget)
    # Equal factors from repeated image factors are adjacent after sorting.
    _allocation[Tuple[SparsePoly[C], Int]](len(factors), budget)
    var distinct = List[Tuple[SparsePoly[C], Int]](capacity=len(factors))
    for i in range(len(factors)):
        budget.charge(1)
        if len(distinct) and _equal(distinct[len(distinct) - 1][0], factors[i][0], budget):
            distinct[len(distinct) - 1][1] = _sum(distinct[len(distinct) - 1][1], 1)
        else:
            var item = (SparsePoly[C].zero(f.order, f.ring, budget), 0)
            swap(item, factors[i])
            distinct.append(item^)
    return MultivariateFactorization[C](_unit(remaining, budget), distinct^)


def _factor[C: CoefficientField](f: SparsePoly[C], mut budget: Budget, seed: UInt64) raises PolyError -> MultivariateFactorization[C]:
    # Strip multiplicities before building images: (x+y)^p can have a tiny
    # sparse block but a very large, highly repeated univariate image.
    var blocks = squarefree_multivariate(f, budget)
    var unit = _copy_coefficient(blocks.unit, budget)
    var factors = List[Tuple[SparsePoly[C], Int]]()
    for entry in blocks.factors:
        budget.charge(1)
        var split = _factor_block(entry[0], budget, seed)
        unit = unit.mul(_pow(split.unit, UInt64(entry[1]), budget), budget)
        _allocation[Tuple[SparsePoly[C], Int]](len(split.factors), budget)
        for i in range(len(split.factors)):
            var exponent = _product(entry[1], split.factors[i][1])
            var p = SparsePoly[C].zero(f.order, f.ring, budget)
            swap(p, split.factors[i][0])
            factors.append((p^, exponent))
    _sort(factors, budget)
    return MultivariateFactorization[C](unit^, factors^)


def factor_multivariate_gaussian(f: SparsePoly[QQi], mut budget: Budget) raises PolyError -> MultivariateFactorization[QQi]:
    """Factor a sparse Gaussian-rational polynomial by exact recombination.

    Args:
        f: The polynomial in any order, including zero and constants.
        budget: The shared allowance, including the dense substitution image.

    Returns:
        A scalar and distinct monic irreducibles with positive multiplicities,
        sorted by total degree, monomials and real-then-imaginary coefficients.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    return _factor(f, budget, 0)


def factor_multivariate_finite_field(f: SparsePoly[ModInt], mut budget: Budget, seed: UInt64 = 0) raises PolyError -> MultivariateFactorization[ModInt]:
    """Factor any sparse prime-field polynomial, including inseparable inputs.

    Args:
        f: The polynomial in any order, including zero and constants.
        budget: The shared allowance, including the dense substitution image.
        seed: The deterministic univariate splitting seed; defaults to zero.

    Returns:
        A scalar and distinct monic irreducibles with positive multiplicities,
        sorted by total degree, monomials and residue coefficients.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    return _factor(f, budget, seed)
