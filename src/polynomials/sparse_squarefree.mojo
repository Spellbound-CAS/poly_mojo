"""Multiplicity blocks from all partial derivatives, with inverse Frobenius."""
from .coefficients import CoefficientField, CoefficientRing, ModInt, _pow
from .sparse import SparsePoly, Term
from .monomial import Monomial, MonomialOrder
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from ._coefficient_ops import _copy as _copy_coefficient
from ._gcd import _gcd_value
from ._tower import _exact_poly
from ._views import _is_unit
from .evaluate import derivative
from .certificates import _equal


@fieldwise_init
struct SparseSquarefreeDecomposition[C: CoefficientRing](Copyable, Movable, Writable):
    """A sparse polynomial's scalar and monic multiplicity blocks.

    Parameters:
        C: The coefficient ring.
    """
    var unit: Self.C
    """The normalizing scalar, or zero for zero input."""
    var factors: List[Tuple[SparsePoly[Self.C], Int]]
    """Nonconstant blocks in increasing multiplicity, preserving the input order."""


def squarefree_multivariate[C: CoefficientField](f: SparsePoly[C], mut budget: Budget) raises PolyError -> SparseSquarefreeDecomposition[C]:
    """Group sparse factors by multiplicity without factoring into irreducibles.

    Parameters:
        C: A characteristic-zero field or ModInt prime field.

    Args:
        f: The polynomial in any stored order, including zero and constants.
        budget: The shared work allowance.

    Returns:
        The normalizing scalar and monic, pairwise coprime squarefree blocks.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    if f.order.kind != MonomialOrder.KERNEL:
        return _squarefree(f, budget)
    # KERNEL is a storage order, not a multiplicative monomial order.
    # Compute in lex, then normalize each returned block in storage order
    # and account for its scalar to the block's full multiplicity.
    var variables = f.variables(budget)
    var lex = MonomialOrder(MonomialOrder.LEX, variables^, budget)
    var result = _squarefree(f.reorder(lex, budget), budget)
    for i in range(len(result.factors)):
        budget.charge(1)
        var block = result.factors[i][0].reorder(f.order, budget)
        var leading = _copy_coefficient(block.terms[len(block.terms) - 1].coefficient, budget)
        result.unit = result.unit.mul(_pow(leading, UInt64(result.factors[i][1]), budget), budget)
        result.factors[i][0] = block.scale(leading.inverse(budget), budget)
    return result^


def _squarefree[C: CoefficientField](f: SparsePoly[C], mut budget: Budget) raises PolyError -> SparseSquarefreeDecomposition[C]:
    f.check_invariants(budget)
    var unit = C.zero(f.ring)
    if not f.is_zero():
        unit = _copy_coefficient(f.terms[len(f.terms) - 1].coefficient, budget)
    var factors = List[Tuple[SparsePoly[C], Int]]()
    if f.is_zero() or _is_unit(f, budget):
        return SparseSquarefreeDecomposition[C](unit^, factors^)
    var variables = f.variables(budget)
    var remaining = f.scale(unit.inverse(budget), budget)
    var scale = 1
    while not _is_unit(remaining, budget):
        budget.charge(1)
        var common = remaining._copy(budget)
        for v in variables:
            budget.charge(1)
            common = _gcd_value(common, derivative(remaining, v, budget), budget)
            if _is_unit(common, budget):
                break
        var visible = _exact_poly(remaining, common, budget)
        var multiplicity = 1
        while not _is_unit(visible, budget):
            budget.charge(1)
            var shared = _gcd_value(visible, common, budget)
            var block = _exact_poly(visible, shared, budget)
            if not _is_unit(block, budget):
                _allocation[Tuple[SparsePoly[C], Int]](1, budget)
                factors.append((block^, _product(scale, multiplicity)))
            visible = shared^
            common = _exact_poly(common, visible, budget)
            multiplicity = _sum(multiplicity, 1)
        if _is_unit(common, budget):
            break
        comptime if C == ModInt:
            var p = rebind[ModInt.Context](f.ring).value
            _allocation[Term[C]](len(common.terms), budget)
            var terms = List[Term[C]](capacity=len(common.terms))
            for term in common.terms:
                _allocation[Tuple[Int64, UInt32]](len(term.monomial.pairs), budget)
                var pairs = List[Tuple[Int64, UInt32]](capacity=len(term.monomial.pairs))
                for pair in term.monomial.pairs:
                    budget.charge(1)
                    if UInt64(pair[1]) % p != 0:
                        raise PolyError(PolyErrorKind.DOMAIN, "inseparable factor has a non-Frobenius exponent")
                    pairs.append((pair[0], UInt32(UInt64(pair[1]) // p)))
                terms.append(Term[C](Monomial(pairs^, budget), _copy_coefficient(term.coefficient, budget)))
            remaining = SparsePoly[C].from_terms(terms^, f.order, f.ring, budget)
            scale = _product(scale, Int(p))
        else:
            raise PolyError(PolyErrorKind.DOMAIN, "nonconstant inseparable factor in characteristic zero")
    for i in range(1, len(factors)):
        var j = i
        while j > 0 and factors[j][1] < factors[j - 1][1]:
            budget.charge(1)
            var held = (SparsePoly[C].zero(f.order, f.ring, budget), 0)
            swap(held, factors[j])
            swap(held, factors[j - 1])
            swap(held, factors[j])
            j -= 1
    return SparseSquarefreeDecomposition[C](unit^, factors^)


def check_sparse_squarefree_product[C: CoefficientRing](f: SparsePoly[C], decomposition: SparseSquarefreeDecomposition[C], mut budget: Budget) raises PolyError -> Bool:
    """Check a sparse multiplicity product without asserting squarefreeness.

    Parameters:
        C: The coefficient ring.

    Args:
        f: The original polynomial.
        decomposition: The proposed scalar, blocks and multiplicities.
        budget: The shared work allowance.

    Returns:
        Whether all multiplicities are positive and the product equals f.

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    var product = SparsePoly[C].constant(decomposition.unit, f.order, budget)
    for entry in decomposition.factors:
        budget.charge(1)
        entry[0].check_invariants(budget)
        if entry[1] <= 0:
            return False
        product = product.mul(entry[0].pow(UInt64(entry[1]), budget), budget)
    return _equal(product, f, budget)
