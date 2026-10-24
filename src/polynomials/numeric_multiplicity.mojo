"""Exact real-root multiplicities with independently checked squarefree evidence.

A Descartes result isolates the roots of g, the squarefree part of f, one
simple root per bracket. When f is u times the product of g_k^(e_k) and the
product of the blocks g_k is a constant times g, the blocks are squarefree
and pairwise coprime, so each root of g is a root of exactly one block: the
one that changes sign on its bracket, or vanishes at an exact root. Its
multiplicity in f is that block's exponent.
"""
from apn_mojo import Rational
from ._root_ops import _sign_at
from .coefficients import CoefficientRing, QQ, ZZ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .certificates import SquarefreeDecomposition, check_squarefree_product
from .squarefree import squarefree
from .roots import RootInterval
from .numeric_options import NumericOptions
from .descartes import DescartesRootResult, descartes_intervals, _exact_input, _check_descartes_roots
from ._images import _integer_dense
from ._integer_ops import _primitive


@fieldwise_init
struct RootMultiplicityCertificate(Copyable, Movable):
    """Prove the multiplicity of every isolated real root through squarefree blocks."""
    var decomposition: SquarefreeDecomposition[QQ]
    """A product decomposition of the exact polynomial into multiplicity blocks."""
    var blocks: List[Int]
    """For each root in increasing order, the index of the block that vanishes there."""
    var multiplicities: List[Int]
    """For each root in increasing order, its multiplicity in the polynomial."""


def _block_vanishes(block: DensePoly[QQ], root: RootInterval, mut budget: Budget) raises PolyError -> Bool:
    # A sign change across an open bracket, or a zero at an exact root.
    if root.low == root.high:
        return _sign_at(block, root.low, budget) == 0
    var left = _sign_at(block, root.low, budget)
    return left != 0 and _sign_at(block, root.high, budget) == -left


def real_root_multiplicities[C: CoefficientRing](poly: DensePoly[C], roots: DescartesRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> RootMultiplicityCertificate:
    """Attach exact multiplicities to a checked complete real-root result.

    Parameters:
        C: ZZ or QQ.

    Args:
        poly: The exact polynomial.
        roots: A complete Descartes result for poly.
        options: The options describing that result.
        budget: The shared work allowance.

    Returns:
        The squarefree decomposition, and for each root its block and multiplicity.

    Raises:
        PolyError: DOMAIN for invalid or incomplete results, LIMIT or CAPACITY on work limits.
    """
    var exact = _exact_input(poly, budget)
    if not roots.complete or not _check_descartes_roots(exact, roots, options, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "multiplicities require a checked complete real isolation")
    var decomposition = squarefree(exact, budget)
    var intervals = descartes_intervals(roots.certificate, budget)
    _allocation[Int](2 * len(intervals), budget)
    var blocks = List[Int](capacity=len(intervals))
    var multiplicities = List[Int](capacity=len(intervals))
    for root in intervals:
        var found = -1
        for k in range(len(decomposition.factors)):
            if _block_vanishes(decomposition.factors[k][0], root, budget):
                found = k
                break
        if found < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "no multiplicity block vanishes at an isolated root")
        blocks.append(found)
        multiplicities.append(decomposition.factors[found][1])
    return RootMultiplicityCertificate(decomposition^, blocks^, multiplicities^)


def check_real_root_multiplicities[C: CoefficientRing](poly: DensePoly[C], roots: DescartesRootResult, options: NumericOptions, certificate: RootMultiplicityCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check root multiplicities without performing squarefree decomposition or gcd.

    Parameters:
        C: ZZ or QQ.

    Args:
        poly: The exact polynomial.
        roots: The claimed complete Descartes result.
        options: The options describing that result.
        certificate: The decomposition and per-root blocks.
        budget: The shared work allowance.

    Returns:
        Whether every per-root multiplicity is proved exactly.

    Raises:
        PolyError: DOMAIN for malformed polynomial data, LIMIT or CAPACITY on work limits.
    """
    var exact = _exact_input(poly, budget)
    if not roots.complete or not _check_descartes_roots(exact, roots, options, budget):
        return False
    if not check_squarefree_product(exact, certificate.decomposition, budget):
        return False
    # The blocks multiply to a constant times the squarefree part g.
    var product = DensePoly[QQ].from_coefficients(exact.variable, NoContext(), [QQ(Rational(1))], budget)
    for entry in certificate.decomposition.factors:
        if entry[0].degree() < 1 or entry[1] < 1:
            return False
        product = product.mul(entry[0], budget)
    var normal = _primitive(_integer_dense(product, budget), budget)
    if normal.coefficients[normal.degree()].value.sign() < 0:
        normal = normal.neg(budget)
    if normal != roots.squarefree:
        return False
    var intervals = descartes_intervals(roots.certificate, budget)
    if len(certificate.blocks) != len(intervals) or len(certificate.multiplicities) != len(intervals):
        return False
    for i in range(len(intervals)):
        budget.charge(1)
        var k = certificate.blocks[i]
        if k < 0 or k >= len(certificate.decomposition.factors):
            return False
        if certificate.multiplicities[i] != certificate.decomposition.factors[k][1]:
            return False
        if not _block_vanishes(certificate.decomposition.factors[k][0], intervals[i], budget):
            return False
    return True
