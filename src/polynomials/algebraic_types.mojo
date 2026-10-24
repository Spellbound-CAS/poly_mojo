"""Algebraic value witnesses and independent comparison checks."""
from std.collections import Optional
from .errors import Budget, PolyError, _allocation
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .certificates import _dense_equal
from .roots import RootInterval
from .descartes import _one_root
from ._dense_ops import _copy as _copy_poly
from ._root_ops import _compare as _compare_q, _copy as _copy_q, _sign_at


@fieldwise_init
struct RealAlgebraic(Copyable, Movable, Writable):
    """An exact real root isolated by Descartes' rule of signs.

    An irrational value's polynomial has nonzero values of opposite sign at
    the endpoints and Descartes bound 1 on the open interval between them,
    which also excludes a repeated root there.
    """
    var polynomial: DensePoly[QQ]
    """The monic defining polynomial in internal variable zero."""
    var low: QQ
    """The lower endpoint, equal to high for a rational value."""
    var high: QQ
    """The upper endpoint of the root interval."""
    var rational: Bool
    """Whether the endpoints specify an exact rational value."""


def _rational(q: QQ, mut budget: Budget) raises PolyError -> RealAlgebraic:
    _allocation[QQ](2, budget)
    var p = DensePoly[QQ].from_coefficients(0, NoContext(), [q.neg(budget), QQ.one(NoContext())], budget)
    return RealAlgebraic(p^, _copy_q(q, budget), _copy_q(q, budget), True)


def _copy_real(a: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    return RealAlgebraic(_copy_poly(a.polynomial, budget), _copy_q(a.low, budget), _copy_q(a.high, budget), a.rational)



@fieldwise_init
struct AlgebraicComparisonCertificate(Copyable, Movable, Writable):
    """Interval separation or a shared-root proof of exact comparison."""
    var relation: Int
    """Negative one, zero, or positive one for the certified comparison."""
    var left: RealAlgebraic
    """A refined witness for the left operand."""
    var right: RealAlgebraic
    """A refined witness for the right operand."""
    var common: DensePoly[QQ]
    """A common factor for irrational equality."""
    var left_cofactor: DensePoly[QQ]
    """The quotient of the left definition by the common factor."""
    var right_cofactor: DensePoly[QQ]
    """The quotient of the right definition by the common factor."""
    var shared: Optional[RootInterval]
    """An open interval inside both refined intervals where the common factor has Descartes bound 1 and changes sign."""


def _check_real(a: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    a.polynomial.check_invariants(budget)
    if a.polynomial.variable != 0 or a.polynomial.degree() < 1:
        return False
    if _compare_q(a.polynomial.coefficients[a.polynomial.degree()], QQ.one(NoContext()), budget) != 0:
        return False
    if a.rational:
        return a.polynomial.degree() == 1 and _compare_q(a.low, a.high, budget) == 0 and _sign_at(a.polynomial, a.low, budget) == 0
    return _one_root(a.polynomial, a.low, a.high, budget)


def _same_root(a: RealAlgebraic, original: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    if not _check_real(a, budget) or not _dense_equal(a.polynomial, original.polynomial, budget):
        return False
    if original.rational:
        return a.rational and _compare_q(a.low, original.low, budget) == 0
    return _compare_q(a.low, original.low, budget) >= 0 and _compare_q(a.high, original.high, budget) <= 0


def _check_comparison(a: RealAlgebraic, b: RealAlgebraic, proof: AlgebraicComparisonCertificate, mut budget: Budget) raises PolyError -> Bool:
    if not _check_real(a, budget) or not _check_real(b, budget):
        return False
    if not _same_root(proof.left, a, budget) or not _same_root(proof.right, b, budget):
        return False
    if proof.relation == -1 or proof.relation == 1:
        var comparison = _compare_q(proof.left.high, proof.right.low, budget) if proof.relation == -1 else _compare_q(proof.right.high, proof.left.low, budget)
        return comparison < 0 or (comparison == 0 and not (proof.left.rational and proof.right.rational))
    if proof.relation != 0:
        return False
    if proof.left.rational or proof.right.rational:
        # Production normalizes rational roots to a linear polynomial, but the
        # checker does not assume irreducibility for adversarial witnesses.
        ref point = proof.left.low if proof.left.rational else proof.right.low
        ref other = proof.right if proof.left.rational else proof.left
        if other.rational:
            return _compare_q(point, other.low, budget) == 0
        return _compare_q(other.low, point, budget) < 0 and _compare_q(point, other.high, budget) < 0 and _sign_at(other.polynomial, point, budget) == 0
    if not proof.shared:
        return False
    # The common factor's single root in the shared interval is a root of both
    # definitions inside both isolating intervals, so it is both values.
    ref shared = proof.shared.value()
    if not _one_root(proof.common, shared.low, shared.high, budget):
        return False
    if not _dense_equal(proof.common.mul(proof.left_cofactor, budget), a.polynomial, budget) or not _dense_equal(proof.common.mul(proof.right_cofactor, budget), b.polynomial, budget):
        return False
    return _compare_q(shared.low, proof.left.low, budget) >= 0 and _compare_q(shared.high, proof.left.high, budget) <= 0 and _compare_q(shared.low, proof.right.low, budget) >= 0 and _compare_q(shared.high, proof.right.high, budget) <= 0
