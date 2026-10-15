"""Independent checks of algebraic arithmetic from polynomial identities."""
from std.collections import Optional
from .errors import Budget, PolyError
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .certificates import ResultantResult, RootCountCertificate, check_root_count, check_resultant_identity, _dense_equal, _equal
from .algebraic_types import RealAlgebraic, _check_real
from ._root_ops import _compare as _compare_q
from ._algebraic_ops import _embed, _binary_definitions, _unary_definition, _operation_bounds, _rational_operation


@fieldwise_init
struct AlgebraicArithmeticCertificate(Copyable, Movable, Writable):
    """An arithmetic result proved by elimination and unique-root selection."""
    var operation: Int
    """Zero for add, one for subtract, two for multiply, three for negate, four for inverse."""
    var left: RealAlgebraic
    """A normalized and refined witness for the left operand."""
    var right: RealAlgebraic
    """A normalized right operand, unused for unary operations."""
    var left_factor: DensePoly[QQ]
    """Multiplier recovering the original left defining polynomial."""
    var right_factor: DensePoly[QQ]
    """Multiplier recovering the original right definition, unused for unary operations."""
    var value: RealAlgebraic
    """The certified arithmetic result."""
    var annihilator: DensePoly[QQ]
    """A nonzero polynomial vanishing at the intended arithmetic result."""
    var value_factor: DensePoly[QQ]
    """Multiplier expressing the annihilator as a multiple of the result definition."""
    var elimination: Optional[ResultantResult[QQ]]
    """The binary elimination identity, absent for direct rational and unary cases."""
    var selection: Optional[RootCountCertificate]
    """One annihilator root enclosing the operand-image bounds and result interval."""


def _representative(value: RealAlgebraic, original: RealAlgebraic, multiplier: DensePoly[QQ], mut budget: Budget) raises PolyError -> Bool:
    if not _check_real(value, budget) or not _dense_equal(value.polynomial.mul(multiplier, budget), original.polynomial, budget):
        return False
    if original.rational:
        return value.rational and _compare_q(value.low, original.low, budget) == 0
    return _compare_q(value.low, original.low, budget) >= 0 and _compare_q(value.high, original.high, budget) <= 0


def check_algebraic_arithmetic(a: RealAlgebraic, b: Optional[RealAlgebraic], certificate: AlgebraicArithmeticCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Verify algebraic arithmetic without computing resultants or isolating roots.

    Args:
        a: The original left operand.
        b: The original right operand, absent for negate and inverse.
        certificate: The operation, result and complete polynomial witnesses.
        budget: The shared work allowance.

    Returns:
        Whether the certificate proves its stated arithmetic result.

    Raises:
        PolyError: DOMAIN for malformed polynomials, LIMIT or CAPACITY on work limits.
    """
    ref proof = certificate
    var op = proof.operation
    if op < 0 or op > 4 or Bool(b) != (op < 3):
        return False
    if not _check_real(a, budget) or not _representative(proof.left, a, proof.left_factor, budget) or not _check_real(proof.value, budget):
        return False
    if op < 3:
        if not _check_real(b.value(), budget) or not _representative(proof.right, b.value(), proof.right_factor, budget):
            return False
    if op == 4 and proof.left.rational and proof.left.low.is_zero():
        return False
    var exact = proof.left.rational and (op >= 3 or proof.right.rational)
    var zero_product = op == 2 and ((proof.left.rational and proof.left.low.is_zero()) or (proof.right.rational and proof.right.low.is_zero()))
    if exact or zero_product:
        var expected = QQ.zero(NoContext())
        if exact:
            expected = _rational_operation(proof.left.low, proof.right.low, op, budget)
        return proof.value.rational and _compare_q(proof.value.low, expected, budget) == 0
    if not proof.selection or proof.annihilator.is_zero() or proof.annihilator.variable != 0:
        return False
    if op < 3:
        if not proof.elimination:
            return False
        var definitions = _binary_definitions(proof.left, proof.right, op, budget)
        ref identity = proof.elimination.value()
        if not _equal(identity.resultant, _embed(proof.annihilator, 0, budget), budget) or not check_resultant_identity(definitions[0], definitions[1], identity, budget):
            return False
    else:
        if op == 4 and proof.left.low.value.sign() * proof.left.high.value.sign() <= 0:
            return False
        if not _dense_equal(proof.annihilator, _unary_definition(proof.left, op, budget), budget):
            return False
    if not _dense_equal(proof.value.polynomial.mul(proof.value_factor, budget), proof.annihilator, budget):
        return False
    ref selection = proof.selection.value()
    if selection.count != 1 or not check_root_count(proof.annihilator, selection, budget):
        return False
    var bounds = _operation_bounds(proof.left, proof.right, op, budget)
    return _compare_q(selection.a, bounds[0], budget) <= 0 and _compare_q(bounds[1], selection.b, budget) <= 0 and _compare_q(selection.a, proof.value.low, budget) <= 0 and _compare_q(proof.value.high, selection.b, budget) <= 0
