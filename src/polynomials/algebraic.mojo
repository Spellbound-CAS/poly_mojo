"""Certified exact real algebraic values and comparisons."""
from .errors import Budget, PolyError, PolyErrorKind
from .coefficients import QQ
from .dense import DensePoly
from .algebraic_types import RealAlgebraic, AlgebraicComparisonCertificate, _check_real, _check_comparison
from ._algebraic import _rational, _root, _comparison
from .algebraic_arithmetic import AlgebraicArithmeticCertificate
from ._algebraic_certify import _arithmetic
from .coefficients import NoContext


def real_rational(value: QQ, mut budget: Budget) raises PolyError -> RealAlgebraic:
    """Embed an exact rational as a real algebraic value.

    Args:
        value: The rational number.
        budget: The shared work allowance.

    Returns:
        A linear defining polynomial and exact rational endpoints.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    return _rational(value, budget)


def real_root(f: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> RealAlgebraic:
    """Select the unique real root in an exact rational interval.

    Args:
        f: A nonconstant rational polynomial, with repetitions allowed.
        low: The lower endpoint, or the exact root when endpoints agree.
        high: The upper endpoint; open intervals require nonroot endpoints.
        budget: The shared work allowance.

    Returns:
        A monic irreducible defining polynomial and a checked root witness.

    Raises:
        PolyError: DOMAIN for invalid root selection, LIMIT or CAPACITY on work limits.
    """
    return _root(f, low, high, budget)


def algebraic_compare(a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicComparisonCertificate:
    """Compare exact real algebraic values with an independent witness.

    Args:
        a: The left root witness, with a possibly nonminimal definition.
        b: The right root witness.
        budget: The shared work allowance.

    Returns:
        A certificate whose relation is negative one, zero or positive one.

    Raises:
        PolyError: DOMAIN for invalid value witnesses, LIMIT or CAPACITY on work limits.
    """
    if not _check_real(a, budget) or not _check_real(b, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic value witness")
    return _comparison(a, b, budget)


def algebraic_sign(value: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicComparisonCertificate:
    """Determine the sign of an exact real algebraic value with a comparison proof.

    Args:
        value: The exact real value.
        budget: The shared work allowance.

    Returns:
        A comparison against rational zero, with relation minus one, zero or one.

    Raises:
        PolyError: DOMAIN for an invalid witness, LIMIT or CAPACITY on work limits.
    """
    return algebraic_compare(value, _rational(QQ.zero(NoContext()), budget), budget)


def check_real_algebraic(value: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    """Verify an algebraic value by exact evaluation or a Descartes root bound.

    Args:
        value: The proposed rational value or unique real root witness.
        budget: The shared work allowance.

    Returns:
        Whether the witness defines exactly one real number.

    Raises:
        PolyError: DOMAIN for malformed polynomials, LIMIT or CAPACITY on work limits.
    """
    return _check_real(value, budget)


def check_algebraic_comparison(a: RealAlgebraic, b: RealAlgebraic, certificate: AlgebraicComparisonCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check interval separation or common-root identities without gcd or isolation.

    Args:
        a: The original left value.
        b: The original right value.
        certificate: The comparison witness to check.
        budget: The shared work allowance.

    Returns:
        Whether the stated relation follows from the witnesses.

    Raises:
        PolyError: DOMAIN for malformed polynomials, LIMIT or CAPACITY on work limits.
    """
    return _check_comparison(a, b, certificate, budget)


def algebraic_add(a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicArithmeticCertificate:
    """Add exact algebraic values with a complete arithmetic certificate.

    Args:
        a: The left operand.
        b: The right operand.
        budget: The shared work allowance.

    Returns:
        The arithmetic result and its independent polynomial witnesses.

    Raises:
        PolyError: DOMAIN for invalid operands or inverse of zero, LIMIT or CAPACITY on work limits.
    """
    return _arithmetic(a, b, 0, budget)


def algebraic_subtract(a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicArithmeticCertificate:
    """Subtract exact algebraic values with a complete arithmetic certificate.

    Args:
        a: The left operand.
        b: The right operand.
        budget: The shared work allowance.

    Returns:
        The arithmetic result and its independent polynomial witnesses.

    Raises:
        PolyError: DOMAIN for invalid operands or inverse of zero, LIMIT or CAPACITY on work limits.
    """
    return _arithmetic(a, b, 1, budget)


def algebraic_multiply(a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicArithmeticCertificate:
    """Multiply exact algebraic values with a complete arithmetic certificate.

    Args:
        a: The left operand.
        b: The right operand.
        budget: The shared work allowance.

    Returns:
        The arithmetic result and its independent polynomial witnesses.

    Raises:
        PolyError: DOMAIN for invalid operands or inverse of zero, LIMIT or CAPACITY on work limits.
    """
    return _arithmetic(a, b, 2, budget)


def algebraic_negate(a: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicArithmeticCertificate:
    """Negate exact algebraic values with a complete arithmetic certificate.

    Args:
        a: The left operand.
        budget: The shared work allowance.

    Returns:
        The arithmetic result and its independent polynomial witnesses.

    Raises:
        PolyError: DOMAIN for invalid operands or inverse of zero, LIMIT or CAPACITY on work limits.
    """
    return _arithmetic(a, _rational(QQ.zero(NoContext()), budget), 3, budget)


def algebraic_inverse(a: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicArithmeticCertificate:
    """Invert exact algebraic values with a complete arithmetic certificate.

    Args:
        a: The left operand.
        budget: The shared work allowance.

    Returns:
        The arithmetic result and its independent polynomial witnesses.

    Raises:
        PolyError: DOMAIN for invalid operands or inverse of zero, LIMIT or CAPACITY on work limits.
    """
    return _arithmetic(a, _rational(QQ.zero(NoContext()), budget), 4, budget)
