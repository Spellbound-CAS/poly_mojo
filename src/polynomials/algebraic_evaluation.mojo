"""Certified Horner evaluation over exact real algebraic coefficients."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .coefficients import QQ, NoContext
from .algebraic_types import RealAlgebraic, _check_real, _copy_real, _same_root
from .algebraic_arithmetic import AlgebraicArithmeticCertificate, check_algebraic_arithmetic
from ._algebraic_certify import _arithmetic_cached
from ._algebraic import _NormalizationCache
from .algebraic_types import _rational


@fieldwise_init
struct AlgebraicEvaluationCertificate(Copyable, Movable, Writable):
    """A complete sequence of Horner multiplication and addition witnesses."""
    var products: List[AlgebraicArithmeticCertificate]
    """Multiplications by the evaluation point, from highest coefficient down."""
    var sums: List[AlgebraicArithmeticCertificate]
    """Additions of the corresponding coefficients after each multiplication."""
    var value: RealAlgebraic
    """The final exact value, or rational zero for an empty coefficient list."""


def algebraic_evaluate(coefficients: List[RealAlgebraic], x: RealAlgebraic, mut budget: Budget) raises PolyError -> AlgebraicEvaluationCertificate:
    """Evaluate an algebraic-coefficient polynomial with complete arithmetic proofs.

    Args:
        coefficients: The coefficient list in ascending degree, including optional trailing zeros.
        x: The exact real evaluation point.
        budget: The shared work allowance.

    Returns:
        The value and every multiplication and addition witness in Horner order.

    Raises:
        PolyError: DOMAIN for invalid value witnesses, LIMIT or CAPACITY on work limits.
    """
    var cache = _NormalizationCache()
    return _evaluate_cached(coefficients, x, cache, budget)


def _evaluate_cached(coefficients: List[RealAlgebraic], x: RealAlgebraic, mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> AlgebraicEvaluationCertificate:
    if not _check_real(x, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic evaluation point")
    _allocation[AlgebraicArithmeticCertificate](len(coefficients), budget)
    _allocation[AlgebraicArithmeticCertificate](len(coefficients), budget)
    var products = List[AlgebraicArithmeticCertificate]()
    var sums = List[AlgebraicArithmeticCertificate]()
    var value = _rational(QQ.zero(NoContext()), budget)
    for i in range(len(coefficients) - 1, -1, -1):
        budget.charge(1)
        var product = _arithmetic_cached(value, x, 2, cache, budget)
        var sum = _arithmetic_cached(product.value, coefficients[i], 0, cache, budget)
        value = _copy_real(sum.value, budget)
        products.append(product^)
        sums.append(sum^)
    return AlgebraicEvaluationCertificate(products^, sums^, value^)


def check_algebraic_evaluation(coefficients: List[RealAlgebraic], x: RealAlgebraic, certificate: AlgebraicEvaluationCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check all Horner steps without computing algebraic arithmetic results.

    Args:
        coefficients: The original ascending coefficient list.
        x: The original exact evaluation point.
        certificate: The proposed arithmetic chain and final value.
        budget: The shared work allowance.

    Returns:
        Whether the chain proves the stated polynomial evaluation.

    Raises:
        PolyError: DOMAIN for malformed polynomials, LIMIT or CAPACITY on work limits.
    """
    if len(certificate.products) != len(coefficients) or len(certificate.sums) != len(coefficients) or not _check_real(x, budget):
        return False
    var value = _rational(QQ.zero(NoContext()), budget)
    for i in range(len(coefficients)):
        budget.charge(1)
        ref product = certificate.products[i]
        ref sum = certificate.sums[i]
        if product.operation != 2 or sum.operation != 0:
            return False
        if not check_algebraic_arithmetic(value, _copy_real(x, budget), product, budget) or not check_algebraic_arithmetic(product.value, _copy_real(coefficients[len(coefficients) - 1 - i], budget), sum, budget):
            return False
        value = _copy_real(sum.value, budget)
    return _same_root(certificate.value, value, budget)
