"""Complete certified real roots for polynomials with algebraic coefficients."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .coefficients import QQ
from .dense import DensePoly
from .certificates import ResultantResult, RootCountCertificate, check_resultant_identity, check_root_count, _dense_equal, _equal
from .algebraic_types import RealAlgebraic, _check_real, _copy_real, _same_root
from .algebraic_arithmetic import _representative
from .algebraic_evaluation import AlgebraicEvaluationCertificate, _evaluate_cached, check_algebraic_evaluation
from ._algebraic_ops import _embed, _norm_expression, _cauchy_bound, _is_zero_real
from ._algebraic import _root_cached, _roots_cached, _NormalizationCache, _refine
from ._dense_ops import _exact
from ._root_ops import _compare as _compare_q
from .roots import count_roots
from .resultant import resultant


@fieldwise_init
struct AlgebraicRootsCertificate(Copyable, Movable, Writable):
    """A nonzero rational annihilator, all its real roots, and exact filtering proofs."""
    var coefficients: List[RealAlgebraic]
    """Normalized representatives of every supplied coefficient, in ascending degree."""
    var coefficient_factors: List[DensePoly[QQ]]
    """Multipliers binding the normalized coefficient definitions to the originals."""
    var eliminations: List[ResultantResult[QQ]]
    """Successive identities eliminating the nonrational coefficient variables."""
    var annihilator: DensePoly[QQ]
    """A nonzero rational polynomial containing every possible real root."""
    var candidates: List[RealAlgebraic]
    """Every distinct real root of the annihilator, in strictly increasing order."""
    var candidate_factors: List[DensePoly[QQ]]
    """Multipliers proving each candidate definition divides the annihilator."""
    var all_roots: RootCountCertificate
    """The total root count on an interval containing the Cauchy root bound."""
    var evaluations: List[AlgebraicEvaluationCertificate]
    """Exact evaluations of the original polynomial at every candidate."""
    var accepted: List[Int]
    """Increasing candidate indices whose evaluations are zero."""
    var roots: List[RealAlgebraic]
    """The selected real roots in strictly increasing order."""


def algebraic_roots(coefficients: List[RealAlgebraic], mut budget: Budget) raises PolyError -> AlgebraicRootsCertificate:
    """Find every distinct real root of an algebraic-coefficient polynomial.

    Args:
        coefficients: Exact real coefficients in ascending degree, including optional trailing zeros.
        budget: The shared work allowance.

    Returns:
        Sorted distinct real roots with elimination, completeness and filtering witnesses.

    Raises:
        PolyError: DOMAIN for zero or invalid value witnesses, LIMIT or CAPACITY on work limits.
    """
    var cache = _NormalizationCache()
    return _roots_with_cache(coefficients, cache, budget)


def _roots_with_cache(coefficients: List[RealAlgebraic], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> AlgebraicRootsCertificate:
    _allocation[RealAlgebraic](len(coefficients), budget)
    _allocation[DensePoly[QQ]](len(coefficients), budget)
    var normalized = List[RealAlgebraic]()
    var coefficient_factors = List[DensePoly[QQ]]()
    for coefficient in coefficients:
        budget.charge(1)
        if not _check_real(coefficient, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic coefficient")
        var value = _root_cached(coefficient.polynomial, coefficient.low, coefficient.high, cache, budget)
        coefficient_factors.append(_exact(coefficient.polynomial, value.polynomial, budget))
        normalized.append(value^)
    var expression = _norm_expression(normalized, budget)
    if expression.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "the zero polynomial has no finite root list")
    var eliminations = List[ResultantResult[QQ]]()
    for i in range(len(normalized)):
        budget.charge(1)
        if normalized[i].rational:
            continue
        var variable = Int64(_sum(i, 1))
        var proof = resultant(expression, _embed(normalized[i].polynomial, variable, budget), variable, budget)
        expression = proof.resultant._copy(budget)
        _allocation[ResultantResult[QQ]](1, budget)
        eliminations.append(proof^)
    var h = DensePoly[QQ].from_sparse(expression, budget, variable=Int64(0))
    if h.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic coefficient norm vanished")
    var bound = _cauchy_bound(h, budget)
    var lower = bound.neg(budget)
    var all_roots = count_roots(h, lower, bound, budget)
    var candidates = _roots_cached(h, cache, budget)
    _allocation[DensePoly[QQ]](len(candidates), budget)
    _allocation[AlgebraicEvaluationCertificate](len(candidates), budget)
    var candidate_factors = List[DensePoly[QQ]]()
    var evaluations = List[AlgebraicEvaluationCertificate]()
    var accepted = List[Int]()
    var roots = List[RealAlgebraic]()
    for i in range(len(candidates)):
        budget.charge(1)
        while _compare_q(candidates[i].low, lower, budget) < 0 or _compare_q(candidates[i].high, bound, budget) > 0:
            budget.charge(1)
            _refine(candidates[i], budget)
        candidate_factors.append(_exact(h, candidates[i].polynomial, budget))
        var evaluation = _evaluate_cached(coefficients, candidates[i], cache, budget)
        if _is_zero_real(evaluation.value, budget):
            _allocation[Int](1, budget)
            _allocation[RealAlgebraic](1, budget)
            accepted.append(i)
            roots.append(_copy_real(candidates[i], budget))
        evaluations.append(evaluation^)
    return AlgebraicRootsCertificate(normalized^, coefficient_factors^, eliminations^, h^, candidates^, candidate_factors^, all_roots^, evaluations^, accepted^, roots^)


def check_algebraic_roots(coefficients: List[RealAlgebraic], certificate: AlgebraicRootsCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Verify root completeness and filtering without elimination or root isolation.

    Args:
        coefficients: The original ascending coefficient list.
        certificate: All elimination, count, candidate and evaluation witnesses.
        budget: The shared work allowance.

    Returns:
        Whether the supplied roots are exactly all distinct real roots.

    Raises:
        PolyError: DOMAIN for malformed polynomials, LIMIT or CAPACITY on work limits.
    """
    ref proof = certificate
    var n = len(coefficients)
    if len(proof.coefficients) != n or len(proof.coefficient_factors) != n:
        return False
    for i in range(n):
        budget.charge(1)
        if not _check_real(coefficients[i], budget) or not _representative(proof.coefficients[i], coefficients[i], proof.coefficient_factors[i], budget):
            return False
    var expression = _norm_expression(proof.coefficients, budget)
    if expression.is_zero():
        return False
    var elimination_index = 0
    for i in range(n):
        budget.charge(1)
        if proof.coefficients[i].rational:
            continue
        if elimination_index == len(proof.eliminations):
            return False
        ref step = proof.eliminations[elimination_index]
        var definition = _embed(proof.coefficients[i].polynomial, Int64(_sum(i, 1)), budget)
        if not check_resultant_identity(expression, definition, step, budget):
            return False
        expression = step.resultant._copy(budget)
        elimination_index += 1
    ref h = proof.annihilator
    if elimination_index != len(proof.eliminations) or h.is_zero() or h.variable != 0 or not _equal(expression, _embed(h, 0, budget), budget):
        return False
    var bound = _cauchy_bound(h, budget)
    if _compare_q(proof.all_roots.a, bound.neg(budget), budget) > 0 or _compare_q(proof.all_roots.b, bound, budget) < 0 or not check_root_count(h, proof.all_roots, budget):
        return False
    var count = len(proof.candidates)
    if proof.all_roots.count != count or len(proof.candidate_factors) != count or len(proof.evaluations) != count or len(proof.accepted) != len(proof.roots):
        return False
    var accepted_index = 0
    for i in range(count):
        budget.charge(1)
        ref candidate = proof.candidates[i]
        if not _check_real(candidate, budget) or not _dense_equal(candidate.polynomial.mul(proof.candidate_factors[i], budget), h, budget):
            return False
        if _compare_q(candidate.low, proof.all_roots.a, budget) < 0 or _compare_q(candidate.high, proof.all_roots.b, budget) > 0:
            return False
        if i:
            var separation = _compare_q(proof.candidates[i - 1].high, candidate.low, budget)
            if separation > 0 or (separation == 0 and proof.candidates[i - 1].rational and candidate.rational):
                return False
        if not check_algebraic_evaluation(coefficients, candidate, proof.evaluations[i], budget):
            return False
        if _is_zero_real(proof.evaluations[i].value, budget):
            if accepted_index == len(proof.accepted) or proof.accepted[accepted_index] != i or not _same_root(proof.roots[accepted_index], candidate, budget):
                return False
            accepted_index += 1
    return accepted_index == len(proof.accepted)
