"""Exact rational-function calculus, partial fractions and local expansions."""
from apn_mojo import Rational
from .coefficients import QQ, NoContext
from .dense import DensePoly, _slice
from .sparse import SparsePoly
from .fractions import RationalFunction
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _charge_size
from .evaluate import evaluate, derivative
from .calculus import substitute_many, taylor_shift
from .factor import factor
from .certificates import Factorization
from .division import divmod
from ._euclidean import _extended
from ._dense_ops import _exact, _remainder, _constant, _copy
from ._toolkit_exact import _ExactOps
from ._toolkit_ops import _divide_values, _mul_values, _values


def _validate(f: RationalFunction, mut budget: Budget) raises PolyError:
    f.numerator.check_invariants(budget)
    f.denominator.check_invariants(budget)
    f.numerator._same(f.denominator, budget)
    if f.denominator.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "zero rational-function denominator")


def _normalized(p: SparsePoly[QQ], q: SparsePoly[QQ], mut budget: Budget) raises PolyError -> RationalFunction:
    var common = SparsePoly[QQ].zero(p.order, p.ring, budget)
    return RationalFunction._reduce(p, q, common, budget)


def _canonical(f: RationalFunction, mut budget: Budget) raises PolyError -> RationalFunction:
    _validate(f, budget)
    return _normalized(f.numerator, f.denominator, budget)


def _univariate(f: RationalFunction, variable: Int64, mut budget: Budget) raises PolyError -> Tuple[DensePoly[QQ], DensePoly[QQ]]:
    _validate(f, budget)
    return (DensePoly[QQ].from_sparse(f.numerator, budget, variable=variable), DensePoly[QQ].from_sparse(f.denominator, budget, variable=variable))


def rational_evaluate(f: RationalFunction, assignments: List[Tuple[Int64, QQ]], mut budget: Budget) raises PolyError -> QQ:
    """Evaluate a rational function at a nonpole exact assignment.

    Args:
        f: The rational function.
        assignments: Values for all occurring variables.
        budget: The shared work allowance.

    Returns:
        The exact rational value after removing common factors.

    Raises:
        PolyError: DOMAIN for poles or invalid assignments, LIMIT or CAPACITY on work limits.
    """
    var reduced = _canonical(f, budget)
    var denominator = evaluate(reduced.denominator, assignments, budget)
    return evaluate(reduced.numerator, assignments, budget).mul(denominator.inverse(budget), budget)


def rational_derivative(f: RationalFunction, variable: Int64, mut budget: Budget) raises PolyError -> RationalFunction:
    """Differentiate a rational function by the quotient rule.

    Args:
        f: The rational function.
        variable: The differentiation variable.
        budget: The shared work allowance.

    Returns:
        The normalized exact derivative.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    _validate(f, budget)
    var numerator = derivative(f.numerator, variable, budget).mul(f.denominator, budget).sub(f.numerator.mul(derivative(f.denominator, variable, budget), budget), budget)
    var denominator = f.denominator.mul(f.denominator, budget)
    return _normalized(numerator, denominator, budget)


def rational_substitute(f: RationalFunction, assignments: List[Tuple[Int64, SparsePoly[QQ]]], mut budget: Budget) raises PolyError -> RationalFunction:
    """Substitute polynomials simultaneously into a rational function.

    Args:
        f: The rational function.
        assignments: Distinct variables and compatible polynomial replacements.
        budget: The shared work allowance.

    Returns:
        The normalized substituted fraction.

    Raises:
        PolyError: DOMAIN for invalid data or an identically zero substituted denominator, LIMIT or CAPACITY on work limits.
    """
    var reduced = _canonical(f, budget)
    var numerator = substitute_many(reduced.numerator, assignments, budget)
    var denominator = substitute_many(reduced.denominator, assignments, budget)
    return _normalized(numerator, denominator, budget)


def _homogenize(poly: DensePoly[QQ], numerator: DensePoly[QQ], denominator: DensePoly[QQ], degree: Int, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    var result = DensePoly[QQ](poly.variable, NoContext())
    var power = _constant(poly, QQ(Rational(1)), budget)
    for k in range(degree, -1, -1):
        result = result.mul(numerator, budget)
        if k < len(poly.coefficients):
            result = result.add(power.scale(poly.coefficients[k], budget), budget)
        if k:
            power = power.mul(denominator, budget)
    return result^


def rational_compose(outer: RationalFunction, inner: RationalFunction, variable: Int64, mut budget: Budget) raises PolyError -> RationalFunction:
    """Compose univariate rational functions through homogeneous evaluation.

    Args:
        outer: The outer rational function.
        inner: The replacement rational function in the same order and variable.
        variable: The sole polynomial variable.
        budget: The shared work allowance.

    Returns:
        The normalized composition.

    Raises:
        PolyError: DOMAIN for incompatible inputs or an identically zero denominator, LIMIT or CAPACITY on work limits.
    """
    var f = _canonical(outer, budget)
    var g = _canonical(inner, budget)
    f.numerator._same(g.numerator, budget)
    var a = _univariate(f, variable, budget)
    var b = _univariate(g, variable, budget)
    var degree = max(a[0].degree(), a[1].degree())
    var numerator = _homogenize(a[0], b[0], b[1], degree, budget)
    var denominator = _homogenize(a[1], b[0], b[1], degree, budget)
    return _normalized(numerator.to_sparse(f.numerator.order, budget), denominator.to_sparse(f.numerator.order, budget), budget)


@fieldwise_init
struct PartialFractionTerm(Copyable, Movable):
    """A proper numerator divided by a positive power of a pole factor."""
    var numerator: DensePoly[QQ]
    """The numerator of degree less than the base factor."""
    var factor: DensePoly[QQ]
    """The monic base denominator factor."""
    var exponent: Int
    """The positive denominator power."""


@fieldwise_init
struct PartialFractions(Copyable, Movable):
    """An exact polynomial part and a list of proper fraction terms."""
    var polynomial: DensePoly[QQ]
    """The polynomial quotient."""
    var terms: List[PartialFractionTerm]
    """Nonzero proper terms grouped by deterministic factor order."""


def partial_fractions(f: RationalFunction, variable: Int64, mut budget: Budget) raises PolyError -> PartialFractions:
    """Decompose a univariate rational function over QQ.

    Args:
        f: The rational function.
        variable: Its sole variable.
        budget: The shared work allowance.

    Returns:
        A polynomial part and proper terms over irreducible rational factors.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    var reduced = _canonical(f, budget)
    var pair = _univariate(reduced, variable, budget)
    var division = divmod(pair[0], pair[1], budget)
    var factors = factor(pair[1], budget)
    var terms = List[PartialFractionTerm]()
    for entry in factors.factors:
        var power = entry[0].pow(UInt64(entry[1]), budget)
        var rest = _exact(pair[1], power, budget)
        var bezout = _extended(rest, power, budget)
        var local = _remainder(division.remainder.mul(bezout[1], budget), power, budget)
        for exponent in range(entry[1], 0, -1):
            var digit = divmod(local, entry[0], budget)
            if not digit.remainder.is_zero():
                _allocation[PartialFractionTerm](1, budget)
                terms.append(PartialFractionTerm(_copy(digit.remainder, budget), _copy(entry[0], budget), exponent))
            swap(local, digit.quotient)
    var polynomial = DensePoly[QQ](variable, NoContext())
    swap(polynomial, division.quotient)
    return PartialFractions(polynomial^, terms^)


def check_partial_fractions(f: RationalFunction, variable: Int64, decomposition: PartialFractions, mut budget: Budget) raises PolyError -> Bool:
    """Check a partial-fraction identity without repeating factorization.

    Args:
        f: The original rational function.
        variable: Its sole variable.
        decomposition: The proposed polynomial and proper terms.
        budget: The shared work allowance.

    Returns:
        Whether all terms are proper and recomposition equals f; irreducibility is not asserted by this identity check.

    Raises:
        PolyError: DOMAIN for malformed inputs, LIMIT or CAPACITY on work limits.
    """
    var pair = _univariate(f, variable, budget)
    decomposition.polynomial.check_invariants(budget)
    if decomposition.polynomial.variable != variable:
        return False
    var total = decomposition.polynomial.mul(pair[1], budget)
    for term in decomposition.terms:
        term.factor.check_invariants(budget)
        term.numerator.check_invariants(budget)
        if term.exponent <= 0 or term.factor.variable != variable or term.numerator.variable != variable or term.factor.degree() < 1 or term.numerator.degree() >= term.factor.degree():
            return False
        _charge_size(term.factor.coefficients[term.factor.degree()].size(), budget)
        if term.factor.coefficients[term.factor.degree()] != QQ(Rational(1)):
            return False
        var power = term.factor.pow(UInt64(term.exponent), budget)
        var division = divmod(pair[1], power, budget)
        if not division.remainder.is_zero():
            return False
        total = total.add(term.numerator.mul(division.quotient, budget), budget)
    for coefficient_index in range(len(total.coefficients)):
        ref coefficient = total.coefficients[coefficient_index]
        _charge_size(coefficient.size(), budget)
    for coefficient_index in range(len(pair[0].coefficients)):
        ref coefficient = pair[0].coefficients[coefficient_index]
        _charge_size(coefficient.size(), budget)
    return total == pair[0]


def rational_poles(f: RationalFunction, variable: Int64, mut budget: Budget) raises PolyError -> Factorization[QQ]:
    """Represent all complex poles by exact rational factors and multiplicities.

    Args:
        f: The rational function.
        variable: Its sole variable.
        budget: The shared work allowance.

    Returns:
        The factorization of its reduced denominator; each factor represents all its complex roots.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    var reduced = _canonical(f, budget)
    var pair = _univariate(reduced, variable, budget)
    return factor(pair[1], budget)


@fieldwise_init
struct LaurentExpansion(Copyable, Movable):
    """Consecutive exact coefficients in powers of x-minus-center."""
    var variable: Int64
    """The polynomial variable."""
    var center: QQ
    """The rational expansion center."""
    var valuation: Int
    """The lowest exponent; zero for the identically zero function."""
    var coefficients: List[QQ]
    """Coefficients beginning at the valuation, with a formal unknown tail."""


def _local(f: RationalFunction, variable: Int64, center: QQ, mut budget: Budget) raises PolyError -> Tuple[List[QQ], List[QQ], Int]:
    var pair = _univariate(f, variable, budget)
    var numerator = taylor_shift(pair[0], center, budget)
    var denominator = taylor_shift(pair[1], center, budget)
    var ni = 0
    var di = 0
    while ni < len(numerator.coefficients):
        _charge_size(numerator.coefficients[ni].size(), budget)
        if not numerator.coefficients[ni].is_zero():
            break
        ni += 1
    while di < len(denominator.coefficients):
        _charge_size(denominator.coefficients[di].size(), budget)
        if not denominator.coefficients[di].is_zero():
            break
        di += 1
    var valuation = ni - di if not numerator.is_zero() else 0
    return (_slice(_view(numerator.coefficients), ni, len(numerator.coefficients), budget), _slice(_view(denominator.coefficients), di, len(denominator.coefficients), budget), valuation)


def laurent_expansion(f: RationalFunction, variable: Int64, center: QQ, terms: Int, mut budget: Budget) raises PolyError -> LaurentExpansion:
    """Expand a rational function formally at an exact rational center.

    Args:
        f: The rational function.
        variable: Its sole variable.
        center: The expansion center, possibly a pole.
        terms: The nonnegative number of coefficients beginning at the valuation.
        budget: The shared work allowance.

    Returns:
        The exact Laurent coefficients and their lowest exponent.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    var local = _local(f, variable, center, budget)
    var valuation = local[2]
    var local_a = List[QQ]()
    var local_b = List[QQ]()
    swap(local_a, local[0])
    swap(local_b, local[1])
    var coefficients = _divide_values[_ExactOps[QQ]](_view(local_a), _view(local_b), terms, NoContext(), budget)
    _charge_size(center.size(), budget)
    return LaurentExpansion(variable, center.copy(), valuation, coefficients^)


def check_laurent_expansion(f: RationalFunction, expansion: LaurentExpansion, mut budget: Budget) raises PolyError -> Bool:
    """Check Laurent coefficients by a truncated multiplication identity.

    Args:
        f: The rational function.
        expansion: The proposed local expansion.
        budget: The shared work allowance.

    Returns:
        Whether the valuation and every supplied coefficient are correct.

    Raises:
        PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
    """
    var local = _local(f, expansion.variable, expansion.center, budget)
    var valuation = local[2]
    var local_a = List[QQ]()
    var local_b = List[QQ]()
    swap(local_a, local[0])
    swap(local_b, local[1])
    if valuation != expansion.valuation:
        return False
    var length = len(expansion.coefficients)
    var product = _mul_values[_ExactOps[QQ]](_view(local_b), _view(expansion.coefficients), length, NoContext(), budget)
    var actual = _values[_ExactOps[QQ]](_view(product), length, NoContext(), budget)
    var expected = _values[_ExactOps[QQ]](_view(local_a), length, NoContext(), budget)
    for i in range(length):
        _charge_size(_sum(actual[i].size(), expected[i].size()), budget)
        if actual[i] != expected[i]:
            return False
    return True


def rational_residue(f: RationalFunction, variable: Int64, center: QQ, mut budget: Budget) raises PolyError -> QQ:
    """Extract the coefficient of (x-center)^(-1) exactly.

    Args:
        f: The rational function.
        variable: Its sole variable.
        center: The exact rational expansion center.
        budget: The shared work allowance.

    Returns:
        The rational residue, including zero at regular points.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    var local = _local(f, variable, center, budget)
    var valuation = local[2]
    var local_a = List[QQ]()
    var local_b = List[QQ]()
    swap(local_a, local[0])
    swap(local_b, local[1])
    if valuation >= 0:
        return QQ(Rational(0))
    var coefficients = _divide_values[_ExactOps[QQ]](_view(local_a), _view(local_b), -valuation, NoContext(), budget)
    return coefficients.pop()
