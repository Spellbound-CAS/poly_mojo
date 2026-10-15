"""Exact polynomial calculus, affine transformations and simultaneous substitution."""
from apn_mojo import Rational
from .coefficients import CoefficientRing, QQ, ZZ
from .dense import DensePoly
from .sparse import SparsePoly, Term
from .monomial import Monomial
from ._assignment import _AssignmentPlan
from ._bounded_products import _PolynomialSum
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum, _product
from .evaluate import evaluate_dense, _derivative_dense, _natural
from ._coefficient_ops import _exact, _copy
from ._images import _rational_image
from ._toolkit_exact import _ExactOps
from ._toolkit_ops import _integral_values, _compose_values, _jet_values


def derivative_dense[C: CoefficientRing](poly: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    """Differentiate a dense polynomial exactly.

    Parameters:
        C: The exact coefficient ring.

    Args:
        poly: The input polynomial.
        budget: The shared work allowance.

    Returns:
        The formal derivative, retaining the variable and coefficient ring.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    return _derivative_dense(poly, budget)


def rational_dense(poly: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    """Promote an integer polynomial explicitly to rational coefficients.

    Args:
        poly: The integer polynomial.
        budget: The shared work allowance.

    Returns:
        The same polynomial over QQ.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    _allocation[QQ](len(poly.coefficients), budget)
    var values = List[QQ](capacity=len(poly.coefficients))
    for coefficient_index in range(len(poly.coefficients)):
        ref value = poly.coefficients[coefficient_index]
        _charge_size(value.size(), budget)
        values.append(QQ(Rational(value.value)))
    return DensePoly[QQ].from_coefficients(poly.variable, poly.ring, values^, budget)


def rational_sparse(poly: SparsePoly[ZZ], mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    """Promote a sparse integer polynomial explicitly to QQ.

    Args:
        poly: The integer polynomial.
        budget: The shared work allowance.

    Returns:
        The same polynomial with rational coefficients.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    return _rational_image(poly, budget)


def integrate_dense[C: CoefficientRing](poly: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    """Integrate exactly with zero constant when coefficient divisions exist.

    Parameters:
        C: An exact ring supporting the required divisions.

    Args:
        poly: The input polynomial.
        budget: The shared work allowance.

    Returns:
        The zero-constant antiderivative in the same ring.

    Raises:
        PolyError: DOMAIN for a nonintegrable term, UNSUPPORTED for missing exact division, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    var values = _integral_values[_ExactOps[C]](_view(poly.coefficients), poly.ring, budget)
    return DensePoly[C].from_coefficients(poly.variable, poly.ring, values^, budget)


def integrate_integer_dense(poly: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    """Integrate an integer polynomial after explicit rational promotion.

    Args:
        poly: The integer polynomial.
        budget: The shared work allowance.

    Returns:
        Its rational antiderivative with zero constant.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    return integrate_dense(rational_dense(poly, budget), budget)


def definite_integral[C: CoefficientRing](poly: DensePoly[C], low: C, high: C, mut budget: Budget) raises PolyError -> C:
    """Evaluate an exact antiderivative difference at two endpoints.

    Parameters:
        C: An exact ring supporting the required divisions.

    Args:
        poly: The polynomial.
        low: The lower endpoint in its coefficient ring.
        high: The upper endpoint in its coefficient ring.
        budget: The shared work allowance.

    Returns:
        F(high) minus F(low), where F is an antiderivative.

    Raises:
        PolyError: DOMAIN for incompatible endpoints or nonintegrable terms, LIMIT or CAPACITY on work limits.
    """
    if low.context() != poly.ring or high.context() != poly.ring:
        raise PolyError(PolyErrorKind.DOMAIN, "integration endpoints in another ring")
    var integral = integrate_dense(poly, budget)
    return evaluate_dense(integral, high, budget).sub(evaluate_dense(integral, low, budget), budget)


def compose_dense[C: CoefficientRing](outer: DensePoly[C], inner: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    """Substitute a dense polynomial into another in the same variable.

    Parameters:
        C: The exact coefficient ring.

    Args:
        outer: The outer polynomial.
        inner: The replacement polynomial.
        budget: The shared work allowance.

    Returns:
        The exact composition outer(inner).

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    outer.check_invariants(budget)
    inner.check_invariants(budget)
    outer._same(inner)
    var count = _sum(_product(max(0, outer.degree()), max(0, inner.degree())), 1)
    var values = _compose_values[_ExactOps[C]](_view(outer.coefficients), _view(inner.coefficients), count, outer.ring, budget)
    return DensePoly[C].from_coefficients(outer.variable, outer.ring, values^, budget)


def affine_transform[C: CoefficientRing](poly: DensePoly[C], scale: C, shift: C, mut budget: Budget) raises PolyError -> DensePoly[C]:
    """Compute p(scale*x + shift) exactly.

    Parameters:
        C: The exact coefficient ring.

    Args:
        poly: The polynomial.
        scale: The variable scale in the coefficient ring.
        shift: The translation in the coefficient ring.
        budget: The shared work allowance.

    Returns:
        The transformed polynomial.

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    var inner = DensePoly[C].from_coefficients(poly.variable, poly.ring, [_copy(shift, budget), _copy(scale, budget)], budget)
    return compose_dense(poly, inner, budget)


def taylor_shift[C: CoefficientRing](poly: DensePoly[C], point: C, mut budget: Budget) raises PolyError -> DensePoly[C]:
    """Compute p(x + point) without factorial division.

    Parameters:
        C: The exact coefficient ring.

    Args:
        poly: The polynomial.
        point: The shift in its coefficient ring.
        budget: The shared work allowance.

    Returns:
        The Taylor polynomial in the shifted variable.

    Raises:
        PolyError: DOMAIN for incompatible input, LIMIT or CAPACITY on work limits.
    """
    var values = jet_evaluate(poly, point, len(poly.coefficients), budget)
    return DensePoly[C].from_coefficients(poly.variable, poly.ring, values^, budget)


def jet_evaluate[C: CoefficientRing](poly: DensePoly[C], point: C, order: Int, mut budget: Budget) raises PolyError -> List[C]:
    """Evaluate Taylor coefficients through degree order minus one.

    Parameters:
        C: The exact coefficient ring.

    Args:
        poly: The polynomial.
        point: The expansion point.
        order: The nonnegative number of coefficients requested.
        budget: The shared work allowance.

    Returns:
        Hasse derivatives, equal to p^(k)(point)/k! in characteristic zero.

    Raises:
        PolyError: DOMAIN for a negative order or mixed rings, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    if point.context() != poly.ring:
        raise PolyError(PolyErrorKind.DOMAIN, "Taylor point in another ring")
    return _jet_values[_ExactOps[C]](_view(poly.coefficients), point, order, poly.ring, budget)


def integrate_sparse[C: CoefficientRing](poly: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Integrate one sparse variable with zero integration constant.

    Parameters:
        C: An exact ring supporting the required coefficient divisions.

    Args:
        poly: The polynomial.
        variable: The integration variable, allowed by its monomial order.
        budget: The shared work allowance.

    Returns:
        An exact antiderivative in that variable.

    Raises:
        PolyError: DOMAIN for unavailable divisions or invalid variables, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    var unit = Monomial([(variable, UInt32(1))], budget)
    poly.order.validate(unit, budget)
    _allocation[Term[C]](len(poly.terms), budget)
    var terms = List[Term[C]](capacity=len(poly.terms))
    for term in poly.terms:
        var exponent = term.monomial.exponent(variable, budget)
        if exponent == UInt32.MAX:
            raise PolyError(PolyErrorKind.CAPACITY, "integration exponent exceeds UInt32")
        var coefficient = _exact(term.coefficient, _natural[C](UInt64(exponent) + 1, poly.ring, budget), budget)
        terms.append(Term[C](term.monomial.mul(unit, budget), coefficient^))
    return SparsePoly[C].from_terms(terms^, poly.order, poly.ring, budget)


def integrate_integer_sparse(poly: SparsePoly[ZZ], variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    """Integrate sparse integer coefficients after promotion to QQ.

    Args:
        poly: The integer polynomial.
        variable: The integration variable.
        budget: The shared work allowance.

    Returns:
        Its rational antiderivative with zero integration constant.

    Raises:
        PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
    """
    return integrate_sparse(rational_sparse(poly, budget), variable, budget)


def substitute_many[C: CoefficientRing](poly: SparsePoly[C], assignments: List[Tuple[Int64, SparsePoly[C]]], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Replace several variables simultaneously, retaining unassigned variables.

    Parameters:
        C: The exact coefficient ring.

    Args:
        poly: The polynomial.
        assignments: Distinct variable IDs and replacements in its ring and order.
        budget: The shared work allowance.

    Returns:
        The simultaneous substitution, with no substitutions inside replacements.

    Raises:
        PolyError: DOMAIN for duplicates or incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    var plan = _AssignmentPlan()
    for i in range(len(assignments)):
        poly._same(assignments[i][1], budget)
        assignments[i][1].check_invariants(budget)
        poly.order.validate(Monomial([(assignments[i][0], UInt32(1))], budget), budget)
        plan.assign(assignments[i][0], i, budget)
    var result = SparsePoly[C].zero(poly.order, poly.ring, budget)
    var powers = List[SparsePoly[C]]()
    var remaining = List[Int]()
    var prepared = List[Bool]()
    var bounded = len(poly.terms) >= 32
    var sums = _PolynomialSum[C]()
    if bounded:
        for term in poly.terms:
            for pair in term.monomial.pairs:
                var slot = plan.power(pair, budget)
                if slot == len(remaining):
                    _allocation[Int](1, budget)
                    _allocation[SparsePoly[C]](1, budget)
                    _allocation[Bool](1, budget)
                    prepared.append(False)
                    remaining.append(0)
                    powers.append(SparsePoly[C].zero(poly.order, poly.ring, budget))
                remaining[slot] = _sum(remaining[slot], 1)
    for term in poly.terms:
        var product = SparsePoly[C].constant(_copy(term.coefficient, budget), poly.order, budget)
        for pair in term.monomial.pairs:
            var slot = plan.power(pair, budget)
            if slot == len(powers) or (bounded and not prepared[slot]):
                var base = SparsePoly[C].variable(pair[0], poly.order, poly.ring, budget)
                var assigned = plan.assignment(pair[0], budget)
                if assigned >= 0:
                    base = assignments[assigned][1]._copy(budget)
                var power = base.pow(UInt64(pair[1]), budget)
                _allocation[SparsePoly[C]](1, budget)
                if slot == len(powers):
                    powers.append(power^)
                else:
                    powers[slot] = power^
                    prepared[slot] = True
            product = product.mul(powers[slot], budget)
            if bounded:
                remaining[slot] -= 1
                if remaining[slot] == 0:
                    powers[slot].terms = List[Term[C]]()
        if bounded:
            sums.add(product^, budget)
        else:
            result = result.add(product, budget)
    return sums.finish(poly, budget) if bounded else result^
