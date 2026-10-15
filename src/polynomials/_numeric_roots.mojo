"""Shared exact bounds, derivative families and resource handling for root searches."""
from apn_mojo import Rational
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from ._root_ops import _compare
from .numeric_bounds import NumericBox, _point, _box_valid, _box_zero, _box_real, _box_point, _box_mul, _box_add
from .numeric_options import NumericOutcome, _norm_lower, _norm_upper


def _root_degree(coefficients: List[NumericBox], mut budget: Budget) raises PolyError -> Optional[Int]:
    var degree = len(coefficients) - 1
    for box in coefficients:
        if not _box_valid(box, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid root coefficient enclosure")
    while degree >= 0 and coefficients[degree] == _point(0):
        budget.charge(1)
        degree -= 1
    if degree < 0:
        raise PolyError(PolyErrorKind.DOMAIN, "zero polynomial has no finite root isolation")
    if _box_zero(coefficients[degree], budget):
        return None
    return degree


def _root_bound(coefficients: List[NumericBox], degree: Int, mut budget: Budget) raises PolyError -> QQ:
    var leading = _norm_lower(coefficients[degree], budget)
    if leading.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "root bound requires a uniformly nonzero leading coefficient")
    var inverse = leading.inverse(budget)
    var maximum = QQ(Rational(0))
    for i in range(degree):
        var ratio = _norm_upper(coefficients[i], budget).mul(inverse, budget)
        if _compare(ratio, maximum, budget) > 0:
            maximum = ratio^
    return maximum.add(QQ(Rational(2)), budget)


def _derivative_boxes(coefficients: List[NumericBox], mut budget: Budget) raises PolyError -> List[NumericBox]:
    _allocation[NumericBox](max(0, len(coefficients) - 1), budget)
    var values = List[NumericBox]()
    for i in range(1, len(coefficients)):
        values.append(_box_mul(coefficients[i], _point(i), budget))
    return values^


def _rational_polynomial(coefficients: List[NumericBox], mut budget: Budget) raises PolyError -> Optional[DensePoly[QQ]]:
    _allocation[QQ](len(coefficients), budget)
    var values = List[QQ]()
    for box in coefficients:
        if not _box_real(box, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "real roots require real coefficients")
        if not _box_point(box, budget):
            return None
        _charge_size(box.real_low.size(), budget)
        values.append(box.real_low.copy())
    return DensePoly[QQ].from_coefficients(0, NoContext(), values^, budget)


def _exact_box_evaluate(coefficients: List[NumericBox], argument: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    var result = _point(0)
    for i in range(len(coefficients) - 1, -1, -1):
        result = _box_add(_box_mul(result, argument, budget), coefficients[i], budget)
    return result^


def _translated_boxes(coefficients: List[NumericBox], center: NumericBox, mut budget: Budget) raises PolyError -> List[NumericBox]:
    # Horner in (w + center), retaining exact rational coefficient enclosures.
    _allocation[NumericBox](len(coefficients), budget)
    var result = List[NumericBox]()
    for i in range(len(coefficients) - 1, -1, -1):
        result.append(_point(0))
        for j in range(len(result) - 1, 0, -1):
            result[j] = _box_add(_box_mul(result[j], center, budget), result[j - 1], budget)
        result[0] = _box_add(_box_mul(result[0], center, budget), coefficients[i], budget)
    return result^


def _strict_inside(inner: NumericBox, outer: NumericBox, mut budget: Budget, real: Bool = False) raises PolyError -> Bool:
    if not _box_valid(inner, budget) or not _box_valid(outer, budget):
        return False
    if _compare(outer.real_low, inner.real_low, budget) >= 0 or _compare(inner.real_high, outer.real_high, budget) >= 0:
        return False
    return real or (_compare(outer.imag_low, inner.imag_low, budget) < 0 and _compare(inner.imag_high, outer.imag_high, budget) < 0)


def _disjoint(a: NumericBox, b: NumericBox, mut budget: Budget, touching: Bool = False) raises PolyError -> Bool:
    var threshold = 0 if touching else -1
    return _compare(a.real_high, b.real_low, budget) <= threshold or _compare(b.real_high, a.real_low, budget) <= threshold or _compare(a.imag_high, b.imag_low, budget) <= threshold or _compare(b.imag_high, a.imag_low, budget) <= threshold


def _require_success(outcome: NumericOutcome, required: Bool) raises PolyError:
    if required and outcome != NumericOutcome.SUCCESS:
        var kind = PolyErrorKind.UNSUPPORTED if outcome == NumericOutcome.UNCERTAIN_DEGREE or outcome == NumericOutcome.INPUT_UNCERTAINTY else PolyErrorKind.LIMIT
        raise PolyError(kind, "numerical root request is inconclusive")
