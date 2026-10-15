"""Exact rectangular arithmetic for numerical division witnesses."""
from apn_mojo import Rational
from .coefficients import QQ
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .numeric_bounds import NumericBox, _box_charge, _box_zero, _box_real, _box_mul, _box_sub, _point
from ._coefficient_ops import _copy
from ._root_ops import _compare


def _copy_boxes(values: List[NumericBox], mut budget: Budget) raises PolyError -> List[NumericBox]:
    _allocation[NumericBox](len(values), budget)
    var result = List[NumericBox](capacity=len(values))
    for box in values:
        _box_charge(box, budget)
        result.append(box.copy())
    return result^


def _square_interval(low: QQ, high: QQ, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    var a = low.mul(low, budget)
    var b = high.mul(high, budget)
    var minimum = QQ(Rational(0))
    if low.value.sign() > 0 or high.value.sign() < 0:
        minimum = _copy(a, budget) if _compare(a, b, budget) < 0 else _copy(b, budget)
    var maximum = _copy(a, budget) if _compare(a, b, budget) > 0 else _copy(b, budget)
    return (minimum^, maximum^)


def _box_divide(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    if _box_zero(b, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "division enclosure contains zero")
    if _box_real(b, budget):
        var inverse = NumericBox.real_interval(b.real_high.inverse(budget), b.real_low.inverse(budget), budget)
        return _box_mul(a, inverse, budget)
    var real = _square_interval(b.real_low, b.real_high, budget)
    var imag = _square_interval(b.imag_low, b.imag_high, budget)
    var norm_low = real[0].add(imag[0], budget)
    var norm_high = real[1].add(imag[1], budget)
    var conjugate = NumericBox(_copy(b.real_low, budget), _copy(b.real_high, budget), b.imag_high.neg(budget), b.imag_low.neg(budget))
    var numerator = _box_mul(a, conjugate, budget)
    var inverse = NumericBox.real_interval(norm_high.inverse(budget), norm_low.inverse(budget), budget)
    return _box_mul(numerator, inverse, budget)


def _residual_boxes(a: List[NumericBox], b: List[NumericBox], q: List[NumericBox], r: List[NumericBox], mut budget: Budget) raises PolyError -> List[NumericBox]:
    var product_count = _sum(len(q), len(b)) - 1 if len(q) and len(b) else 0
    var count = max(len(a), max(len(r), product_count))
    _allocation[NumericBox](count, budget)
    var result = List[NumericBox](length=count, fill=_point(0))
    for i in range(count):
        if i < len(a):
            _box_charge(a[i], budget)
            result[i] = a[i].copy()
        if i < len(r):
            result[i] = _box_sub(result[i], r[i], budget)
    for i in range(len(q)):
        for j in range(len(b)):
            result[i + j] = _box_sub(result[i + j], _box_mul(q[i], b[j], budget), budget)
    return result^
