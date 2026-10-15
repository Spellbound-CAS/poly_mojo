"""Exact-only dispatch and Newton updates for formal logarithms and exponentials."""
from .coefficients import CoefficientField, ModInt
from .coefficient_buffer import _CoefficientView, _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from ._coefficient_ops import _copy
from ._toolkit_exact import _ExactOps, _use_newton
from ._toolkit_ops import _mul_values
from ._series_ops import _log_values, _exp_values, _series_constant


def _exact_log[C: CoefficientField](a: _CoefficientView[C, _], length: Int, context: C.Context, mut budget: Budget) raises PolyError -> List[C]:
    if _use_newton[C](max(0, length - 1), len(a)):
        return _log_values[_ExactOps[C], True](_view(a), length, context, budget)
    return _log_values[_ExactOps[C]](_view(a), length, context, budget)


def _exact_exp[C: CoefficientField](a: _CoefficientView[C, _], length: Int, context: C.Context, mut budget: Budget) raises PolyError -> List[C]:
    comptime threshold = 1024 if C == ModInt else 128
    if length < threshold or len(a) < 16:
        return _exp_values[_ExactOps[C]](_view(a), length, context, budget)
    return _newton_exp(_view(a), length, context, budget)


def _newton_exp[C: CoefficientField](a: _CoefficientView[C, _], length: Int, context: C.Context, mut budget: Budget) raises PolyError -> List[C]:
    if length < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "exponential needs a known constant coefficient")
    _series_constant[_ExactOps[C]](_view(a), 0, context, budget)
    var known = min(length, 16)
    var result = _exp_values[_ExactOps[C]](_view(a), known, context, budget)
    while known < length:
        var next = _sum(known, min(known, length - known))
        # g_new = g * (1 + a - log(g)). The error starts at x^known,
        # so only its high block and a truncated product are needed.
        # _exact_log checks every required integer divisor, including when
        # the numerator vanishes in a field of positive characteristic.
        var logarithm = _exact_log(_view(result), next, context, budget)
        _allocation[C](next - known, budget)
        var high = List[C](capacity=next - known)
        for k in range(known, next):
            var value = _copy(a[k], budget) if k < len(a) else C.zero(context)
            high.append(value.sub(logarithm[k], budget))
        var correction = _mul_values[_ExactOps[C]](_view(result), _view(high), next - known, context, budget)
        _allocation[C](next - known, budget)
        result.reserve(next)
        for k in range(next - known):
            budget.charge(1)
            result.append(_copy(correction[k], budget) if k < len(correction) else C.zero(context))
        known = next
    return result^
