"""Explicit numerical scalar operations for shared coefficient-vector algorithms."""
from .numeric_bounds import numeric_bounds, _point, _box_charge
from ._toolkit_ops import _ScalarOps, _convolve_values
from .numeric_domains import NumericDomain, _numeric_integer, _numeric_apply, _numeric_cost, _numeric_zero, _context_cost
from .coefficient_buffer import _CoefficientView, _view
from .errors import Budget, PolyError, _charge_size


def _validate[D: NumericDomain](values: List[D.Scalar], context: D.Context, mut budget: Budget) raises PolyError:
    _ = _context_cost[D](context)
    for value in values:
        _charge_size(_numeric_cost[D](value), budget)


struct _NumericOps[D: NumericDomain](_ScalarOps):
    comptime Scalar = Self.D.Scalar
    comptime Context = Self.D.Context
    comptime EXACT = False

    @staticmethod
    def multiply(a: _CoefficientView[Self.Scalar, _], b: _CoefficientView[Self.Scalar, _], length: Int, context: Self.Context, mut budget: Budget) raises PolyError -> List[Self.Scalar]:
        return _convolve_values[Self](_view(a), _view(b), length, context, budget)

    @staticmethod
    def integer(value: Int, context: Self.Context, mut budget: Budget) raises PolyError -> Self.Scalar:
        return _numeric_integer[Self.D](value, context, budget)

    @staticmethod
    def apply(a: Self.Scalar, b: Self.Scalar, operation: Int, context: Self.Context, mut budget: Budget) raises PolyError -> Self.Scalar:
        return _numeric_apply[Self.D](a, b, operation, context, budget)

    @staticmethod
    def copy(value: Self.Scalar, mut budget: Budget) raises PolyError -> Self.Scalar:
        _charge_size(_numeric_cost[Self.D](value), budget)
        return value.copy()

    @staticmethod
    def matches_integer(value: Self.Scalar, integer: Int, context: Self.Context, mut budget: Budget) raises PolyError -> Bool:
        var bounds = numeric_bounds[Self.D](value, budget)
        _box_charge(bounds, budget)
        return bounds == _point(integer)

    @staticmethod
    def zero(value: Self.Scalar) -> Bool:
        return _numeric_zero[Self.D](value)
