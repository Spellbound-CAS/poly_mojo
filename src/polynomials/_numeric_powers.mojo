"""Call-local shared binary squares preserving each original power schedule."""
from std.collections import Dict
from .errors import Budget, PolyError, _allocation, _charge_size
from .numeric_domains import NumericDomain, _numeric_integer, _numeric_apply, _numeric_cost


struct _NumericPowerCache[D: NumericDomain](Movable):
    var _indices: Dict[Tuple[Int, Int], Int]
    var _values: List[Self.D.Scalar]

    def __init__(out self):
        self._indices = Dict[Tuple[Int, Int], Int]()
        self._values = List[Self.D.Scalar]()

    def power(mut self, assigned: Int, value: Self.D.Scalar, exponent: UInt64, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        var result = _numeric_integer[Self.D](1, context, budget)
        var base = value.copy()
        var n = exponent
        var bit = 0
        while n:
            budget.charge(1)
            if n & 1:
                result = _numeric_apply[Self.D](result, base, 2, context, budget)
            n >>= 1
            bit += 1
            if n:
                budget.charge(1)
                var key = (assigned, bit)
                var slot = self._indices.get(key, -1)
                if slot < 0:
                    var square = _numeric_apply[Self.D](base, base, 2, context, budget)
                    _allocation[Tuple[Int, Int]](1, budget)
                    _allocation[Self.D.Scalar](1, budget)
                    self._indices[key] = len(self._values)
                    self._values.append(square.copy())
                    base = square^
                else:
                    _charge_size(_numeric_cost[Self.D](self._values[slot]), budget)
                    base = self._values[slot].copy()
        return result^
