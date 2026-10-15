"""Call-scoped exact powers, sharing squares and repeated exponents."""
from .coefficients import CoefficientRing
from ._coefficient_ops import _copy
from .errors import Budget, PolyError, _allocation, _charge_size


@fieldwise_init
struct _Powers[C: CoefficientRing](Movable):
    var squares: List[Self.C]
    var exponents: List[UInt64]
    var values: List[Self.C]

    def power(mut self, exponent: UInt64, mut budget: Budget) raises PolyError -> Self.C:
        var lo = 0
        var hi = len(self.exponents)
        while lo < hi:
            budget.charge(1)
            var mid = lo + (hi - lo) // 2
            if self.exponents[mid] < exponent:
                lo = mid + 1
            else:
                hi = mid
        if lo < len(self.exponents) and self.exponents[lo] == exponent:
            return _copy(self.values[lo], budget)
        var result = Self.C.one(self.squares[0].context())
        var remaining = exponent
        var bit = 0
        while remaining:
            budget.charge(1)
            if bit == len(self.squares):
                var square = self.squares[bit - 1].mul(self.squares[bit - 1], budget)
                _allocation[Self.C](1, budget)
                self.squares.append(square^)
            if remaining & 1:
                result = result.mul(self.squares[bit], budget)
            remaining >>= 1
            bit += 1
        var saved = _copy(result, budget)
        _allocation[UInt64](1, budget)
        _allocation[Self.C](1, budget)
        _charge_size(len(self.exponents) - lo, budget)
        _charge_size(len(self.values) - lo, budget)
        # Every raising operation precedes the paired cache mutation.
        self.exponents.insert(lo, exponent)
        self.values.insert(lo, saved^)
        return result^


def _powers[C: CoefficientRing](base: C, mut budget: Budget) raises PolyError -> _Powers[C]:
    _allocation[C](1, budget)
    var squares = List[C](capacity=1)
    squares.append(_copy(base, budget))
    return _Powers[C](squares^, List[UInt64](), List[C]())
