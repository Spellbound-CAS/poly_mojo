"""Call-scoped field modulus and in-place reduction of owned products."""
from .coefficients import CoefficientField
from .dense import DensePoly, _slice, _trim
from .coefficient_buffer import CoefficientBuffer, _view
from .errors import Budget, PolyError, PolyErrorKind


@fieldwise_init
struct _DenseReducer[C: CoefficientField](Movable):
    var modulus: DensePoly[Self.C]
    var inverse: Self.C

    @staticmethod
    def _reduce(modulus: DensePoly[Self.C], inverse: Self.C, mut coefficients: CoefficientBuffer[Self.C], mut budget: Budget) raises PolyError:
        var size = len(modulus.coefficients)
        while len(coefficients) >= size:
            budget.charge(1)
            var shift = len(coefficients) - size
            var factor = coefficients[len(coefficients) - 1].mul(inverse, budget)
            for i in range(size):
                budget.charge(1)
                coefficients[shift + i] = coefficients[shift + i].sub(factor.mul(modulus.coefficients[i], budget), budget)
            _trim(coefficients, budget)

    def remainder(self, var f: DensePoly[Self.C], mut budget: Budget) raises PolyError -> DensePoly[Self.C]:
        self.modulus._same(f)
        Self._reduce(self.modulus, self.inverse, f.coefficients, budget)
        return f^

    def product(self, a: DensePoly[Self.C], b: DensePoly[Self.C], mut budget: Budget) raises PolyError -> DensePoly[Self.C]:
        # Multiplication already returns an owned vector. Reducing it directly
        # avoids both scratch initialization and copying the bounded remainder.
        return self.remainder(a.mul(b, budget), budget)


def _dense_reducer[C: CoefficientField](g: DensePoly[C], mut budget: Budget) raises PolyError -> _DenseReducer[C]:
    g.check_invariants(budget)
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "reduction by zero polynomial")
    var inverse = g.coefficients[g.degree()].inverse(budget)
    var modulus = DensePoly[C]._from_canonical(g.variable, g.ring, _slice(_view(g.coefficients), 0, len(g.coefficients), budget))
    return _DenseReducer[C](modulus^, inverse^)
