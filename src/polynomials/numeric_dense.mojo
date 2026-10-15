"""Dense numerical polynomials with explicit precision and stable operation order."""
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _charge_size
from .numeric_domains import NumericDomain, _numeric_cost, _numeric_zero, _numeric_nonzero, _numeric_integer, _numeric_apply, _context_cost
from ._toolkit_numeric import _NumericOps
from ._toolkit_ops import _mul_values


struct NumericDensePoly[D: NumericDomain](Copyable, Movable):
    """A finite univariate numerical polynomial in ascending coefficient order.

    Parameters:
        D: An apn_mojo numerical domain adapter.
    """
    var variable: Int64
    """The variable identifier, also retained by constants."""
    var coefficients: List[Self.D.Scalar]
    """Ascending coefficients with only certain trailing zeros removed."""

    def __init__(out self, variable: Int64 = 0):
        """Construct the zero polynomial.

        Args:
            variable: The retained variable identifier.
        """
        self.variable = variable
        self.coefficients = List[Self.D.Scalar]()

    @staticmethod
    def from_coefficients(variable: Int64, var coefficients: List[Self.D.Scalar], mut budget: Budget) raises PolyError -> Self:
        """Validate finite coefficients and remove certain trailing zeros.

        Args:
            variable: The variable identifier.
            coefficients: The owned ascending coefficient vector.
            budget: The shared work allowance.

        Returns:
            A numerical polynomial retaining uncertain coefficients.

        Raises:
            PolyError: DOMAIN for nonfinite values, LIMIT or CAPACITY on work limits.
        """
        var result = Self(variable)
        for value in coefficients:
            _charge_size(_numeric_cost[Self.D](value), budget)
        result.coefficients = coefficients^
        while len(result.coefficients) and _numeric_zero[Self.D](result.coefficients[len(result.coefficients) - 1]):
            budget.charge(1)
            _ = result.coefficients.pop()
        return result^

    def degree_bound(self) -> Int:
        """Return the degree upper bound, or minus one for zero.

        Returns:
            The largest stored power.
        """
        return len(self.coefficients) - 1

    def certified_degree(self, mut budget: Budget) raises PolyError -> Optional[Int]:
        """Determine a degree shared by every represented polynomial.

        Args:
            budget: The shared work allowance.

        Returns:
            The degree, minus one for zero, or None for an uncertain leading term.

        Raises:
            PolyError: DOMAIN for invalid coefficients or LIMIT on exhaustion.
        """
        if not len(self.coefficients):
            return -1
        _charge_size(_numeric_cost[Self.D](self.coefficients[len(self.coefficients) - 1]), budget)
        return Optional(self.degree_bound()) if _numeric_nonzero[Self.D](self.coefficients[len(self.coefficients) - 1]) else None

    def _same(self, other: Self) raises PolyError:
        if self.variable != other.variable:
            raise PolyError(PolyErrorKind.DOMAIN, "different numerical polynomial variables")

    def _combine(self, other: Self, operation: Int, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        self._same(other)
        var count = max(len(self.coefficients), len(other.coefficients))
        _allocation[Self.D.Scalar](count, budget)
        var zero = _numeric_integer[Self.D](0, context, budget)
        var values = List[Self.D.Scalar](capacity=count)
        for i in range(count):
            if i < len(self.coefficients):
                _charge_size(_numeric_cost[Self.D](self.coefficients[i]), budget)
            if i < len(other.coefficients):
                _charge_size(_numeric_cost[Self.D](other.coefficients[i]), budget)
            var a = self.coefficients[i].copy() if i < len(self.coefficients) else zero.copy()
            var b = other.coefficients[i].copy() if i < len(other.coefficients) else zero.copy()
            values.append(_numeric_apply[Self.D](a, b, operation, context, budget))
        return Self.from_coefficients(self.variable, values^, budget)

    def add(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Add with the requested output precision.

        Args:
            other: The polynomial in the same variable.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing sum.

        Raises:
            PolyError: DOMAIN for incompatible or invalid values, LIMIT or CAPACITY on work limits.
        """
        return self._combine(other, 0, context, budget)

    def sub(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Subtract with the requested output precision.

        Args:
            other: The polynomial in the same variable.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing difference.

        Raises:
            PolyError: DOMAIN for incompatible or invalid values, LIMIT or CAPACITY on work limits.
        """
        return self._combine(other, 1, context, budget)

    def neg(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Negate every coefficient at the requested precision.

        Args:
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing negation.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        return Self(self.variable)._combine(self, 1, context, budget)

    def mul(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Multiply by convolution in increasing left and right index order.

        Args:
            other: The polynomial in the same variable.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing product.

        Raises:
            PolyError: DOMAIN for incompatible or invalid values, LIMIT or CAPACITY on work limits.
        """
        self._same(other)
        _ = _context_cost[Self.D](context)
        if not len(self.coefficients) or not len(other.coefficients):
            return Self(self.variable)
        var count = _sum(len(self.coefficients), len(other.coefficients)) - 1
        var values = _mul_values[_NumericOps[Self.D]](_view(self.coefficients), _view(other.coefficients), count, context, budget)
        return Self.from_coefficients(self.variable, values^, budget)

    def scale(self, scalar: Self.D.Scalar, context: Self.D.Context, mut budget: Budget, *, divide: Bool = False) raises PolyError -> Self:
        """Multiply coefficients by a scalar or divide by a certified nonzero scalar.

        Args:
            scalar: The finite scalar or bounded enclosure.
            context: Explicit output settings.
            budget: The shared work allowance.
            divide: Select division, requiring the scalar to exclude zero.

        Returns:
            The rounded or enclosing scaled polynomial.

        Raises:
            PolyError: DOMAIN for invalid or uncertain divisors, LIMIT or CAPACITY on work limits.
        """
        _charge_size(_numeric_cost[Self.D](scalar), budget)
        _ = _context_cost[Self.D](context)
        if divide and not _numeric_nonzero[Self.D](scalar):
            raise PolyError(PolyErrorKind.DOMAIN, "scalar divisor is not certified nonzero")
        _allocation[Self.D.Scalar](len(self.coefficients), budget)
        var values = List[Self.D.Scalar](capacity=len(self.coefficients))
        for a in self.coefficients:
            values.append(_numeric_apply[Self.D](a, scalar, 3 if divide else 2, context, budget))
        return Self.from_coefficients(self.variable, values^, budget)

    def pow(self, exponent: UInt64, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Compute a nonnegative power by deterministic binary powering.

        Args:
            exponent: The nonnegative power.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing power, including one for exponent zero.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        var result = Self.from_coefficients(self.variable, [_numeric_integer[Self.D](1, context, budget)], budget)
        _charge_size(len(self.coefficients), budget)
        for coefficient in self.coefficients:
            _charge_size(_numeric_cost[Self.D](coefficient), budget)
        var base = self.copy()
        var n = exponent
        while n:
            budget.charge(1)
            if n & 1:
                result = result.mul(base, context, budget)
            n >>= 1
            if n:
                base = base.mul(base, context, budget)
        return result^

    def derivative(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Differentiate with coefficient products rounded at the requested precision.

        Args:
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing derivative.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        _ = _context_cost[Self.D](context)
        _allocation[Self.D.Scalar](max(0, len(self.coefficients) - 1), budget)
        var values = List[Self.D.Scalar]()
        for i in range(1, len(self.coefficients)):
            var factor = _numeric_integer[Self.D](i, context, budget)
            values.append(_numeric_apply[Self.D](self.coefficients[i], factor, 2, context, budget))
        return Self.from_coefficients(self.variable, values^, budget)

    def evaluate(self, value: Self.D.Scalar, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Evaluate by descending Horner steps without implicit reassociation.

        Args:
            value: The finite point or bounded argument enclosure.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded value or rigorous enclosure.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        _charge_size(_numeric_cost[Self.D](value), budget)
        var result = _numeric_integer[Self.D](0, context, budget)
        for i in range(len(self.coefficients) - 1, -1, -1):
            result = _numeric_apply[Self.D](result, value, 2, context, budget)
            result = _numeric_apply[Self.D](result, self.coefficients[i], 0, context, budget)
        return result^

    def evaluate_with_derivative(self, value: Self.D.Scalar, context: Self.D.Context, mut budget: Budget) raises PolyError -> Tuple[Self.D.Scalar, Self.D.Scalar]:
        """Evaluate the value and derivative together by extended Horner steps.

        Args:
            value: The finite point or bounded argument enclosure.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The value and derivative, each rounded or enclosed.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        _charge_size(_numeric_cost[Self.D](value), budget)
        var result = _numeric_integer[Self.D](0, context, budget)
        var slope = _numeric_integer[Self.D](0, context, budget)
        for i in range(len(self.coefficients) - 1, -1, -1):
            slope = _numeric_apply[Self.D](_numeric_apply[Self.D](slope, value, 2, context, budget), result, 0, context, budget)
            result = _numeric_apply[Self.D](_numeric_apply[Self.D](result, value, 2, context, budget), self.coefficients[i], 0, context, budget)
        return (result^, slope^)
