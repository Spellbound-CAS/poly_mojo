"""Formal truncated series with explicit precision and domain conditions."""
from .numeric_domains import NumericDomain
from .numeric_dense import NumericDensePoly
from ._toolkit_numeric import _NumericOps, _validate
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _sum
from ._toolkit_ops import _values, _add_values, _mul_values, _divide_values, _compose_values, _derivative_values, _integral_values
from ._series_ops import _series_shape, _series_same_order, _series_truncation, _series_composition, _log_values, _exp_values, _revert_values


@fieldwise_init
struct NumericSeries[D: NumericDomain](Copyable, Movable):
    """Store coefficients known modulo an explicit power of x.

    Parameters:
        D: The numerical domain; enclosures describe coefficients, not analytic tails.
    """
    var polynomial: NumericDensePoly[Self.D]
    """The canonical polynomial containing known coefficients."""
    var order: Int
    """The nonnegative exponent of the unknown remainder O(x^order)."""

    @staticmethod
    def from_polynomial(poly: NumericDensePoly[Self.D], order: Int, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Truncate a polynomial to a specified formal precision.

        Args:
            poly: The polynomial supplying known coefficients.
            order: The nonnegative number of known coefficients.
            context: Explicit numerical arithmetic settings.
            budget: The shared work allowance.

        Returns:
            A formal series with zeros filled from the polynomial below the order.

        Raises:
            PolyError: DOMAIN for invalid inputs, LIMIT or CAPACITY on work limits.
        """
        if order < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "negative series order")
        _validate[Self.D](poly.coefficients, context, budget)
        var values = _values[_NumericOps[Self.D]](_view(poly.coefficients), min(order, len(poly.coefficients)), context, budget)
        return Self(NumericDensePoly[Self.D].from_coefficients(poly.variable, values^, budget), order)

    def _validate(self, context: Self.D.Context, mut budget: Budget) raises PolyError:
        _validate[Self.D](self.polynomial.coefficients, context, budget)
        _series_shape(self.order, len(self.polynomial.coefficients), budget)

    def _same(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError:
        self._validate(context, budget)
        other._validate(context, budget)
        self.polynomial._same(other.polynomial)
        _series_same_order(self.order, other.order)

    def truncate(self, order: Int, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Reduce precision without inventing unknown coefficients.

        Args:
            order: The requested lower precision.
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            A series at the requested order.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        _series_truncation(self.order, order)
        return Self.from_polynomial(self.polynomial, order, context, budget)

    def add(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Add two series at equal precision.

        Args:
            other: The second series in the same variable, domain and order.
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The truncated sum.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._same(other, context, budget)
        var values = _add_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), _view(other.polynomial.coefficients), self.order, False, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def sub(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Subtract two series at equal precision.

        Args:
            other: The second series in the same variable, domain and order.
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The truncated difference.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._same(other, context, budget)
        var values = _add_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), _view(other.polynomial.coefficients), self.order, True, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def mul(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Multiply modulo the common truncation order.

        Args:
            other: The second series in the same variable, domain and order.
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The truncated product.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._same(other, context, budget)
        var values = _mul_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), _view(other.polynomial.coefficients), self.order, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def divide(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Divide by a series with invertible constant coefficient.

        Args:
            other: The second series in the same variable, domain and order.
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The quotient modulo the common order.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._same(other, context, budget)
        var values = _divide_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), _view(other.polynomial.coefficients), self.order, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def inverse(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Invert a series with invertible constant coefficient.

        Args:
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The multiplicative inverse at the original precision.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        var values = _divide_values[_NumericOps[Self.D]](_view([_NumericOps[Self.D].integer(1, context, budget)]), _view(self.polynomial.coefficients), self.order, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def compose(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Compose with a series whose constant is exactly zero.

        Args:
            other: The second series in the same variable, domain and order.
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The composition at the common precision.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._same(other, context, budget)
        _series_composition[_NumericOps[Self.D]](_view(other.polynomial.coefficients), self.order, context, budget)
        var values = _compose_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), _view(other.polynomial.coefficients), self.order, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def derivative(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Differentiate and decrease known precision by one.

        Args:
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The formal derivative, with order max(0, order-1).

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        var values = _derivative_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, max(0, self.order - 1))

    def integrate(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Integrate with zero constant and increase precision by one.

        Args:
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The antiderivative when required coefficient divisions exist.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        var values = _integral_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, _sum(self.order, 1))

    def log(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Compute a formal logarithm with constant coefficient one.

        Args:
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The zero-constant logarithm at the original precision.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        var values = _log_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), self.order, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def exp(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Compute a formal exponential with constant coefficient zero.

        Args:
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The unit-constant exponential at the original precision.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        var values = _exp_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), self.order, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)

    def revert(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Invert composition for zero constant and invertible linear term.

        Args:
            context: Explicit rounding or enclosing arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The reversion at the original precision, requiring order at least two.

        Raises:
            PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        var values = _revert_values[_NumericOps[Self.D]](_view(self.polynomial.coefficients), self.order, context, budget)
        var poly = NumericDensePoly[Self.D].from_coefficients(self.polynomial.variable, values^, budget)
        return Self(poly^, self.order)
