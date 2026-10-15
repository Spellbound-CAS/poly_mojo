"""Exact Gaussian rational coefficients with component-wise work accounting."""
from apn_mojo import ExactComplex
from .coefficients import CoefficientField, NoContext, QQ
from .errors import Budget, PolyError, PolyErrorKind


@fieldwise_init
struct QQi(CoefficientField):
    """An exact Gaussian rational coefficient a + b*i, with rational a and b."""

    var value: ExactComplex
    """The canonical apn_mojo value with two exact rational components."""

    comptime Context = NoContext
    """The unique field context; no precision or rounding mode is needed."""

    @staticmethod
    def zero(ctx: NoContext) -> Self:
        """Construct zero.

        Args:
            ctx: The field context.

        Returns:
            The additive identity.
        """
        return Self(ExactComplex(0))

    @staticmethod
    def one(ctx: NoContext) -> Self:
        """Construct one.

        Args:
            ctx: The field context.

        Returns:
            The multiplicative identity.
        """
        return Self(ExactComplex(1))

    def context(self) -> NoContext:
        """Report the unique field context.

        Returns:
            The field's context.
        """
        return NoContext()

    def is_zero(self) -> Bool:
        """Test whether both exact components are zero.

        Returns:
            Whether the value is zero.
        """
        return self.value.is_real() and self.value.real().sign() == 0

    def size(self) -> Int:
        """Report the sum of the rational component sizes in bytes.

        Returns:
            The positive size used by the work model.
        """
        return QQ(self.value.real()).size() + QQ(self.value.imag()).size()

    def __eq__(self, other: Self) -> Bool:
        """Compare exact values under the fixed-signature protocol exception.

        Args:
            other: The value to compare.

        Returns:
            Whether both canonical rational components agree.
        """
        return self.value == other.value

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Add the exact coefficients with budgeted rational arithmetic.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact sum.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        var real = QQ(self.value.real()).add(QQ(other.value.real()), budget)
        var imag = QQ(self.value.imag()).add(QQ(other.value.imag()), budget)
        return Self(ExactComplex(real.value, imag.value))

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Subtract the exact coefficients with budgeted rational arithmetic.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact difference.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        var real = QQ(self.value.real()).sub(QQ(other.value.real()), budget)
        var imag = QQ(self.value.imag()).sub(QQ(other.value.imag()), budget)
        return Self(ExactComplex(real.value, imag.value))

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Multiply exactly, charging each rational intermediate operation.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact product.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        var a = QQ(self.value.real())
        var b = QQ(self.value.imag())
        var c = QQ(other.value.real())
        var d = QQ(other.value.imag())
        var real = a.mul(c, budget).sub(b.mul(d, budget), budget)
        var imag = a.mul(d, budget).add(b.mul(c, budget), budget)
        return Self(ExactComplex(real.value, imag.value))

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Negate both rational components after charging their sizes.

        Args:
            budget: The shared work allowance.

        Returns:
            The exact additive inverse.

        Raises:
            PolyError: LIMIT if the allowance is insufficient.
        """
        var real = QQ(self.value.real()).neg(budget)
        var imag = QQ(self.value.imag()).neg(budget)
        return Self(ExactComplex(real.value, imag.value))

    def inverse(self, mut budget: Budget) raises PolyError -> Self:
        """Invert a nonzero coefficient using its conjugate and squared norm.

        Args:
            budget: The shared work allowance.

        Returns:
            The exact reciprocal (a - b*i)/(a*a + b*b).

        Raises:
            PolyError: DOMAIN on zero, LIMIT or CAPACITY on work limits.
        """
        if self.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "inverse of zero")
        var a = QQ(self.value.real())
        var b = QQ(self.value.imag())
        var reciprocal_norm = a.mul(a, budget).add(b.mul(b, budget), budget).inverse(budget)
        var real = a.mul(reciprocal_norm, budget)
        var imag = b.neg(budget).mul(reciprocal_norm, budget)
        return Self(ExactComplex(real.value, imag.value))
