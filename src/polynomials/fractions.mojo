"""Canonical exact rational functions over QQ."""
from .errors import Budget, PolyError, PolyErrorKind
from .coefficients import QQ
from .sparse import SparsePoly
from .certificates import CancellationCertificate
from .gcd import gcd
from ._views import _zero


@fieldwise_init
struct RationalFunction(Copyable, Movable, Equatable, Writable):
    """A reduced polynomial fraction with a monic nonzero denominator."""
    var numerator: SparsePoly[QQ]
    """The canonical numerator; zero has denominator one."""
    var denominator: SparsePoly[QQ]
    """The monic denominator, coprime to the numerator."""

    @staticmethod
    def _reduce(p: SparsePoly[QQ], q: SparsePoly[QQ], mut common: SparsePoly[QQ], mut budget: Budget) raises PolyError -> Self:
        p._same(q, budget)
        if q.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "zero rational-function denominator")
        var result = gcd(p, q, budget)
        ref scalar = result.cofactors[1].terms[len(result.cofactors[1].terms) - 1].coefficient
        var inverse = scalar.inverse(budget)
        common = result.gcd.scale(scalar, budget)
        return Self(result.cofactors[0].scale(inverse, budget), result.cofactors[1].scale(inverse, budget))

    @staticmethod
    def normalize(p: SparsePoly[QQ], q: SparsePoly[QQ], mut budget: Budget) raises PolyError -> Tuple[Self, CancellationCertificate]:
        """Reduce a polynomial fraction and produce its cancellation identity.

        Args:
            p: The numerator.
            q: The nonzero denominator in the same ring and order.
            budget: The shared work allowance.

        Returns:
            The unique normalized fraction and its cancellation certificate.

        Raises:
            PolyError: DOMAIN for a zero denominator or incompatible inputs, LIMIT or CAPACITY on work limits.
        """
        var common = _zero(p, budget)
        var fraction = Self._reduce(p, q, common, budget)
        var certificate = CancellationCertificate(common^, (fraction.numerator._copy(budget), fraction.denominator._copy(budget)))
        return (fraction^, certificate^)

    def __eq__(self, other: Self) -> Bool:
        """Compare canonical fractions under the fixed-signature protocol exception.

        Args:
            other: Another normalized fraction.

        Returns:
            Whether both canonical polynomials and their contexts agree.
        """
        return self.numerator == other.numerator and self.denominator == other.denominator

    def _sum(self, other: Self, subtract: Bool, mut budget: Budget) raises PolyError -> Self:
        self.numerator._same(other.numerator, budget)
        var right = other.numerator.mul(self.denominator, budget)
        var numerator = self.numerator.mul(other.denominator, budget)
        numerator = numerator.sub(right, budget) if subtract else numerator.add(right, budget)
        var denominator = self.denominator.mul(other.denominator, budget)
        var common = _zero(numerator, budget)
        return Self._reduce(numerator, denominator, common, budget)

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute a normalized sum.

        Args:
            other: The second fraction in the same ring and order.
            budget: The shared work allowance.

        Returns:
            The canonical sum.

        Raises:
            PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
        """
        return self._sum(other, False, budget)

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute a normalized difference.

        Args:
            other: The subtracted fraction in the same ring and order.
            budget: The shared work allowance.

        Returns:
            The canonical difference.

        Raises:
            PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
        """
        return self._sum(other, True, budget)

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute a normalized product.

        Args:
            other: The second fraction in the same ring and order.
            budget: The shared work allowance.

        Returns:
            The canonical product.

        Raises:
            PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
        """
        var numerator = self.numerator.mul(other.numerator, budget)
        var denominator = self.denominator.mul(other.denominator, budget)
        var common = _zero(numerator, budget)
        return Self._reduce(numerator, denominator, common, budget)

    def div(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute a normalized quotient by a nonzero fraction.

        Args:
            other: The divisor in the same ring and order.
            budget: The shared work allowance.

        Returns:
            The canonical quotient.

        Raises:
            PolyError: DOMAIN for a zero divisor or incompatible inputs, LIMIT or CAPACITY on work limits.
        """
        if other.numerator.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "division by zero rational function")
        var numerator = self.numerator.mul(other.denominator, budget)
        var denominator = self.denominator.mul(other.numerator, budget)
        var common = _zero(numerator, budget)
        return Self._reduce(numerator, denominator, common, budget)

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Negate the fraction while retaining its canonical denominator.

        Args:
            budget: The shared work allowance.

        Returns:
            The canonical additive inverse.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        return Self(self.numerator.neg(budget), self.denominator._copy(budget))
