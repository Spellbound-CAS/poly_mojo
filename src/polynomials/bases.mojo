"""Newton, Chebyshev and Bernstein polynomial representations."""
from .coefficients import CoefficientField
from .dense import DensePoly
from .interpolation import _validate
from ._toolkit_exact import _ExactOps
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind
from ._toolkit_ops import _values
from ._basis_ops import _basis_shape, _basis_values, _basis_evaluate, _newton_coefficients, _chebyshev_coefficients, _bernstein_coefficients, _bernstein_split


@fieldwise_init
struct BasisPolynomial[C: CoefficientField](Copyable, Movable):
    """Represent a polynomial in one of three structured bases.

    Parameters:
        C: The exact coefficient field.
    """
    comptime NEWTON = UInt8(0)
    """Products of successive x-minus-node factors."""
    comptime CHEBYSHEV = UInt8(1)
    """First-kind Chebyshev polynomials with T0=1 and T1=x."""
    comptime BERNSTEIN = UInt8(2)
    """Bernstein polynomials of the retained degree on [0,1]."""
    var variable: Int64
    """The polynomial variable."""
    var kind: UInt8
    """The basis tag."""
    var ring: Self.C.Context
    """The exact coefficient field context."""
    var coefficients: List[Self.C]
    """Ascending basis coefficients, retaining representation degree."""
    var nodes: List[Self.C]
    """The Newton nodes; empty for the other bases."""

    def _validate(self, mut budget: Budget) raises PolyError:
        _basis_shape(self.kind, len(self.coefficients), len(self.nodes), budget)
        _validate(self.coefficients, self.ring, budget)
        _validate(self.nodes, self.ring, budget)

    @staticmethod
    def newton(poly: DensePoly[Self.C], nodes: List[Self.C], mut budget: Budget) raises PolyError -> Self:
        """Convert monomial coefficients to a Newton basis.

        Args:
            poly: The monomial polynomial.
            nodes: Newton nodes, with count at least the polynomial degree; repetitions are allowed.
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing basis representation.

        Raises:
            PolyError: DOMAIN for invalid data or noninvertible conversion factors, LIMIT or CAPACITY on work limits.
        """
        poly.check_invariants(budget)
        _validate(nodes, poly.ring, budget)
        var coefficients = _newton_coefficients[_ExactOps[Self.C]](_view(poly.coefficients), _view(nodes), poly.ring, budget)
        var copied_nodes = _values[_ExactOps[Self.C]](_view(nodes), len(nodes), poly.ring, budget)
        return Self(poly.variable, Self.NEWTON, poly.ring.copy(), coefficients^, copied_nodes^)

    @staticmethod
    def chebyshev(poly: DensePoly[Self.C], mut budget: Budget) raises PolyError -> Self:
        """Convert monomial coefficients to a Chebyshev basis.

        Args:
            poly: The monomial polynomial.
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing basis representation.

        Raises:
            PolyError: DOMAIN for invalid data or noninvertible conversion factors, LIMIT or CAPACITY on work limits.
        """
        poly.check_invariants(budget)
        var coefficients = _chebyshev_coefficients[_ExactOps[Self.C]](_view(poly.coefficients), poly.ring, budget)
        var copied_nodes = List[Self.C]()
        return Self(poly.variable, Self.CHEBYSHEV, poly.ring.copy(), coefficients^, copied_nodes^)

    @staticmethod
    def bernstein(poly: DensePoly[Self.C], degree: Int, mut budget: Budget) raises PolyError -> Self:
        """Convert monomial coefficients to a retained Bernstein degree.

        Args:
            poly: The monomial polynomial.
            degree: The retained degree, at least the polynomial degree and zero.
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing basis representation.

        Raises:
            PolyError: DOMAIN for invalid data or noninvertible conversion factors, LIMIT or CAPACITY on work limits.
        """
        poly.check_invariants(budget)
        var coefficients = _bernstein_coefficients[_ExactOps[Self.C]](_view(poly.coefficients), degree, poly.ring, budget)
        var copied_nodes = List[Self.C]()
        return Self(poly.variable, Self.BERNSTEIN, poly.ring.copy(), coefficients^, copied_nodes^)

    def to_polynomial(self, mut budget: Budget) raises PolyError -> DensePoly[Self.C]:
        """Expand the basis into monomial coefficients.

        Args:
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing monomial polynomial.

        Raises:
            PolyError: DOMAIN for invalid representations, LIMIT or CAPACITY on work limits.
        """
        self._validate(budget)
        var values = _basis_values[_ExactOps[Self.C]](self.kind, _view(self.coefficients), _view(self.nodes), self.ring, budget)
        return DensePoly[Self.C].from_coefficients(self.variable, self.ring, values^, budget)

    def evaluate(self, point: Self.C, mut budget: Budget) raises PolyError -> Self.C:
        """Evaluate using nested Newton, Clenshaw or de Casteljau arithmetic.

        Args:
            point: The evaluation point or subdivision parameter.
            budget: The shared work allowance.

        Returns:
            The exact value, rounded value or uniform enclosure.

        Raises:
            PolyError: DOMAIN for invalid representations or incompatible points, LIMIT or CAPACITY on work limits.
        """
        self._validate(budget)
        if point.context() != self.ring:
            raise PolyError(PolyErrorKind.DOMAIN, "basis point in another field")
        return _basis_evaluate[_ExactOps[Self.C]](self.kind, _view(self.coefficients), _view(self.nodes), point, self.ring, budget)

    def subdivide(self, point: Self.C, mut budget: Budget) raises PolyError -> Tuple[Self, Self]:
        """Subdivide a Bernstein representation at an explicit parameter.

        Args:
            point: The evaluation point or subdivision parameter.
            budget: The shared work allowance.

        Returns:
            Representations of p(point*x) and p(point+(1-point)*x), retaining degree.

        Raises:
            PolyError: DOMAIN for invalid representations or incompatible points, LIMIT or CAPACITY on work limits.
        """
        self._validate(budget)
        if point.context() != self.ring:
            raise PolyError(PolyErrorKind.DOMAIN, "basis point in another field")
        if self.kind != Self.BERNSTEIN:
            raise PolyError(PolyErrorKind.DOMAIN, "subdivision requires a Bernstein representation")
        var split = _bernstein_split[_ExactOps[Self.C]](_view(self.coefficients), point, self.ring, budget)
        var left = List[Self.C]()
        var right = List[Self.C]()
        swap(left, split[0])
        swap(right, split[1])
        return (Self(self.variable, self.kind, self.ring.copy(), left^, List[Self.C]()), Self(self.variable, self.kind, self.ring.copy(), right^, List[Self.C]()))
