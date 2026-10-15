"""Newton, Chebyshev and Bernstein polynomial representations."""
from .numeric_domains import NumericDomain
from .numeric_dense import NumericDensePoly
from ._toolkit_numeric import _NumericOps, _validate
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind
from ._toolkit_ops import _values
from ._basis_ops import _basis_shape, _basis_values, _basis_evaluate, _newton_coefficients, _chebyshev_coefficients, _bernstein_coefficients, _bernstein_split


@fieldwise_init
struct NumericBasisPolynomial[D: NumericDomain](Copyable, Movable):
    """Represent a polynomial in one of three structured bases.

    Parameters:
        D: The numerical domain.
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
    var coefficients: List[Self.D.Scalar]
    """Ascending basis coefficients, retaining representation degree."""
    var nodes: List[Self.D.Scalar]
    """The Newton nodes; empty for the other bases."""

    def _validate(self, context: Self.D.Context, mut budget: Budget) raises PolyError:
        _basis_shape(self.kind, len(self.coefficients), len(self.nodes), budget)
        _validate[Self.D](self.coefficients, context, budget)
        _validate[Self.D](self.nodes, context, budget)

    @staticmethod
    def newton(poly: NumericDensePoly[Self.D], nodes: List[Self.D.Scalar], context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Convert monomial coefficients to a Newton basis.

        Args:
            poly: The monomial polynomial.
            nodes: Newton nodes, with count at least the polynomial degree; repetitions are allowed.
            context: Explicit numerical arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing basis representation.

        Raises:
            PolyError: DOMAIN for invalid data or noninvertible conversion factors, LIMIT or CAPACITY on work limits.
        """
        _validate[Self.D](poly.coefficients, context, budget)
        _validate[Self.D](nodes, context, budget)
        var coefficients = _newton_coefficients[_NumericOps[Self.D]](_view(poly.coefficients), _view(nodes), context, budget)
        var copied_nodes = _values[_NumericOps[Self.D]](_view(nodes), len(nodes), context, budget)
        return Self(poly.variable, Self.NEWTON, coefficients^, copied_nodes^)

    @staticmethod
    def chebyshev(poly: NumericDensePoly[Self.D], context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Convert monomial coefficients to a Chebyshev basis.

        Args:
            poly: The monomial polynomial.
            context: Explicit numerical arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing basis representation.

        Raises:
            PolyError: DOMAIN for invalid data or noninvertible conversion factors, LIMIT or CAPACITY on work limits.
        """
        _validate[Self.D](poly.coefficients, context, budget)
        var coefficients = _chebyshev_coefficients[_NumericOps[Self.D]](_view(poly.coefficients), context, budget)
        var copied_nodes = List[Self.D.Scalar]()
        return Self(poly.variable, Self.CHEBYSHEV, coefficients^, copied_nodes^)

    @staticmethod
    def bernstein(poly: NumericDensePoly[Self.D], degree: Int, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Convert monomial coefficients to a retained Bernstein degree.

        Args:
            poly: The monomial polynomial.
            degree: The retained degree, at least the polynomial degree and zero.
            context: Explicit numerical arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing basis representation.

        Raises:
            PolyError: DOMAIN for invalid data or noninvertible conversion factors, LIMIT or CAPACITY on work limits.
        """
        _validate[Self.D](poly.coefficients, context, budget)
        var coefficients = _bernstein_coefficients[_NumericOps[Self.D]](_view(poly.coefficients), degree, context, budget)
        var copied_nodes = List[Self.D.Scalar]()
        return Self(poly.variable, Self.BERNSTEIN, coefficients^, copied_nodes^)

    def to_polynomial(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[Self.D]:
        """Expand the basis into monomial coefficients.

        Args:
            context: Explicit numerical arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The exact, rounded or enclosing monomial polynomial.

        Raises:
            PolyError: DOMAIN for invalid representations, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        var values = _basis_values[_NumericOps[Self.D]](self.kind, _view(self.coefficients), _view(self.nodes), context, budget)
        return NumericDensePoly[Self.D].from_coefficients(self.variable, values^, budget)

    def evaluate(self, point: Self.D.Scalar, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Evaluate using nested Newton, Clenshaw or de Casteljau arithmetic.

        Args:
            point: The evaluation point or subdivision parameter.
            context: Explicit numerical arithmetic settings.
            budget: The shared work allowance.

        Returns:
            The exact value, rounded value or uniform enclosure.

        Raises:
            PolyError: DOMAIN for invalid representations or incompatible points, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        return _basis_evaluate[_NumericOps[Self.D]](self.kind, _view(self.coefficients), _view(self.nodes), point, context, budget)

    def subdivide(self, point: Self.D.Scalar, context: Self.D.Context, mut budget: Budget) raises PolyError -> Tuple[Self, Self]:
        """Subdivide a Bernstein representation at an explicit parameter.

        Args:
            point: The evaluation point or subdivision parameter.
            context: Explicit numerical arithmetic settings.
            budget: The shared work allowance.

        Returns:
            Representations of p(point*x) and p(point+(1-point)*x), retaining degree.

        Raises:
            PolyError: DOMAIN for invalid representations or incompatible points, LIMIT or CAPACITY on work limits.
        """
        self._validate(context, budget)
        if self.kind != Self.BERNSTEIN:
            raise PolyError(PolyErrorKind.DOMAIN, "subdivision requires a Bernstein representation")
        var split = _bernstein_split[_NumericOps[Self.D]](_view(self.coefficients), point, context, budget)
        var left = List[Self.D.Scalar]()
        var right = List[Self.D.Scalar]()
        swap(left, split[0])
        swap(right, split[1])
        return (Self(self.variable, self.kind, left^, List[Self.D.Scalar]()), Self(self.variable, self.kind, right^, List[Self.D.Scalar]()))
