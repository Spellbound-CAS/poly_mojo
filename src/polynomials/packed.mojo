"""Sparse polynomials in packed-term storage, with products, powers and exact quotients that stay packed."""
from .coefficients import CoefficientRing, ZZ, ModInt
from .errors import Budget, PolyError, PolyErrorKind, _charge_size
from .sparse import SparsePoly
from .division import divide_exact
from ._packed_poly import _PackedPoly, _KEY_LIMIT
from ._packed_products import _packed_poly_product, _packed_poly_power, _packed_poly_divide


struct PackedPoly[C: CoefficientRing](Copyable, Movable, Equatable, Sized):
    """A sparse polynomial stored as one exponent word per term and a contiguous coefficient list.

    The order is LEX over the variables 0..n-1 in that order. Each monomial is
    one word by mixed-radix substitution with radix_i = (largest exponent of
    variable i) + 1, variable 0 most significant, so word order is LEX order;
    the product of the radices is below 2**62. Products, powers and exact
    quotients emit their words directly, with no per-term monomial lists until
    `to_sparse`.

    Parameters:
        C: The coefficient ring type.
    """
    var _packed: _PackedPoly[Self.C]

    def __init__(out self, *, var _packed: _PackedPoly[Self.C]):
        """Wrap packed storage whose radix is the exact exponent bound.

        Args:
            _packed: Canonical packed storage.
        """
        self._packed = _packed^

    @staticmethod
    def from_sparse(p: SparsePoly[Self.C], mut budget: Budget) raises PolyError -> Optional[Self]:
        """Pack a sparse polynomial, or return no value when it does not fit.

        Args:
            p: A polynomial in LEX order over the variables 0..n-1 in that order.
            budget: The shared work allowance.

        Returns:
            The packed polynomial, or None for another order or a key range of 2**62 or more.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        var packed = _PackedPoly[Self.C].from_sparse(p, budget)
        if not packed:
            return None
        return Self(_packed=packed.take())

    def to_sparse(self, mut budget: Budget) raises PolyError -> SparsePoly[Self.C]:
        """Materialize the canonical sparse form.

        Args:
            budget: The shared work allowance.

        Returns:
            The same polynomial as a `SparsePoly`.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        return self._packed.to_sparse(budget)

    def is_zero(self) -> Bool:
        """Report whether the polynomial has no terms.

        Returns:
            Whether the polynomial is zero.
        """
        return not len(self._packed)

    def __len__(self) -> Int:
        """Count the nonzero terms.

        Returns:
            The number of terms.
        """
        return len(self._packed)

    def __eq__(self, other: Self) -> Bool:
        """Compare rings, orders and terms.

        Args:
            other: The other polynomial.

        Returns:
            Whether both are the same polynomial in the same ring and order.
        """
        return self._packed.ring == other._packed.ring and self._packed == other._packed

    def _same(self, other: Self, mut budget: Budget) raises PolyError:
        _charge_size(len(self._packed.order.variables), budget)
        if self._packed.order != other._packed.order or self._packed.ring != other._packed.ring:
            raise PolyError(PolyErrorKind.DOMAIN, "different polynomial rings or orders")

    def _top(self) -> List[UInt64]:
        var top = List[UInt64](capacity=len(self._packed.radix))
        for r in self._packed.radix:
            top.append(r - 1)
        return top^


def _fits(top: List[UInt128]) -> Bool:
    """Whether radices top_i + 1 keep the key range below 2**62."""
    var total = UInt128(1)
    for e in top:
        if e >= UInt128(_KEY_LIMIT):
            return False
        total *= e + 1
        if total >= UInt128(_KEY_LIMIT):
            return False
    return True


def packed_mul[C: CoefficientRing](a: PackedPoly[C], b: PackedPoly[C], mut budget: Budget) raises PolyError -> Optional[PackedPoly[C]]:
    """Multiply packed polynomials, keeping the product packed.

    ZZ coefficients of at most 63 bits and prime-field residues take the packed
    product kernel; other operands multiply as `SparsePoly.mul` does and are
    packed again.

    Parameters:
        C: The coefficient ring type.

    Args:
        a: The left factor.
        b: The right factor in the same ring and order.
        budget: The shared work allowance.

    Returns:
        The product, or None when its key range would reach 2**62.

    Raises:
        PolyError: DOMAIN for incompatible operands, LIMIT or CAPACITY on work limits.
    """
    a._same(b, budget)
    var ta = a._top()
    var tb = b._top()
    var top = List[UInt128](capacity=len(ta))
    for i in range(len(ta)):
        top.append(UInt128(ta[i]) + UInt128(tb[i]))
    if not _fits(top):
        return None
    comptime if C == ZZ or C == ModInt:
        if len(top):
            var direct = _packed_poly_product(a._packed, b._packed, budget)
            if direct:
                return PackedPoly[C](_packed=direct.take())
    var product = a.to_sparse(budget).mul(b.to_sparse(budget), budget)
    return PackedPoly[C].from_sparse(product, budget)


def packed_pow[C: CoefficientRing](f: PackedPoly[C], n: UInt64, mut budget: Budget) raises PolyError -> Optional[PackedPoly[C]]:
    """Raise a packed polynomial to a natural power, keeping the power packed.

    Few-term ZZ bases take the power recurrence where `SparsePoly.pow` would;
    otherwise the power squares from the exponent's leading bit down with
    packed products. Zero power equals one.

    Parameters:
        C: The coefficient ring type.

    Args:
        f: The base.
        n: The exponent.
        budget: The shared work allowance.

    Returns:
        The power, or None when its key range would reach 2**62.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    var tf = f._top()
    var top = List[UInt128](capacity=len(tf))
    for e in tf:
        top.append(UInt128(e) * UInt128(n))
    if not _fits(top):
        return None
    comptime if C == ZZ or C == ModInt:
        if len(top) and n >= 1 and n <= UInt64(1 << 62):
            var direct = _packed_poly_power(f._packed, Int(n), budget)
            if direct:
                return PackedPoly[C](_packed=direct.take())
    var power = f.to_sparse(budget).pow(n, budget)
    return PackedPoly[C].from_sparse(power, budget)


def packed_divide_exact[C: CoefficientRing](f: PackedPoly[C], g: PackedPoly[C], mut budget: Budget) raises PolyError -> Optional[PackedPoly[C]]:
    """Return the exact packed quotient f / g, or no value when g does not divide f.

    ZZ divisors with coefficients of at most 63 bits and prime-field divisors
    take chunked packed division; other operands divide as `divide_exact`
    does. The quotient is packed in its own radix.

    Parameters:
        C: ZZ or a coefficient field.

    Args:
        f: The dividend.
        g: The nonzero divisor in the same ring and order.
        budget: The shared work allowance.

    Returns:
        The quotient if g divides f, otherwise None.

    Raises:
        PolyError: DOMAIN for invalid operands, UNSUPPORTED for another ring, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "division by zero polynomial")
    comptime if C == ZZ or C == ModInt:
        var direct = _packed_poly_divide(f._packed, g._packed, budget)
        if direct:
            var quotient = direct.take()
            if not quotient:
                return None
            return PackedPoly[C](_packed=quotient.take())
    var sparse = divide_exact(f.to_sparse(budget), g.to_sparse(budget), budget)
    if not sparse:
        return None
    return PackedPoly[C].from_sparse(sparse.value(), budget)
