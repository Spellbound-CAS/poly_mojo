"""Budgeted conversions between the public dense and sparse representations."""
from std.iter import Iterable
from polynomials import DensePoly, Term, Budget, PolyError, PolyErrorKind
from . import dense, sparse
from ._common import _Element, _Kind, _Ring, _allocation, _product


def to_sparse[I: Iterable](values: I, mut budget: Budget) raises PolyError -> sparse.Array[Int, _Kind[_Element[I]]]:
    """Convert ascending coefficients to nonzero (power, coefficient) terms.

    The coefficient domain and modulus are retained. Dense powers that exceed
    UInt32.MAX raise CAPACITY. Input preparation and conversion share the budget.

    Parameters:
        I: The inferred input iterable type.

    Args:
        values: The coefficient sequence, sparse terms, or prepared array.
        budget: The caller-owned work allowance, shared through preparation and computation.

    Returns:
        The univariate nonzero terms, preserving domain and modulus.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    comptime k = _Kind[_Element[I]]
    var p = dense.asarray(values, budget)
    var order = sparse._order[Int](budget)
    # The core conversion grows its output geometrically. Account for term
    # storage and one exponent pair per possible nonconstant term in advance.
    _allocation[Term[_Ring[k]]](_product(len(p), 4), budget)
    _allocation[Tuple[Int64, UInt32]](len(p), budget)
    return sparse.Array[Int, k](p._core[].to_sparse(order, budget), budget)


def to_sparse[I: Iterable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> sparse.Array[Int, _Kind[_Element[I]]]:
    """Convert dense coefficients to sparse terms with an optional work limit.

    Parameters:
        I: The inferred input iterable type.

    Args:
        values: The coefficient sequence, sparse terms, or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The univariate nonzero terms, preserving domain and modulus.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    var budget = Budget(max_steps)
    return to_sparse(values, budget)


def to_dense[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, mut budget: Budget) raises PolyError -> dense.Array[_Kind[C]] where _Element[I] == Tuple[P, C]:
    """Fill missing powers of a univariate sparse polynomial with exact zeros.

    Integer powers and one-element power tuples are accepted. The required
    dense length is charged before allocation; use a finite budget for large
    sparse powers. Multivariate power shapes raise DOMAIN even when constant.

    Parameters:
        I: The inferred input iterable type.
        P: The power shape: Int or a tuple of Int powers.
        C: The inferred input coefficient type.

    Args:
        values: The coefficient sequence, sparse terms, or prepared array.
        budget: The caller-owned work allowance, shared through preparation and computation.

    Returns:
        The ascending coefficients, including zeros for missing powers.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    if sparse._width[P]() != 1:
        raise PolyError(PolyErrorKind.DOMAIN, "dense conversion needs a one-variable power shape")
    comptime k = _Kind[C]
    var p = sparse.asarray[I, P, C](values, budget)
    return dense.Array[k](DensePoly[_Ring[k]].from_sparse(p._core[], budget, variable=Optional(Int64(0))), budget)


def to_dense[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> dense.Array[_Kind[C]] where _Element[I] == Tuple[P, C]:
    """Convert univariate sparse terms to dense coefficients with a work limit.

    Parameters:
        I: The inferred input iterable type.
        P: The power shape: Int or a tuple of Int powers.
        C: The inferred input coefficient type.

    Args:
        values: The coefficient sequence, sparse terms, or prepared array.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The ascending coefficients, including zeros for missing powers.

    Raises:
        PolyError: LIMIT on exhausted work, CAPACITY on unrepresentable sizes, or DOMAIN for incompatible input.
    """
    var budget = Budget(max_steps)
    return to_dense[I, P, C](values, budget)


from polynomials import NumericDomain, NumericSparsePoly, MonomialOrder
from ._interop import _variables


def to_sparse[D: NumericDomain](values: dense.NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> sparse.NumericArray[Int, D]:
    """Convert numerical dense storage to sparse terms with explicit zero arithmetic.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        budget: The shared work allowance, including preparation, computation, and result storage.
        context: Explicit output precision, rounding settings, or ball working precision.

    Returns:
        The result described above.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var order = MonomialOrder(MonomialOrder.LEX, _variables(1, List[Int64](), budget), budget)
    return sparse.NumericArray[Int, D](NumericSparsePoly[D].from_dense(values._core[], order, context, budget), budget)


def to_sparse[D: NumericDomain](values: dense.NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> sparse.NumericArray[Int, D]:
    """Convert numerical dense coefficients to sparse terms.

    Parameters:
        D: The APN numerical domain associated with the stored scalars.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        context: Explicit output precision, rounding settings, or ball working precision.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The result described above.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return to_sparse[D](values, budget, context=context)


def to_dense[P: Copyable & Deinitable, D: NumericDomain](values: sparse.NumericArray[P, D], mut budget: Budget, *, context: D.Context) raises PolyError -> dense.NumericArray[D]:
    """Materialize a univariate numerical polynomial, charging gap zeros before allocation.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.
        D: The APN numerical domain associated with the stored scalars.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        budget: The shared work allowance, including preparation, computation, and result storage.
        context: Explicit output precision, rounding settings, or ball working precision.

    Returns:
        The result described above.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert sparse._width[P]() == 1, "dense conversion requires a univariate power shape"
    return dense.NumericArray[D](values._core[].to_dense(0, context, budget), budget)


def to_dense[P: Copyable & Deinitable, D: NumericDomain](values: sparse.NumericArray[P, D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> dense.NumericArray[D]:
    """Convert numerical sparse terms to ascending coefficients.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.
        D: The APN numerical domain associated with the stored scalars.

    Args:
        values: Input coefficients, sparse terms, or prepared values as required by this overload.
        context: Explicit output precision, rounding settings, or ball working precision.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The result described above.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return to_dense[P, D](values, budget, context=context)
