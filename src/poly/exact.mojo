"""Functional access and independent checks for retained exact-algebra results."""
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink, BallDomain, QQ
from . import _determinant
from ._determinant import Certificate as DeterminantCertificate
comptime _Determinant[k: Int] = ProjectionResultantCertificate if k == 1 else DeterminantCertificate[_Ring[k]]
from std.memory import ArcPointer
from std.builtin.rebind import downcast
from polynomials import Budget, PolyError, PolyErrorKind, CoefficientRing, CoefficientField, DensePoly, SparsePoly, GcdResult, CoprimalityCertificate, ResultantResult, DivisionResult, ListDivision, Factorization, MultivariateFactorization, SquarefreeDecomposition
from polynomials import check_common_divisor, check_bezout, check_coprimality, check_resultant_identity, check_division, check_list_division, check_factorization, check_multivariate_factorization, check_squarefree_product, check_interpolation, check_hermite_interpolation
from polynomials import SparseSquarefreeDecomposition, check_sparse_squarefree_product
from polynomials import ProjectionResultantCertificate, check_projection_resultant
from . import dense, sparse
from ._common import _Ring, _Scalar, _allocation, _charge_size, _unwrap
from ._interop import _dense_sparse, _dense_copy

comptime _Array[k: Int, P: Copyable & Deinitable, is_sparse: Bool] = sparse.Array[P, k] if is_sparse else dense.Array[k]


@fieldwise_init
struct _Data[C: CoefficientRing, E: Copyable & Deinitable](Copyable, Movable):
    var input: SparsePoly[Self.C]
    var others: List[SparsePoly[Self.C]]
    var proof: Self.E
    var axis: Int
    var algebraic_images: List[Optional[Factorization[Self.C]]]


@fieldwise_init
struct Result[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable]( _Archive, Copyable, Movable):
    """An exact result retaining its original operands and immutable evidence.

    Use the functional value, quotient, remainder, factors and check helpers.
    The checker’s claim depends on the operation that produced the result.
    Multivariate algebraic factorization also retains image certificates internally;
    keep the complete Result for its irreducibility check.

    Parameters:
        k: The inferred exact coefficient kind.
        P: The original sparse power shape, or Int for dense arrays.
        is_sparse: Whether output polynomials use sparse arrays.
        E: The retained core evidence type.
    """
    var _data: ArcPointer[_Data[_Ring[Self.k], Self.E]]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] Self.E:
        """Borrow the evidence for an independent core checker.

        This is the core product certificate for sparse algebraic factorization.
        Auxiliary irreducibility certificates remain in the complete Result
        consumed by exact.check.

        Returns:
            Immutable evidence; copy it explicitly for mutable interoperability.
        """
        return self._data[].proof



    @staticmethod
    def _archive_shape() -> Int:
        from .sparse import _width
        return _width[Self.P]()

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Self.E
        _result_io._depth(depth, budget)
        _result_io._write_node[_Data[_Ring[Self.k], Self.E], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Self.E
        _result_io._depth(depth, budget)
        var field654 = _result_io._read_node[_Data[_Ring[Self.k], Self.E], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_Data[_Ring[Self.k], Self.E]](1, budget)
        var field653 = ArcPointer(field654^)
        return Self(field653^)

@fieldwise_init
struct _Interpolation[C: CoefficientRing](Copyable, Movable):
    var points: List[Self.C]
    var values: List[Self.C]
    var jets: List[List[Self.C]]
    var derivatives: Bool
    var hermite: Bool


def _retain[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](var input: SparsePoly[_Ring[k]], var others: List[SparsePoly[_Ring[k]]], var proof: E, mut budget: Budget, *, axis: Int = 0) raises PolyError -> Result[k, P, is_sparse, E]:
    _allocation[_Data[_Ring[k], E]](1, budget)
    return Result[k, P, is_sparse, E](ArcPointer(_Data[_Ring[k], E](input^, others^, proof^, axis, List[Optional[Factorization[_Ring[k]]]]())))


def _output[k: Int, P: Copyable & Deinitable, is_sparse: Bool](value: SparsePoly[_Ring[k]], mut budget: Budget) raises PolyError -> _Array[k, P, is_sparse]:
    comptime if is_sparse:
        from ._interop import _sparse_copy, _variables
        return rebind_var[_Array[k, P, is_sparse]](sparse.Array[P, k](_sparse_copy(value, _variables(sparse._width[P](), List[Int64](), budget), False, budget), budget))
    else:
        return rebind_var[_Array[k, P, is_sparse]](dense._from_sparse[k](value, budget))


def _dense_output[k: Int, P: Copyable & Deinitable, is_sparse: Bool](value: DensePoly[_Ring[k]], mut budget: Budget) raises PolyError -> _Array[k, P, is_sparse]:
    comptime if is_sparse:
        var p = _dense_sparse(value, budget)
        return _output[k, P, is_sparse](p, budget)
    else:
        return rebind_var[_Array[k, P, is_sparse]](dense.Array[k](_dense_copy(value, 0, budget), budget))


def check[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> Bool:
    """Verify the exact relation claimed by a retained algebra result.

    GCD checks establish maximality up to a unit. Resultant checks establish an ideal-membership identity or the full determinant when produced by resultant_certified. Factorization checks over AA and QQbar also establish normalized irreducible factors with exact multiplicities. Other factorization checks establish product identities only; squarefree checks do not establish squarefreeness.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Whether the retained evidence establishes the documented relation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    ref data = result._data[]
    comptime if E == GcdResult[C]:
        if len(data.others) != 1:
            return False
        ref proof = rebind[GcdResult[C]](data.proof)
        if not check_common_divisor(data.input, data.others[0], proof, budget):
            return False
        if proof.gcd.is_zero():
            return data.input.is_zero() and data.others[0].is_zero()
        if proof.bezout:
            return check_bezout(data.input, data.others[0], proof, budget)
        if proof.coprimality:
            return check_coprimality(proof.cofactors[0], proof.cofactors[1], proof.coprimality.value(), budget)
        return False
    elif E == CoprimalityCertificate[C]:
        return len(data.others) == 1 and check_coprimality(data.input, data.others[0], rebind[CoprimalityCertificate[C]](data.proof), budget)
    elif E == ResultantResult[C]:
        return len(data.others) == 1 and check_resultant_identity(data.input, data.others[0], rebind[ResultantResult[C]](data.proof), budget)
    elif E == DeterminantCertificate[C]:
        return len(data.others) == 1 and Int64(data.axis) == rebind[DeterminantCertificate[C]](data.proof).variable and _determinant._check(data.input, data.others[0], rebind[DeterminantCertificate[C]](data.proof), budget)
    elif E == ProjectionResultantCertificate:
        comptime assert k == 1, "determinant resultants use rational coefficients"
        return len(data.others) == 1 and Int64(data.axis) == rebind[ProjectionResultantCertificate](data.proof).variable and check_projection_resultant(rebind[SparsePoly[QQ]](data.input), rebind[SparsePoly[QQ]](data.others[0]), rebind[ProjectionResultantCertificate](data.proof), budget)
    elif E == DivisionResult[C]:
        if len(data.others) != 1:
            return False
        var a = DensePoly[C].from_sparse(data.input, budget, variable=Int64(0))
        var b = DensePoly[C].from_sparse(data.others[0], budget, variable=Int64(0))
        return check_division(a, b, rebind[DivisionResult[C]](data.proof), budget)
    elif E == ListDivision[C]:
        return check_list_division(data.input, data.others, rebind[ListDivision[C]](data.proof), budget)
    elif E == Factorization[C]:
        var p = DensePoly[C].from_sparse(data.input, budget, variable=Int64(0))
        comptime if k == 4 or k == 5:
            from . import _algebraic_factor
            return _algebraic_factor._check(p, rebind[Factorization[C]](data.proof), budget)
        else:
            return check_factorization(p, rebind[Factorization[C]](data.proof), budget)
    elif E == MultivariateFactorization[C]:
        comptime if k == 4 or k == 5:
            from . import _algebraic_factor
            return _algebraic_factor._check_sparse_evidence(data.input, rebind[MultivariateFactorization[C]](data.proof), data.algebraic_images, budget)
        else:
            return check_multivariate_factorization(data.input, rebind[MultivariateFactorization[C]](data.proof), budget)
    elif E == SparseSquarefreeDecomposition[C]:
        return check_sparse_squarefree_product(data.input, rebind[SparseSquarefreeDecomposition[C]](data.proof), budget)
    elif E == SquarefreeDecomposition[C]:
        var p = DensePoly[C].from_sparse(data.input, budget, variable=Int64(0))
        return check_squarefree_product(p, rebind[SquarefreeDecomposition[C]](data.proof), budget)
    else:
        comptime assert E == _Interpolation[C], "unsupported exact result evidence"
        comptime F = downcast[C, CoefficientField]
        ref proof = rebind[_Interpolation[C]](data.proof)
        var p = DensePoly[C].from_sparse(data.input, budget, variable=Int64(0))
        if proof.hermite:
            return check_hermite_interpolation(rebind[DensePoly[F]](p), rebind[List[F]](proof.points), rebind[List[List[F]]](proof.jets), budget, derivatives=proof.derivatives)
        return check_interpolation(rebind[DensePoly[F]](p), rebind[List[F]](proof.points), rebind[List[F]](proof.values), budget)


def check[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the exact relation claimed by a retained algebra result.

    GCD checks establish maximality up to a unit. Resultant checks establish an ideal-membership identity or the full determinant when produced by resultant_certified. Factorization checks over AA and QQbar also establish normalized irreducible factors with exact multiplicities. Other factorization checks establish product identities only; squarefree checks do not establish squarefreeness.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Whether the retained evidence establishes the documented relation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def value[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> _Array[k, P, is_sparse]:
    """Copy the polynomial value of a GCD, resultant, or interpolation result.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        An independent array in the original dense or sparse representation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    comptime if E == GcdResult[C]:
        return _output[k, P, is_sparse](rebind[GcdResult[C]](result._data[].proof).gcd, budget)
    elif E == ResultantResult[C]:
        return _output[k, P, is_sparse](rebind[ResultantResult[C]](result._data[].proof).resultant, budget)
    elif E == DeterminantCertificate[C]:
        return _output[k, P, is_sparse](rebind[DeterminantCertificate[C]](result._data[].proof).resultant, budget)
    elif E == ProjectionResultantCertificate:
        comptime assert k == 1, "determinant resultants use rational coefficients"
        return _output[k, P, is_sparse](rebind[SparsePoly[C]](rebind[ProjectionResultantCertificate](result._data[].proof).resultant), budget)
    else:
        comptime assert E == _Interpolation[C], "use quotient, remainder or factors for this result"
        return _output[k, P, is_sparse](result._data[].input, budget)


def value[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Array[k, P, is_sparse]:
    """Copy the polynomial value of a GCD, resultant, or interpolation result.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        An independent array in the original dense or sparse representation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return value(result, budget)


def quotient[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> _Array[k, P, is_sparse]:
    """Copy the quotient of a checked division.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        An independent coefficient array.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    comptime if E == DivisionResult[C]:
        return _dense_output[k, P, is_sparse](rebind[DivisionResult[C]](result._data[].proof).quotient, budget)
    else:
        comptime assert E == ListDivision[C], "this result does not contain a division"
        ref proof = rebind[ListDivision[C]](result._data[].proof)
        if len(proof.quotients) != 1:
            raise PolyError(PolyErrorKind.DOMAIN, "use quotients for division by several polynomials")
        return _output[k, P, is_sparse](proof.quotients[0], budget)


def quotient[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Array[k, P, is_sparse]:
    """Copy the quotient of a checked division.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        An independent coefficient array.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return quotient(result, budget)


def remainder[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> _Array[k, P, is_sparse]:
    """Copy the remainder of a checked division.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        An independent coefficient array.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    comptime if E == DivisionResult[C]:
        return _dense_output[k, P, is_sparse](rebind[DivisionResult[C]](result._data[].proof).remainder, budget)
    else:
        comptime assert E == ListDivision[C], "this result does not contain a division"
        ref proof = rebind[ListDivision[C]](result._data[].proof)
        return _output[k, P, is_sparse](proof.remainder, budget)


def remainder[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Array[k, P, is_sparse]:
    """Copy the remainder of a checked division.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        An independent coefficient array.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return remainder(result, budget)


def quotients[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> List[_Array[k, P, is_sparse]]:
    """Copy all quotients in the original divisor order.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Independent quotient arrays in divisor order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    comptime assert E == ListDivision[C], "quotients requires a division-by-list result"
    ref proof = rebind[ListDivision[C]](result._data[].proof)
    _allocation[_Array[k, P, is_sparse]](len(proof.quotients), budget)
    var result_values = List[_Array[k, P, is_sparse]](capacity=len(proof.quotients))
    for q in proof.quotients:
        result_values.append(_output[k, P, is_sparse](q, budget))
    return result_values^


def quotients[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[_Array[k, P, is_sparse]]:
    """Copy all quotients in the original divisor order.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Independent quotient arrays in divisor order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return quotients(result, budget)


def cofactors[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> Tuple[_Array[k, P, is_sparse], _Array[k, P, is_sparse]]:
    """Copy GCD cofactors or resultant identity multipliers.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        The GCD input quotients, or the multipliers u and v in u*a + v*b = resultant.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    comptime assert E == GcdResult[C] or E == ResultantResult[C], "this result has no cofactor pair"
    comptime if E == GcdResult[C]:
        ref proof = rebind[GcdResult[C]](result._data[].proof)
        return _output[k, P, is_sparse](proof.cofactors[0], budget), _output[k, P, is_sparse](proof.cofactors[1], budget)
    else:
        ref proof = rebind[ResultantResult[C]](result._data[].proof)
        return _output[k, P, is_sparse](proof.cofactors[0], budget), _output[k, P, is_sparse](proof.cofactors[1], budget)


def cofactors[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Tuple[_Array[k, P, is_sparse], _Array[k, P, is_sparse]]:
    """Copy GCD cofactors or resultant identity multipliers.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        The GCD input quotients, or the multipliers u and v in u*a + v*b = resultant.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return cofactors(result, budget)


def unit[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> _Scalar[k]:
    """Copy the factorization scalar.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        The scalar multiplier.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    comptime if E == Factorization[C]:
        ref proof = rebind[Factorization[C]](result._data[].proof)
        _charge_size(proof.unit.size(), budget)
        return _unwrap[k](proof.unit)
    elif E == MultivariateFactorization[C]:
        ref proof = rebind[MultivariateFactorization[C]](result._data[].proof)
        _charge_size(proof.unit.size(), budget)
        return _unwrap[k](proof.unit)
    elif E == SparseSquarefreeDecomposition[C]:
        ref proof = rebind[SparseSquarefreeDecomposition[C]](result._data[].proof)
        _charge_size(proof.unit.size(), budget)
        return _unwrap[k](proof.unit)
    elif E == SquarefreeDecomposition[C]:
        ref proof = rebind[SquarefreeDecomposition[C]](result._data[].proof)
        _charge_size(proof.unit.size(), budget)
        return _unwrap[k](proof.unit)
    else:
        comptime assert False, "this result has no factorization"


def unit[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Scalar[k]:
    """Copy the factorization scalar.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        The scalar multiplier.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return unit(result, budget)


def factors[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], mut budget: Budget) raises PolyError -> List[Tuple[_Array[k, P, is_sparse], Int]]:
    """Copy factor polynomials and their positive multiplicities.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Independent factor arrays in deterministic core order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    comptime C = _Ring[k]
    comptime if E == Factorization[C]:
        ref proof = rebind[Factorization[C]](result._data[].proof)
        _allocation[Tuple[_Array[k, P, is_sparse], Int]](len(proof.factors), budget)
        var values = List[Tuple[_Array[k, P, is_sparse], Int]](capacity=len(proof.factors))
        for factor in proof.factors:
            values.append((_dense_output[k, P, is_sparse](factor[0], budget), factor[1]))
        return values^
    elif E == MultivariateFactorization[C]:
        ref proof = rebind[MultivariateFactorization[C]](result._data[].proof)
        _allocation[Tuple[_Array[k, P, is_sparse], Int]](len(proof.factors), budget)
        var values = List[Tuple[_Array[k, P, is_sparse], Int]](capacity=len(proof.factors))
        for factor in proof.factors:
            values.append((_output[k, P, is_sparse](factor[0], budget), factor[1]))
        return values^
    elif E == SparseSquarefreeDecomposition[C]:
        ref proof = rebind[SparseSquarefreeDecomposition[C]](result._data[].proof)
        _allocation[Tuple[_Array[k, P, is_sparse], Int]](len(proof.factors), budget)
        var values = List[Tuple[_Array[k, P, is_sparse], Int]](capacity=len(proof.factors))
        for factor in proof.factors:
            values.append((_output[k, P, is_sparse](factor[0], budget), factor[1]))
        return values^
    elif E == SquarefreeDecomposition[C]:
        ref proof = rebind[SquarefreeDecomposition[C]](result._data[].proof)
        _allocation[Tuple[_Array[k, P, is_sparse], Int]](len(proof.factors), budget)
        var values = List[Tuple[_Array[k, P, is_sparse], Int]](capacity=len(proof.factors))
        for factor in proof.factors:
            values.append((_dense_output[k, P, is_sparse](factor[0], budget), factor[1]))
        return values^
    else:
        comptime assert False, "this result has no factorization"


def factors[k: Int, P: Copyable & Deinitable, is_sparse: Bool, E: Copyable & Deinitable](result: Result[k, P, is_sparse, E], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> List[Tuple[_Array[k, P, is_sparse], Int]]:
    """Copy factor polynomials and their positive multiplicities.

    Parameters:
        k: Inferred from the input types unless explicitly supplied.
        P: Inferred from the input types unless explicitly supplied.
        is_sparse: Inferred from the input types unless explicitly supplied.
        E: Inferred from the input types unless explicitly supplied.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Independent factor arrays in deterministic core order.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return factors(result, budget)
