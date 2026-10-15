"""Normalized rational functions from coefficient arrays and sparse terms."""
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink, BallDomain, QQ
from ._common import _grow
from ._io_meta import _Record, _Fraction
from std.iter import Iterable, iter, next, StopIteration
from std.memory import ArcPointer
from polynomials import ModInt, Modulus
from ._common import _WithModulus
from polynomials import Budget, PolyError, PolyErrorKind, NoContext, DensePoly, SparsePoly, RationalFunction, CancellationCertificate, CoprimalityCertificate, coprimality_certificate, check_coprimality
from ._rational_field import rational_derivative, rational_compose, rational_substitute, rational_evaluate, partial_fractions as _partial_fractions, check_partial_fractions, PartialFractions as CorePartialFractions, PartialFractionTerm, laurent_expansion, check_laurent_expansion, LaurentExpansion, rational_residue, rational_poles as _rational_poles
from . import dense, sparse, _rational_field as _field
from ._algebraic_value import Real, Complex
from polynomials import QQi, CoefficientField
from ._common import _Scalar, _Ring, _CoreKind, _FieldKind, _Join, _TupleKind, _wrap, _unwrap, _domain

comptime _Field[k: Int] = Complex if k == 5 else Real if k == 4 else ModInt if k == 3 else QQi if k == 2 else QQ
comptime _Proof[k: Int] = CancellationCertificate if k == 1 else _field.Cancellation[_Field[k]]
comptime _Export[k: Int] = RationalFunction if k == 1 else Tuple[SparsePoly[_Field[k]], SparsePoly[_Field[k]]]
from ._common import _Kind, _Element, _Nested, _convert, _allocation, _charge_size, _compact, _sequence_coefficients, _values_from_core, Values
from ._interop import _dense_sparse, _dense_copy, _sparse_copy, _variables, _Shape, _take_dense, _take_sparse


@fieldwise_init
struct _Data[k: Int](Copyable, Movable):
    var value: _field.RationalFunction[_Field[Self.k]]
    var original_numerator: SparsePoly[_Field[Self.k]]
    var original_denominator: SparsePoly[_Field[Self.k]]
    var cancellation: _Proof[Self.k]
    var coprimality: CoprimalityCertificate[_Field[Self.k]]


@fieldwise_init
struct Function[P: Copyable & Deinitable, k: Int = 1](Copyable, Movable, _Record, _WithModulus):
    """A normalized rational function sharing immutable exact storage.

    Parameters:
        P: The inferred sparse power shape; Int for dense univariate input.
        k: The inferred coefficient field; users normally leave this implicit.
    """
    comptime _Format = _Fraction
    comptime _Item = _Scalar[Self.k]
    comptime _Power = Self.P
    var _data: ArcPointer[_Data[Self.k]]

    def ndim(self) -> Int:
        """Return the number of polynomial axes.

        Returns:
            The retained power-shape width, including unused axes.
        """
        return sparse._width[Self.P]()

    def domain(self) -> String:
        """Return the coefficient domain.

        Returns:
            The inferred field: QQ, QQi, GF(p), AA, or QQbar.
        """
        return _domain[Self.k](self.modulus())

    def modulus(self) -> UInt64:
        """Return the coefficient modulus.

        Returns:
            The finite-field prime, or zero in characteristic zero.
        """
        comptime if Self.k == 3:
            return rebind[Modulus](self._data[].value.numerator.ring).value
        else:
            return 0

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].cancellation))] _Proof[Self.k]:
        """Borrow the cancellation identity retained during construction.

        Returns:
            Immutable evidence relating the stored input pair to the reduced pair.
        """
        return self._data[].cancellation

    def coprimality_evidence(self) -> ref[ImmOrigin(origin_of(self._data[].coprimality))] CoprimalityCertificate[_Field[Self.k]]:
        """Borrow the proof that the returned numerator and denominator are coprime.

        Returns:
            Immutable core evidence for check_coprimality on the reduced pair.
        """
        return self._data[].coprimality


def _normalize[P: Copyable & Deinitable, k: Int](var numerator: SparsePoly[_Ring[k]], var denominator: SparsePoly[_Ring[k]], mut budget: Budget) raises PolyError -> Function[P, k]:
    comptime assert k >= 1 and k <= 5, "rational functions require exact field coefficients"
    var p = rebind_var[SparsePoly[_Field[k]]](numerator^)
    var q = rebind_var[SparsePoly[_Field[k]]](denominator^)
    var common = SparsePoly[_Field[k]].zero(p.order, p.ring, budget)
    var value = _field.RationalFunction[_Field[k]]._reduce(p, q, common, budget)
    _compact(common.terms, budget)
    var a = value.numerator.reorder(value.numerator.order, budget)
    var b = value.denominator.reorder(value.denominator.order, budget)
    _compact(a.terms, budget)
    _compact(b.terms, budget)
    var proof: _Proof[k]
    comptime if k == 1:
        proof = rebind_var[_Proof[k]](CancellationCertificate(rebind_var[SparsePoly[QQ]](common^), (rebind_var[SparsePoly[QQ]](a^), rebind_var[SparsePoly[QQ]](b^))))
    else:
        proof = rebind_var[_Proof[k]](_field.Cancellation[_Field[k]](common^, (a^, b^)))
    _compact(p.terms, budget)
    _compact(q.terms, budget)
    _compact(value.numerator.terms, budget)
    _compact(value.denominator.terms, budget)
    var coprimality = coprimality_certificate(value.numerator, value.denominator, budget)
    if not coprimality:
        raise PolyError(PolyErrorKind.DOMAIN, "normalized rational function still has a common factor")
    _allocation[_Data[k]](1, budget)
    return Function[P, k](ArcPointer(_Data[k](value^, p^, q^, proof^, coprimality.take())))


def _from_value[P: Copyable & Deinitable, k: Int](var value: _field.RationalFunction[_Field[k]], mut budget: Budget) raises PolyError -> Function[P, k]:
    return _normalize[P, k](rebind_var[SparsePoly[_Ring[k]]](_take_sparse(value.numerator)), rebind_var[SparsePoly[_Ring[k]]](_take_sparse(value.denominator)), budget)


def _promote[k: Int, P: Copyable & Deinitable, source: Int](value: Function[P, source], mut budget: Budget) raises PolyError -> Function[P, k]:
    comptime if k == source:
        budget.charge(1)
        return rebind_var[Function[P, k]](value.copy())
    elif k == 3 or source == 3:
        raise PolyError(PolyErrorKind.DOMAIN, "finite-field rational functions cannot mix with characteristic-zero functions")
    else:
        var p = sparse._promote[k](numerator(value, budget), rebind_var[_Ring[k].Context](NoContext()), budget)
        var q = sparse._promote[k](denominator(value, budget), rebind_var[_Ring[k].Context](NoContext()), budget)
        return _normalize[P, k](p._core[].reorder(p._core[].order, budget), q._core[].reorder(q._core[].order, budget), budget)


def asarray[A: Iterable, B: Iterable](numerator: A, denominator: B, mut budget: Budget) raises PolyError -> Function[Int, _FieldKind[_Join[_Kind[_Element[A]], _Kind[_Element[B]]]]]:
    """Construct and normalize a fraction from two ascending coefficient iterables.

    Parameters:
        A: The numerator iterable type.
        B: The denominator iterable type.

    Args:
        numerator: Ascending exact characteristic-zero coefficients.
        denominator: Ascending coefficients of a nonzero polynomial.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A univariate rational function over the inferred field with a monic nonzero denominator.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime k = _FieldKind[_Join[_Kind[_Element[A]], _Kind[_Element[B]]]]
    comptime assert k >= 1 and k <= 5, "rational functions require exact field coefficients"
    var pair = dense._pair[k](numerator, denominator, budget)
    return _normalize[Int, k](_dense_sparse(pair[0]._core[], budget), _dense_sparse(pair[1]._core[], budget), budget)


def asarray[A: Iterable, B: Iterable](numerator: A, denominator: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[Int, _FieldKind[_Join[_Kind[_Element[A]], _Kind[_Element[B]]]]]:
    """Construct and normalize a fraction from two ascending coefficient iterables.

    Parameters:
        A: The numerator iterable type.
        B: The denominator iterable type.

    Args:
        numerator: Ascending exact characteristic-zero coefficients.
        denominator: Ascending coefficients of a nonzero polynomial.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A univariate rational function over the inferred field with a monic nonzero denominator.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return asarray(numerator, denominator, budget)


def asarray[I: Iterable](values: I, mut budget: Budget) raises PolyError -> Function[Int, _FieldKind[_Kind[_Element[I]]]]:
    """Embed a coefficient polynomial as a rational function with denominator one.

    Parameters:
        I: The coefficient iterable type.

    Args:
        values: Ascending exact characteristic-zero coefficients.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A normalized univariate rational function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    return rebind_var[Function[Int, _FieldKind[_Kind[_Element[I]]]]](asarray(values, [1], budget))


def asarray[I: Iterable](values: I, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[Int, _FieldKind[_Kind[_Element[I]]]]:
    """Embed a coefficient polynomial as a rational function with denominator one.

    Parameters:
        I: The coefficient iterable type.

    Args:
        values: Ascending exact characteristic-zero coefficients.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A normalized univariate rational function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return asarray(values, budget)


def asarray[P: Copyable & Deinitable, k: Int](values: Function[P, k], mut budget: Budget) raises PolyError -> Function[P, k]:
    """Share an already prepared rational function.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        values: The prepared rational function.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A value sharing the immutable data and cancellation evidence.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    budget.charge(1)
    return values.copy()


def asarray[P: Copyable & Deinitable, k: Int](values: Function[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, k]:
    """Share an already prepared rational function.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        values: The prepared rational function.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A value sharing the immutable data and cancellation evidence.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return asarray(values, budget)


def from_terms[A: Iterable, B: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, D: Copyable & Deinitable](numerator: A, denominator: B, mut budget: Budget) raises PolyError -> Function[P, _FieldKind[_Join[_Kind[C], _Kind[D]]]] where _Element[A] == Tuple[P, C] and _Element[B] == Tuple[P, D]:
    """Construct a rational function from sparse numerator and denominator terms.

    Parameters:
        A: The numerator iterable type.
        B: The denominator iterable type.
        P: The common power shape.
        C: The numerator coefficient type.
        D: The denominator coefficient type.

    Args:
        numerator: Exact (power, coefficient) terms.
        denominator: Nonzero denominator terms with the same power shape.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A normalized fraction retaining the input power shape.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime k = _FieldKind[_Join[_Kind[C], _Kind[D]]]
    comptime assert k >= 1 and k <= 5, "rational functions require exact field coefficients"
    var pair = sparse._pair[k, A, B, P, C, D](numerator, denominator, budget)
    return _normalize[P, k](pair[0]._core[].reorder(pair[0]._core[].order, budget), pair[1]._core[].reorder(pair[1]._core[].order, budget), budget)


def from_terms[A: Iterable, B: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, D: Copyable & Deinitable](numerator: A, denominator: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, _FieldKind[_Join[_Kind[C], _Kind[D]]]] where _Element[A] == Tuple[P, C] and _Element[B] == Tuple[P, D]:
    """Construct a rational function from sparse numerator and denominator terms.

    Parameters:
        A: The numerator iterable type.
        B: The denominator iterable type.
        P: The common power shape.
        C: The numerator coefficient type.
        D: The denominator coefficient type.

    Args:
        numerator: Exact (power, coefficient) terms.
        denominator: Nonzero denominator terms with the same power shape.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A normalized fraction retaining the input power shape.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return from_terms[A, B, P, C, D](numerator, denominator, budget)


def numerator[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget) raises PolyError -> sparse.Array[P, k]:
    """Extract the normalized numerator as a sparse array.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The prepared rational function.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        An immutable polynomial with the original power shape.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    ref p = value._data[].value.numerator
    return sparse.Array[P, k](rebind_var[SparsePoly[_Ring[k]]](p.reorder(p.order, budget)), budget)


def numerator[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> sparse.Array[P, k]:
    """Extract the normalized numerator as a sparse array.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The prepared rational function.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An immutable polynomial with the original power shape.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return numerator(value, budget)


def denominator[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget) raises PolyError -> sparse.Array[P, k]:
    """Extract the normalized denominator as a sparse array.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The prepared rational function.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        An immutable polynomial with the original power shape.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    ref p = value._data[].value.denominator
    return sparse.Array[P, k](rebind_var[SparsePoly[_Ring[k]]](p.reorder(p.order, budget)), budget)


def denominator[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> sparse.Array[P, k]:
    """Extract the normalized denominator as a sparse array.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The prepared rational function.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An immutable polynomial with the original power shape.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return denominator(value, budget)


def add[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], mut budget: Budget) raises PolyError -> Function[P, _Join[k, l]]:
    """Add rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _Join[k, l]
    var left = _promote[joined](a, budget)
    var right = _promote[joined](b, budget)
    ref p = left._data[].value
    ref q = right._data[].value
    var numerator = p.numerator.mul(q.denominator, budget).add(q.numerator.mul(p.denominator, budget), budget)
    var denominator = p.denominator.mul(q.denominator, budget)
    return _normalize[P, joined](rebind_var[SparsePoly[_Ring[joined]]](numerator^), rebind_var[SparsePoly[_Ring[joined]]](denominator^), budget)


def add[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, _Join[k, l]]:
    """Add rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return add(a, b, budget)


def subtract[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], mut budget: Budget) raises PolyError -> Function[P, _Join[k, l]]:
    """Subtract rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _Join[k, l]
    var left = _promote[joined](a, budget)
    var right = _promote[joined](b, budget)
    ref p = left._data[].value
    ref q = right._data[].value
    var numerator = p.numerator.mul(q.denominator, budget).sub(q.numerator.mul(p.denominator, budget), budget)
    var denominator = p.denominator.mul(q.denominator, budget)
    return _normalize[P, joined](rebind_var[SparsePoly[_Ring[joined]]](numerator^), rebind_var[SparsePoly[_Ring[joined]]](denominator^), budget)


def subtract[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, _Join[k, l]]:
    """Subtract rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return subtract(a, b, budget)


def multiply[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], mut budget: Budget) raises PolyError -> Function[P, _Join[k, l]]:
    """Multiply rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _Join[k, l]
    var left = _promote[joined](a, budget)
    var right = _promote[joined](b, budget)
    ref p = left._data[].value
    ref q = right._data[].value
    var numerator = p.numerator.mul(q.numerator, budget)
    var denominator = p.denominator.mul(q.denominator, budget)
    return _normalize[P, joined](rebind_var[SparsePoly[_Ring[joined]]](numerator^), rebind_var[SparsePoly[_Ring[joined]]](denominator^), budget)


def multiply[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, _Join[k, l]]:
    """Multiply rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return multiply(a, b, budget)


def divide[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], mut budget: Budget) raises PolyError -> Function[P, _Join[k, l]]:
    """Divide rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _Join[k, l]
    var left = _promote[joined](a, budget)
    var right = _promote[joined](b, budget)
    ref p = left._data[].value
    ref q = right._data[].value
    var numerator = p.numerator.mul(q.denominator, budget)
    var denominator = p.denominator.mul(q.numerator, budget)
    return _normalize[P, joined](rebind_var[SparsePoly[_Ring[joined]]](numerator^), rebind_var[SparsePoly[_Ring[joined]]](denominator^), budget)


def divide[P: Copyable & Deinitable, k: Int, l: Int](a: Function[P, k], b: Function[P, l], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, _Join[k, l]]:
    """Divide rational functions and normalize the resulting fraction.

    Parameters:
        P: The shared polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        a: The left rational function.
        b: The right rational function; division requires a nonzero numerator.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new fraction retaining its cancellation identity.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return divide(a, b, budget)


def negative[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget) raises PolyError -> Function[P, k]:
    """Negate a rational function.

    Parameters:
        P: The inferred polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The input rational function.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new normalized rational function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    ref p = value._data[].value
    return _from_value[P, k](_field.RationalFunction[_Field[k]](p.numerator.neg(budget), p.denominator.reorder(p.denominator.order, budget)), budget)


def negative[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, k]:
    """Negate a rational function.

    Parameters:
        P: The inferred polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The input rational function.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new normalized rational function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return negative(value, budget)


def reciprocal[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget) raises PolyError -> Function[P, k]:
    """Invert a nonzero rational function.

    Parameters:
        P: The inferred polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The input rational function.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new normalized rational function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    ref p = value._data[].value
    return _from_value[P, k](_field.RationalFunction[_Field[k]](p.denominator.reorder(p.denominator.order, budget), p.numerator.reorder(p.numerator.order, budget)), budget)


def reciprocal[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[P, k]:
    """Invert a nonzero rational function.

    Parameters:
        P: The inferred polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The input rational function.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new normalized rational function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return reciprocal(value, budget)


def derivative[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget, *, axis: Int = 0) raises PolyError -> Function[P, k]:
    """Differentiate a rational function along one array axis.

    Parameters:
        P: The inferred polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The input rational function.
        budget: The shared work allowance, including preparation and result storage.
        axis: The zero-based differentiation axis.

    Returns:
        The normalized exact derivative.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    sparse._axis[P](axis)
    return _from_value[P, k](rational_derivative(value._data[].value, Int64(axis), budget), budget)


def derivative[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX, axis: Int = 0) raises PolyError -> Function[P, k]:
    """Differentiate a rational function along one array axis.

    Parameters:
        P: The inferred polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The input rational function.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axis: The zero-based differentiation axis.

    Returns:
        The normalized exact derivative.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return derivative(value, budget, axis=axis)


def compose[P: Copyable & Deinitable, k: Int, l: Int](outer: Function[P, k], inner: Function[P, l], mut budget: Budget, *, axis: Int = 0) raises PolyError -> Function[P, _Join[k, l]]:
    """Compose rational functions that depend only on the selected axis.

    Parameters:
        P: The common power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        outer: The outer univariate rational function.
        inner: The univariate replacement.
        budget: The shared work allowance, including preparation and result storage.
        axis: The sole occurring variable axis.

    Returns:
        The normalized composition; other occurring axes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _Join[k, l]
    var left = _promote[joined](outer, budget)
    var right = _promote[joined](inner, budget)
    sparse._axis[P](axis)
    return _from_value[P, joined](rational_compose(left._data[].value, right._data[].value, Int64(axis), budget), budget)


def compose[P: Copyable & Deinitable, k: Int, l: Int](outer: Function[P, k], inner: Function[P, l], *, max_steps: UInt64 = UInt64.MAX, axis: Int = 0) raises PolyError -> Function[P, _Join[k, l]]:
    """Compose rational functions that depend only on the selected axis.

    Parameters:
        P: The common power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        l: The second function's inferred coefficient field.

    Args:
        outer: The outer univariate rational function.
        inner: The univariate replacement.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axis: The sole occurring variable axis.

    Returns:
        The normalized composition; other occurring axes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return compose(outer, inner, budget, axis=axis)


def evaluate[P: Copyable & Deinitable, k: Int, *Ts: Copyable & Deinitable](value: Function[P, k], point: Tuple[*Ts], mut budget: Budget) raises PolyError -> _Scalar[_FieldKind[_Join[k, _TupleKind[Tuple[*Ts]]]]]:
    """Evaluate the normalized function at exact coordinates.

    Parameters:
        P: The polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        Ts: The inferred coordinate types.

    Args:
        value: The input rational function.
        point: One exact coordinate per array axis; coordinates can promote the result field.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        An exact value; a remaining pole raises DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _FieldKind[_Join[k, _TupleKind[Tuple[*Ts]]]]
    comptime assert joined >= 1 and joined <= 5, "rational evaluation requires exact field coordinates"
    comptime assert Tuple[*Ts].Ts.length == sparse._width[P](), "rational evaluation point has the wrong number of coordinates"
    var prepared = _promote[joined](value, budget)
    _allocation[Tuple[Int64, _Field[joined]]](Tuple[*Ts].Ts.length, budget)
    var assignments = List[Tuple[Int64, _Field[joined]]](capacity=Tuple[*Ts].Ts.length)
    comptime for i in range(Tuple[*Ts].Ts.length):
        assignments.append((Int64(i), rebind_var[_Field[joined]](_wrap[joined](_convert[joined](point[i], rebind_var[_Ring[joined].Context](prepared._data[].value.numerator.ring.copy()), budget)))))
    var result = rational_evaluate(prepared._data[].value, assignments, budget)
    _charge_size(result.size(), budget)
    return _unwrap[joined](rebind[_Ring[joined]](result))


def evaluate[P: Copyable & Deinitable, k: Int, *Ts: Copyable & Deinitable](value: Function[P, k], point: Tuple[*Ts], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Scalar[_FieldKind[_Join[k, _TupleKind[Tuple[*Ts]]]]]:
    """Evaluate the normalized function at exact coordinates.

    Parameters:
        P: The polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        Ts: The inferred coordinate types.

    Args:
        value: The input rational function.
        point: One exact coordinate per array axis; coordinates can promote the result field.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An exact value; a remaining pole raises DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return evaluate[P, k, *Ts](value, point, budget)


def evaluate[P: Copyable & Deinitable, k: Int, T: Copyable & Deinitable](value: Function[P, k], point: T, mut budget: Budget) raises PolyError -> _Scalar[_FieldKind[_Join[k, _Kind[T]]]] where _Kind[T] >= 0:
    """Evaluate a univariate rational function at an exact point.

    Parameters:
        P: The one-variable power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        T: The exact scalar point type.

    Args:
        value: The input rational function.
        point: The scalar evaluation point.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        An exact value; removable factors have already been cancelled.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime assert sparse._width[P]() == 1, "multivariate rational evaluation needs a coordinate tuple"
    return rebind_var[_Scalar[_FieldKind[_Join[k, _Kind[T]]]]](evaluate(value, (point.copy(),), budget))


def evaluate[P: Copyable & Deinitable, k: Int, T: Copyable & Deinitable](value: Function[P, k], point: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Scalar[_FieldKind[_Join[k, _Kind[T]]]] where _Kind[T] >= 0:
    """Evaluate a univariate rational function at an exact point.

    Parameters:
        P: The one-variable power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        T: The exact scalar point type.

    Args:
        value: The input rational function.
        point: The scalar evaluation point.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An exact value; removable factors have already been cancelled.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return evaluate(value, point, budget)


def substitute[P: Copyable & Deinitable, k: Int, I: Iterable, C: Copyable & Deinitable](value: Function[P, k], replacements: I, mut budget: Budget, *, axes: List[Int] = List[Int]()) raises PolyError -> Function[P, _FieldKind[_Join[k, _Kind[C]]]] where _Element[_Nested[I]] == Tuple[P, C]:
    """Substitute sparse polynomials simultaneously into numerator and denominator.

    Parameters:
        P: The common power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        I: The iterable of replacement polynomial iterables.
        C: The replacement coefficient type.

    Args:
        value: The original rational function.
        replacements: A finite iterable of sparse polynomial inputs with the same power shape.
        budget: The shared work allowance, including preparation and result storage.
        axes: Distinct axes corresponding to replacements; empty assigns every axis.

    Returns:
        A new normalized function; an identically zero denominator raises DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _FieldKind[_Join[k, _Kind[C]]]
    comptime assert joined >= 1 and joined <= 5, "rational substitutions require exact field coefficients"
    var source = _promote[joined](value, budget)
    var original = numerator(source, budget)
    var prepared = sparse._many[joined, sparse.Array[P, joined], I, P, _Scalar[joined], C](original, replacements, budget)
    var count = len(prepared[1])
    if (len(axes) and len(axes) != count) or (not len(axes) and count != sparse._width[P]()):
        raise PolyError(PolyErrorKind.DOMAIN, "provide one replacement per axis, or matching explicit axes")
    _allocation[Tuple[Int64, SparsePoly[_Field[joined]]]](count, budget)
    var assignments = List[Tuple[Int64, SparsePoly[_Field[joined]]]](capacity=count)
    for i in range(count):
        budget.charge(1)
        var axis = axes[i] if len(axes) else i
        sparse._axis[P](axis)
        for j in range(i):
            budget.charge(1)
            if (axes[j] if len(axes) else j) == axis:
                raise PolyError(PolyErrorKind.DOMAIN, "substitution axes must be distinct")
        ref p = prepared[1][i]._core[]
        assignments.append((Int64(axis), rebind_var[SparsePoly[_Field[joined]]](p.reorder(p.order, budget))))
    return _from_value[P, joined](rational_substitute(source._data[].value, assignments, budget), budget)


def substitute[P: Copyable & Deinitable, k: Int, I: Iterable, C: Copyable & Deinitable](value: Function[P, k], replacements: I, *, max_steps: UInt64 = UInt64.MAX, axes: List[Int] = List[Int]()) raises PolyError -> Function[P, _FieldKind[_Join[k, _Kind[C]]]] where _Element[_Nested[I]] == Tuple[P, C]:
    """Substitute sparse polynomials simultaneously into numerator and denominator.

    Parameters:
        P: The common power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        I: The iterable of replacement polynomial iterables.
        C: The replacement coefficient type.

    Args:
        value: The original rational function.
        replacements: A finite iterable of sparse polynomial inputs with the same power shape.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axes: Distinct axes corresponding to replacements; empty assigns every axis.

    Returns:
        A new normalized function; an identically zero denominator raises DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return substitute[P, k, I, C](value, replacements, budget, axes=axes)


def check[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget) raises PolyError -> Bool:
    """Verify complete cancellation and the monic denominator of the returned fraction.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The fraction and its retained normalization evidence.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        Whether the retained identity, coprimality proof, and denominator normalization all hold.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    ref data = value._data[]
    if not check_coprimality(data.value.numerator, data.value.denominator, data.coprimality, budget):
        return False
    comptime if k == 1:
        ref proof = rebind[CancellationCertificate](data.cancellation)
        return _check_fraction(data.value, data.original_numerator, data.original_denominator, rebind[SparsePoly[_Field[k]]](proof.common), rebind[SparsePoly[_Field[k]]](proof.cofactors[0]), rebind[SparsePoly[_Field[k]]](proof.cofactors[1]), budget)
    else:
        ref proof = rebind[_field.Cancellation[_Field[k]]](data.cancellation)
        return _check_fraction(data.value, data.original_numerator, data.original_denominator, rebind[SparsePoly[_Field[k]]](proof.common), rebind[SparsePoly[_Field[k]]](proof.cofactors[0]), rebind[SparsePoly[_Field[k]]](proof.cofactors[1]), budget)



def check[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify complete cancellation and the monic denominator of the returned fraction.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The fraction and its retained normalization evidence.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the retained identity, coprimality proof, and denominator normalization all hold.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return check(value, budget)


def to_core[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget, *, variables: List[Int64] = List[Int64]()) raises PolyError -> _Export[k]:
    """Export independent mutable core numerator and denominator storage.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The fraction to export.
        budget: The shared work allowance, including preparation and result storage.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.

    Returns:
        A core RationalFunction for QQ; otherwise a pair of public SparsePoly values using the selected core variable IDs.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    ref p = value._data[].value
    var mapping = _variables(sparse._width[P](), variables, budget)
    comptime if k == 1:
        return rebind_var[_Export[k]](RationalFunction(rebind_var[SparsePoly[QQ]](_sparse_copy(p.numerator, mapping, True, budget)), rebind_var[SparsePoly[QQ]](_sparse_copy(p.denominator, mapping, True, budget))))
    else:
        return rebind_var[_Export[k]]((_sparse_copy(p.numerator, mapping, True, budget), _sparse_copy(p.denominator, mapping, True, budget)))


def to_core[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX, variables: List[Int64] = List[Int64]()) raises PolyError -> _Export[k]:
    """Export independent mutable core numerator and denominator storage.

    Parameters:
        P: The inferred power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The fraction to export.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.

    Returns:
        A core RationalFunction for QQ; otherwise a pair of public SparsePoly values using the selected core variable IDs.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return to_core(value, budget, variables=variables)


def from_core[ndim: Int](value: RationalFunction, mut budget: Budget, *, variables: List[Int64] = List[Int64]()) raises PolyError -> Function[_Shape[ndim]]:
    """Validate, remap, and normalize a mutable core rational function.

    Parameters:
        ndim: The explicit array power-shape width.

    Args:
        value: The core numerator and denominator.
        budget: The shared work allowance, including preparation and result storage.
        variables: Core variable IDs in array-axis order; empty uses 0 through ndim-1.

    Returns:
        An independent fraction retaining every specified array axis.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime assert ndim >= 0, "a rational function needs at least one array axis"
    var mapping = _variables(ndim, variables, budget)
    return _normalize[_Shape[ndim], 1](_sparse_copy(value.numerator, mapping, False, budget), _sparse_copy(value.denominator, mapping, False, budget), budget)


def from_core[ndim: Int](value: RationalFunction, *, max_steps: UInt64 = UInt64.MAX, variables: List[Int64] = List[Int64]()) raises PolyError -> Function[_Shape[ndim]]:
    """Validate, remap, and normalize a mutable core rational function.

    Parameters:
        ndim: The explicit array power-shape width.

    Args:
        value: The core numerator and denominator.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        variables: Core variable IDs in array-axis order; empty uses 0 through ndim-1.

    Returns:
        An independent fraction retaining every specified array axis.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return from_core[ndim](value, budget, variables=variables)


@fieldwise_init
struct PartialFractions[P: Copyable & Deinitable, k: Int = 1]( _Archive, Copyable, Movable):
    """Polynomial part and proper fraction terms, with the original function retained.

    Parameters:
        P: The original polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
    """
    var _source: Function[Self.P, Self.k]
    var _axis: Int
    var _factors: ArcPointer[_field._Factors[_Field[Self.k]]]
    var polynomial: dense.Array[Self.k]
    """The polynomial part in the selected variable."""
    var terms: Values[Tuple[dense.Array[Self.k], dense.Array[Self.k], Int], Self.k]
    """Each entry is (numerator, denominator factor, positive exponent)."""



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
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[Function[Self.P, Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._source, writer, sink, budget, depth + 1)
        _result_io._write_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._axis, writer, sink, budget, depth + 1)
        _result_io._write_node[_field._Factors[_Field[Self.k]], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._factors[], writer, sink, budget, depth + 1)
        _result_io._write_node[dense.Array[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self.polynomial, writer, sink, budget, depth + 1)
        writer.line(String(len(self.terms)), sink, budget)
        for field714 in self.terms:
            _result_io._write_node[dense.Array[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](field714[0], writer, sink, budget, depth + 1)
            _result_io._write_node[dense.Array[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](field714[1], writer, sink, budget, depth + 1)
            _result_io._write_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](field714[2], writer, sink, budget, depth + 1)
        writer.line(String(self.terms.modulus()), sink, budget)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field715 = _result_io._read_node[Function[Self.P, Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field716 = _result_io._read_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field718 = _result_io._read_node[_field._Factors[_Field[Self.k]], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_field._Factors[_Field[Self.k]]](1, budget)
        var field717 = ArcPointer(field718^)
        var field719 = _result_io._read_node[dense.Array[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field722 = reader.count(budget)
        var field721 = List[Tuple[dense.Array[Self.k], dense.Array[Self.k], Int]]()
        for _ in range(field722):
            var field724 = _result_io._read_node[dense.Array[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            var field725 = _result_io._read_node[dense.Array[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            var field726 = _result_io._read_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            var field723 = (field724^, field725^, field726)
            _grow(field721, budget)
            field721.append(field723^)
        var field720 = _common.Values[Tuple[dense.Array[Self.k], dense.Array[Self.k], Int], Self.k](field721^, reader.number(UInt64.MAX, budget), budget)
        return Self(field715^, field716, field717^, field719^, field720^)

def partial_fractions[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget, *, axis: Int = 0) raises PolyError -> PartialFractions[P, k]:
    """Decompose a univariate fraction into a polynomial part and proper terms.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: A fraction depending only on the selected axis.
        budget: The shared work allowance, including preparation and result storage.
        axis: The sole occurring variable axis.

    Returns:
        Exact terms over distinct monic irreducible factors, retaining evidence for check.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    sparse._axis[P](axis)
    var proof = _partial_fractions(value._data[].value, Int64(axis), budget)
    _allocation[Tuple[dense.Array[k], dense.Array[k], Int]](len(proof.terms), budget)
    var terms = List[Tuple[dense.Array[k], dense.Array[k], Int]](capacity=len(proof.terms))
    for i in range(len(proof.terms)):
        var p = _take_dense(proof.terms[i].numerator)
        var q = _take_dense(proof.terms[i].factor)
        p.variable = 0
        q.variable = 0
        terms.append((dense.Array[k](rebind_var[DensePoly[_Ring[k]]](p^), budget), dense.Array[k](rebind_var[DensePoly[_Ring[k]]](q^), budget), proof.terms[i].exponent))
    var polynomial = _take_dense(proof.polynomial)
    polynomial.variable = 0
    return PartialFractions[P, k](value.copy(), axis, proof.factors.copy(), dense.Array[k](rebind_var[DensePoly[_Ring[k]]](polynomial^), budget), Values[Tuple[dense.Array[k], dense.Array[k], Int], k](terms^, value.modulus(), budget))


def partial_fractions[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX, axis: Int = 0) raises PolyError -> PartialFractions[P, k]:
    """Decompose a univariate fraction into a polynomial part and proper terms.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: A fraction depending only on the selected axis.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axis: The sole occurring variable axis.

    Returns:
        Exact terms over distinct monic irreducible factors, retaining evidence for check.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return partial_fractions(value, budget, axis=axis)


def check[P: Copyable & Deinitable, k: Int](value: PartialFractions[P, k], mut budget: Budget) raises PolyError -> Bool:
    """Check every displayed partial-fraction term against the retained function.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The decomposition to verify.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        Whether cancellation is complete and the proper terms use irreducible factors and recompose correctly.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    sparse._axis[P](value._axis)
    if not check(value._source, budget):
        return False
    _allocation[PartialFractionTerm[_Field[k]]](len(value.terms), budget)
    var terms = List[PartialFractionTerm[_Field[k]]](capacity=len(value.terms))
    for entry in value.terms:
        terms.append(PartialFractionTerm[_Field[k]](rebind_var[DensePoly[_Field[k]]](_dense_copy(entry[0]._core[], Int64(value._axis), budget)), rebind_var[DensePoly[_Field[k]]](_dense_copy(entry[1]._core[], Int64(value._axis), budget)), entry[2]))
    var proof = CorePartialFractions[_Field[k]](rebind_var[DensePoly[_Field[k]]](_dense_copy(value.polynomial._core[], Int64(value._axis), budget)), terms^, value._factors.copy())
    return check_partial_fractions(value._source._data[].value, Int64(value._axis), proof, budget)


def check[P: Copyable & Deinitable, k: Int](value: PartialFractions[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Check every displayed partial-fraction term against the retained function.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The decomposition to verify.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether cancellation is complete and the proper terms use irreducible factors and recompose correctly.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return check(value, budget)


def poles[P: Copyable & Deinitable, k: Int](value: Function[P, k], mut budget: Budget, *, axis: Int = 0) raises PolyError -> Tuple[_Scalar[k], List[Tuple[dense.Array[k], Int]]]:
    """Represent poles by factors of the reduced denominator.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: A fraction univariate in the selected axis.
        budget: The shared work allowance, including preparation and result storage.
        axis: The sole occurring variable axis.

    Returns:
        The denominator scalar and (factor, multiplicity) entries; each factor represents its zeros over an algebraic closure of the coefficient field.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    sparse._axis[P](axis)
    var proof = _rational_poles(value._data[].value, Int64(axis), budget)
    _allocation[Tuple[dense.Array[k], Int]](len(proof.factors), budget)
    var factors = List[Tuple[dense.Array[k], Int]](capacity=len(proof.factors))
    for i in range(len(proof.factors)):
        var p = _take_dense(proof.factors[i][0])
        p.variable = 0
        factors.append((dense.Array[k](rebind_var[DensePoly[_Ring[k]]](p^), budget), proof.factors[i][1]))
    _charge_size(proof.unit.size(), budget)
    return _unwrap[k](rebind[_Ring[k]](proof.unit)), factors^


def poles[P: Copyable & Deinitable, k: Int](value: Function[P, k], *, max_steps: UInt64 = UInt64.MAX, axis: Int = 0) raises PolyError -> Tuple[_Scalar[k], List[Tuple[dense.Array[k], Int]]]:
    """Represent poles by factors of the reduced denominator.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: A fraction univariate in the selected axis.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axis: The sole occurring variable axis.

    Returns:
        The denominator scalar and (factor, multiplicity) entries; each factor represents its zeros over an algebraic closure of the coefficient field.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return poles(value, budget, axis=axis)


@fieldwise_init
struct Expansion[P: Copyable & Deinitable, k: Int = 1]( _Archive, Copyable, Movable):
    """Consecutive Laurent coefficients with their center and lowest exponent.

    Parameters:
        P: The original polynomial power shape.
        k: The inferred coefficient field; users normally leave this implicit.
    """
    var _source: Function[Self.P, Self.k]
    var _axis: Int
    var center: _Scalar[Self.k]
    """The exact expansion center."""
    var valuation: Int
    """The exponent of the first coefficient; negative at a pole."""
    var coefficients: Values[_Scalar[Self.k], Self.k]
    """Known coefficients in ascending powers beginning at valuation."""



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
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[Function[Self.P, Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._source, writer, sink, budget, depth + 1)
        _result_io._write_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._axis, writer, sink, budget, depth + 1)
        _result_io._write_node[_Scalar[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](self.center, writer, sink, budget, depth + 1)
        _result_io._write_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self.valuation, writer, sink, budget, depth + 1)
        writer.line(String(len(self.coefficients)), sink, budget)
        for field727 in self.coefficients:
            _result_io._write_node[_Scalar[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](field727, writer, sink, budget, depth + 1)
        writer.line(String(self.coefficients.modulus()), sink, budget)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = _Ring[Self.k]
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field728 = _result_io._read_node[Function[Self.P, Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field729 = _result_io._read_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field730 = _result_io._read_node[_Scalar[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field731 = _result_io._read_node[Int, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        var field734 = reader.count(budget)
        var field733 = List[_Scalar[Self.k]]()
        for _ in range(field734):
            var field735 = _result_io._read_node[_Scalar[Self.k], _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            _grow(field733, budget)
            field733.append(field735^)
        var field732 = _common.Values[_Scalar[Self.k], Self.k](field733^, reader.number(UInt64.MAX, budget), budget)
        return Self(field728^, field729, field730^, field731, field732^)

def laurent[P: Copyable & Deinitable, k: Int, T: Copyable & Deinitable](value: Function[P, k], center: T, mut budget: Budget, *, terms: Int, axis: Int = 0) raises PolyError -> Expansion[P, _FieldKind[_Join[k, _Kind[T]]]]:
    """Expand a rational function at an exact center, including a pole.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        T: The exact center type.

    Args:
        value: A fraction univariate in the selected axis.
        center: An exact expansion center.
        budget: The shared work allowance, including preparation and result storage.
        terms: Nonnegative number of coefficients beginning at the valuation.
        axis: The sole occurring variable axis.

    Returns:
        The center, valuation, and requested number of consecutive coefficients.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _FieldKind[_Join[k, _Kind[T]]]
    comptime assert joined >= 1 and joined <= 5, "expansion centers require exact field coefficients"
    var source = _promote[joined](value, budget)
    sparse._axis[P](axis)
    if terms < 0:
        raise PolyError(PolyErrorKind.DOMAIN, "Laurent term count must be nonnegative")
    var c = rebind_var[_Field[joined]](_wrap[joined](_convert[joined](center, rebind_var[_Ring[joined].Context](source._data[].value.numerator.ring.copy()), budget)))
    var proof = laurent_expansion(source._data[].value, Int64(axis), c, terms, budget)
    _charge_size(proof.center.size(), budget)
    var coefficients = List[_Field[joined]]()
    swap(coefficients, proof.coefficients)
    return Expansion[P, joined](source.copy(), axis, _unwrap[joined](rebind[_Ring[joined]](proof.center)), proof.valuation, _values_from_core[joined](rebind_var[List[_Ring[joined]]](coefficients^), rebind_var[_Ring[joined].Context](source._data[].value.numerator.ring.copy()), budget))


def laurent[P: Copyable & Deinitable, k: Int, T: Copyable & Deinitable](value: Function[P, k], center: T, *, max_steps: UInt64 = UInt64.MAX, terms: Int, axis: Int = 0) raises PolyError -> Expansion[P, _FieldKind[_Join[k, _Kind[T]]]]:
    """Expand a rational function at an exact center, including a pole.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        T: The exact center type.

    Args:
        value: A fraction univariate in the selected axis.
        center: An exact expansion center.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        terms: Nonnegative number of coefficients beginning at the valuation.
        axis: The sole occurring variable axis.

    Returns:
        The center, valuation, and requested number of consecutive coefficients.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return laurent(value, center, budget, terms=terms, axis=axis)


def check[P: Copyable & Deinitable, k: Int](value: Expansion[P, k], mut budget: Budget) raises PolyError -> Bool:
    """Check displayed Laurent coefficients by a truncated multiplication identity.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The proposed Laurent expansion.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        Whether the center, valuation, and every coefficient agree with the retained function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    sparse._axis[P](value._axis)
    if not check(value._source, budget):
        return False
    var center = rebind_var[_Field[k]](_wrap[k](_convert[k](value.center, rebind_var[_Ring[k].Context](value._source._data[].value.numerator.ring.copy()), budget)))
    var proof = LaurentExpansion[_Field[k]](Int64(value._axis), center^, value.valuation, rebind_var[List[_Field[k]]](_sequence_coefficients[k](value.coefficients, rebind_var[_Ring[k].Context](value._source._data[].value.numerator.ring.copy()), budget)))
    return check_laurent_expansion(value._source._data[].value, proof, budget)


def check[P: Copyable & Deinitable, k: Int](value: Expansion[P, k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Check displayed Laurent coefficients by a truncated multiplication identity.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.

    Args:
        value: The proposed Laurent expansion.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the center, valuation, and every coefficient agree with the retained function.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return check(value, budget)


def residue[P: Copyable & Deinitable, k: Int, T: Copyable & Deinitable](value: Function[P, k], center: T, mut budget: Budget, *, axis: Int = 0) raises PolyError -> _Scalar[_FieldKind[_Join[k, _Kind[T]]]]:
    """Return the coefficient of (x-center)**-1, including zero at a regular point.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        T: The center type.

    Args:
        value: A fraction univariate in the selected axis.
        center: An exact expansion center.
        budget: The shared work allowance, including preparation and result storage.
        axis: The sole occurring variable axis.

    Returns:
        The exact residue.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime joined = _FieldKind[_Join[k, _Kind[T]]]
    comptime assert joined >= 1 and joined <= 5, "expansion centers require exact field coefficients"
    var source = _promote[joined](value, budget)
    sparse._axis[P](axis)
    var result = rational_residue(source._data[].value, Int64(axis), rebind_var[_Field[joined]](_wrap[joined](_convert[joined](center, rebind_var[_Ring[joined].Context](source._data[].value.numerator.ring.copy()), budget))), budget)
    _charge_size(result.size(), budget)
    return _unwrap[joined](rebind[_Ring[joined]](result))


def residue[P: Copyable & Deinitable, k: Int, T: Copyable & Deinitable](value: Function[P, k], center: T, *, max_steps: UInt64 = UInt64.MAX, axis: Int = 0) raises PolyError -> _Scalar[_FieldKind[_Join[k, _Kind[T]]]]:
    """Return the coefficient of (x-center)**-1, including zero at a regular point.

    Parameters:
        P: The original power shape.
        k: The inferred coefficient field; users normally leave this implicit.
        T: The center type.

    Args:
        value: A fraction univariate in the selected axis.
        center: An exact expansion center.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        axis: The sole occurring variable axis.

    Returns:
        The exact residue.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return residue(value, center, budget, axis=axis)


def _check_fraction[C: CoefficientField](value: _field.RationalFunction[C], p: SparsePoly[C], q: SparsePoly[C], common: SparsePoly[C], a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    _field._same(p, q, budget)
    _field._same(value.numerator, value.denominator, budget)
    if q.is_zero() or value.denominator.is_zero():
        return False
    if not value.denominator.terms[len(value.denominator.terms) - 1].coefficient.sub(C.one(value.denominator.ring), budget).is_zero():
        return False
    if not value.numerator.sub(a, budget).is_zero() or not value.denominator.sub(b, budget).is_zero():
        return False
    return p.sub(common.mul(a, budget), budget).is_zero() and q.sub(common.mul(b, budget), budget).is_zero()


def from_core[ndim: Int, C: CoefficientField](value: Tuple[SparsePoly[C], SparsePoly[C]], mut budget: Budget, *, variables: List[Int64] = List[Int64]()) raises PolyError -> Function[_Shape[ndim], _CoreKind[C]]:
    """Import a pair of core polynomials over QQ, QQi, GF(p), AA, or QQbar.

    Parameters:
        ndim: The retained number of array axes, including unused axes.
        C: The inferred coefficient field.

    Args:
        value: Mutable core numerator and denominator polynomials.
        budget: The shared work allowance.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.

    Returns:
        An independent, normalized fraction with cancellation evidence.

    Raises:
        PolyError: DOMAIN for invalid inputs or a zero denominator, or LIMIT/CAPACITY on resource limits.
    """
    comptime k = _CoreKind[C]
    comptime assert k >= 1 and k <= 5, "rational functions require exact field coefficients"
    comptime assert ndim >= 0, "array axis count must be nonnegative"
    var mapping = _variables(ndim, variables, budget)
    return _normalize[_Shape[ndim], k](rebind_var[SparsePoly[_Ring[k]]](_sparse_copy(value[0], mapping, False, budget)), rebind_var[SparsePoly[_Ring[k]]](_sparse_copy(value[1], mapping, False, budget)), budget)


def from_core[ndim: Int, C: CoefficientField](value: Tuple[SparsePoly[C], SparsePoly[C]], *, variables: List[Int64] = List[Int64](), max_steps: UInt64 = UInt64.MAX) raises PolyError -> Function[_Shape[ndim], _CoreKind[C]]:
    """Import a pair of exact core numerator and denominator polynomials.

    Parameters:
        ndim: The retained number of array axes, including unused axes.
        C: The inferred coefficient field.

    Args:
        value: Mutable core numerator and denominator polynomials.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.
        max_steps: The complete call's work allowance.

    Returns:
        An independent, normalized fraction with cancellation evidence.

    Raises:
        PolyError: DOMAIN for invalid inputs or a zero denominator, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return from_core[ndim](value, budget, variables=variables)


from ._numeric_fraction import NumericFunction, Enclosure as FractionEnclosure
from . import _numeric_fraction as _nf, _numeric as _num, evaluation as _ev
from ._conversion import to_sparse as _to_sparse
from polynomials import NumericDomain
from apn_mojo import BallContext, ExactComplex
comptime _NumericDomain[Context: AnyType, A: AnyType, B: AnyType] = _num._Domain[Context, ExactComplex if _num._IsComplex[A] or _num._IsComplex[B] else A]
comptime _PointScalar[X: AnyType] = X
comptime _EnclosureDomain[C: AnyType, X: AnyType] = _ev._Domain[C, _PointScalar[X]]


def asarray[A: Iterable, B: Iterable, Context: Copyable & Deinitable](numerator: A, denominator: B, mut budget: Budget, *, context: Context) raises PolyError -> NumericFunction[Int, _NumericDomain[Context, _Element[A], _Element[B]]]:
    """Construct a numerical fraction without approximate factor cancellation.

    Parameters:
        A: Inferred from the input types.
        B: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        numerator: Ascending numerator coefficients.
        denominator: Ascending denominator coefficients.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        The computed value, retaining the original polynomial axes.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    comptime D = _NumericDomain[Context, _Element[A], _Element[B]]
    var c = rebind[D.Context](context).copy()
    return _nf._make(_to_sparse(_num._dense[D](numerator, c, budget), budget, context=c), _to_sparse(_num._dense[D](denominator, c, budget), budget, context=c), budget)

def asarray[A: Iterable, B: Iterable, Context: Copyable & Deinitable](numerator: A, denominator: B, *, context: Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[Int, _NumericDomain[Context, _Element[A], _Element[B]]]:
    """Construct a numerical fraction without approximate factor cancellation.

    Parameters:
        A: Inferred from the input types.
        B: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        numerator: Ascending numerator coefficients.
        denominator: Ascending denominator coefficients.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A numerical fraction in the context-selected coefficient domain.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return asarray(numerator, denominator, budget, context=context)

def asarray[P: Copyable & Deinitable, k: Int, Context: Copyable & Deinitable](value: Function[P, k], mut budget: Budget, *, context: Context) raises PolyError -> NumericFunction[P, _num._Domain[Context, _Scalar[k]]]:
    """Convert an exact fraction using an explicit numerical context.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        value: An exact characteristic-zero fraction.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A numerical fraction in the context-selected coefficient domain.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert k != 3, "finite-field fractions cannot convert to numerical values"
    comptime D = _num._Domain[Context, _Scalar[k]]
    return _nf._make(_num._sparse[D, sparse.Array[P, k], P, _Scalar[k]](numerator(value, budget), rebind[D.Context](context), budget), _num._sparse[D, sparse.Array[P, k], P, _Scalar[k]](denominator(value, budget), rebind[D.Context](context), budget), budget)

def asarray[P: Copyable & Deinitable, k: Int, Context: Copyable & Deinitable](value: Function[P, k], *, context: Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, _num._Domain[Context, _Scalar[k]]]:
    """Convert an exact fraction using an explicit numerical context.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        value: An exact characteristic-zero fraction.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A numerical fraction in the context-selected coefficient domain.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return asarray(value, budget, context=context)

def from_terms[A: Iterable, B: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, E: Copyable & Deinitable, Context: Copyable & Deinitable](numerator: A, denominator: B, mut budget: Budget, *, context: Context) raises PolyError -> NumericFunction[P, _NumericDomain[Context, C, E]] where _Element[A] == Tuple[P, C] and _Element[B] == Tuple[P, E]:
    """Construct a numerical fraction from sparse terms.

    Parameters:
        A: Inferred from the input types.
        B: Inferred from the input types.
        P: Inferred from the input types.
        C: Inferred from the input types.
        E: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        numerator: Sparse numerator terms.
        denominator: Sparse denominator terms.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A numerical fraction retaining every polynomial axis.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    comptime D = _NumericDomain[Context, C, E]
    var c = rebind[D.Context](context).copy()
    return _nf._make(_num._sparse[D, A, P, C](numerator, c, budget), _num._sparse[D, B, P, E](denominator, c, budget), budget)

def from_terms[A: Iterable, B: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable, E: Copyable & Deinitable, Context: Copyable & Deinitable](numerator: A, denominator: B, *, context: Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, _NumericDomain[Context, C, E]] where _Element[A] == Tuple[P, C] and _Element[B] == Tuple[P, E]:
    """Construct a numerical fraction from sparse terms.

    Parameters:
        A: Inferred from the input types.
        B: Inferred from the input types.
        P: Inferred from the input types.
        C: Inferred from the input types.
        E: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        numerator: Sparse numerator terms.
        denominator: Sparse denominator terms.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A numerical fraction retaining every polynomial axis.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return from_terms[A, B, P, C, E, Context](numerator, denominator, budget, context=context)

def numerator[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], mut budget: Budget) raises PolyError -> sparse.NumericArray[P, D]:
    """Return the stored numerical numerator.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        budget: The shared work allowance.

    Returns:
        The numerical numerator as an immutable sparse array.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    budget.charge(1)
    return value._numerator.copy()

def numerator[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> sparse.NumericArray[P, D]:
    """Return the stored numerical numerator.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        max_steps: The complete call's work allowance.

    Returns:
        The numerical numerator as an immutable sparse array.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return numerator(value, budget)

def denominator[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], mut budget: Budget) raises PolyError -> sparse.NumericArray[P, D]:
    """Return the stored numerical denominator.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        budget: The shared work allowance.

    Returns:
        The numerical denominator as an immutable sparse array.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    budget.charge(1)
    return value._denominator.copy()

def denominator[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> sparse.NumericArray[P, D]:
    """Return the stored numerical denominator.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        max_steps: The complete call's work allowance.

    Returns:
        The numerical denominator as an immutable sparse array.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return denominator(value, budget)

def add[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D]:
    """Add numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        The sum under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var p = sparse.add(sparse.multiply(a._numerator, b._denominator, budget, context=context), sparse.multiply(b._numerator, a._denominator, budget, context=context), budget, context=context)
    var q = sparse.multiply(a._denominator, b._denominator, budget, context=context)
    return _nf._make(p^, q^, budget)

def add[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Add numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        The sum under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return add(a, b, budget, context=context)

def subtract[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D]:
    """Subtract numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        The difference under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var p = sparse.subtract(sparse.multiply(a._numerator, b._denominator, budget, context=context), sparse.multiply(b._numerator, a._denominator, budget, context=context), budget, context=context)
    var q = sparse.multiply(a._denominator, b._denominator, budget, context=context)
    return _nf._make(p^, q^, budget)

def subtract[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Subtract numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        The difference under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return subtract(a, b, budget, context=context)

def multiply[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D]:
    """Multiply numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        The product under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var p = sparse.multiply(a._numerator, b._numerator, budget, context=context)
    var q = sparse.multiply(a._denominator, b._denominator, budget, context=context)
    return _nf._make(p^, q^, budget)

def multiply[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Multiply numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        The product under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return multiply(a, b, budget, context=context)

def divide[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D]:
    """Divide numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        The quotient under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var p = sparse.multiply(a._numerator, b._denominator, budget, context=context)
    var q = sparse.multiply(a._denominator, b._numerator, budget, context=context)
    return _nf._make(p^, q^, budget)

def divide[P: Copyable & Deinitable, D: NumericDomain](a: NumericFunction[P, D], b: NumericFunction[P, D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Divide numerical fractions with the stated rounding schedule.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        a: The first numerical fraction.
        b: The second numerical fraction.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        The quotient under the supplied context, without approximate cancellation.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return divide(a, b, budget, context=context)

def negative[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D]:
    """Return the negative of a numerical fraction.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A fraction with its numerator negated under the supplied context.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    return _nf._make(sparse.negative(value._numerator, budget, context=context), value._denominator.copy(), budget)

def negative[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Return the negative of a numerical fraction.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A fraction with its numerator negated under the supplied context.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return negative(value, budget, context=context)

def reciprocal[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D]:
    """Return the reciprocal of a numerical fraction.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A fraction with numerator and denominator exchanged.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    _num._context[D](context, budget)
    return _nf._make(value._denominator.copy(), value._numerator.copy(), budget)

def reciprocal[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Return the reciprocal of a numerical fraction.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.

    Args:
        value: The numerical fraction.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A fraction with numerator and denominator exchanged.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return reciprocal(value, budget, context=context)

def evaluate[P: Copyable & Deinitable, D: NumericDomain, X: Copyable & Deinitable](value: NumericFunction[P, D], point: X, mut budget: Budget, *, context: D.Context) raises PolyError -> D.Scalar where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Evaluate a numerical fraction; a denominator containing zero raises DOMAIN.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple in the stored numerical domain.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        The rounded or enclosing quotient in the stored numerical domain.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var a = sparse.evaluate(value._numerator, point, budget, context=context)
    var b = sparse.evaluate(value._denominator, point, budget, context=context)
    return _nf._divide[D](a, b, context, budget)

def evaluate[P: Copyable & Deinitable, D: NumericDomain, X: Copyable & Deinitable](value: NumericFunction[P, D], point: X, *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> D.Scalar where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Evaluate a numerical fraction; a denominator containing zero raises DOMAIN.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple in the stored numerical domain.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        The rounded or enclosing quotient in the stored numerical domain.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return evaluate(value, point, budget, context=context)

def enclose[P: Copyable & Deinitable, k: Int, X: Copyable & Deinitable](value: Function[P, k], point: X, mut budget: Budget, *, context: BallContext) raises PolyError -> FractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], X], False, True], _EnclosureDomain[_Scalar[k], X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[_Scalar[k]] >= 0 and k != 3, "enclosure evaluation requires characteristic-zero coefficients"
    var a = sparse.enclose(numerator(value, budget), point, budget, context=context)
    var b = sparse.enclose(denominator(value, budget), point, budget, context=context)
    return rebind_var[FractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], X], False, True], _EnclosureDomain[_Scalar[k], X]]](_nf._enclose(a^, b^, context, budget))

def enclose[P: Copyable & Deinitable, k: Int, X: Copyable & Deinitable](value: Function[P, k], point: X, *, context: BallContext, max_steps: UInt64 = UInt64.MAX) raises PolyError -> FractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], X], False, True], _EnclosureDomain[_Scalar[k], X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return enclose(value, point, budget, context=context)

def enclose[P: Copyable & Deinitable, D: NumericDomain, X: Copyable & Deinitable](value: NumericFunction[P, D], point: X, mut budget: Budget, *, context: BallContext) raises PolyError -> FractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, X], True, True], _EnclosureDomain[D.Scalar, X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var a = sparse.enclose(numerator(value, budget), point, budget, context=context)
    var b = sparse.enclose(denominator(value, budget), point, budget, context=context)
    return _nf._enclose(a^, b^, context, budget)

def enclose[P: Copyable & Deinitable, D: NumericDomain, X: Copyable & Deinitable](value: NumericFunction[P, D], point: X, *, context: BallContext, max_steps: UInt64 = UInt64.MAX) raises PolyError -> FractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, X], True, True], _EnclosureDomain[D.Scalar, X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return enclose(value, point, budget, context=context)

def check[C: AnyType, D: NumericDomain, numerical: Bool](value: FractionEnclosure[_ev.Enclosure[C, D, numerical, True], D], mut budget: Budget) raises PolyError -> Bool:
    """Verify both polynomial evaluations and the exact interval division bound.

    Parameters:
        C: Inferred from the input types.
        D: Inferred from the input types.
        numerical: Inferred from the input types.

    Args:
        value: The retained quotient enclosure.
        budget: The shared work allowance.

    Returns:
        Whether the retained evidence is valid.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    return _nf._check(value, budget)

def check[C: AnyType, D: NumericDomain, numerical: Bool](value: FractionEnclosure[_ev.Enclosure[C, D, numerical, True], D], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify both polynomial evaluations and the exact interval division bound.

    Parameters:
        C: Inferred from the input types.
        D: Inferred from the input types.
        numerical: Inferred from the input types.

    Args:
        value: The retained quotient enclosure.
        max_steps: The complete call's work allowance.

    Returns:
        Whether the retained evidence is valid.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(value, budget)


def evaluate[P: Copyable & Deinitable, D: NumericDomain, *Ts: Copyable & Deinitable](value: NumericFunction[P, D], point: Tuple[*Ts], mut budget: Budget, *, context: D.Context) raises PolyError -> D.Scalar:
    """Evaluate a numerical fraction; a denominator containing zero raises DOMAIN.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple in the stored numerical domain.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        The rounded or enclosing quotient in the stored numerical domain.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var a = sparse.evaluate(value._numerator, point, budget, context=context)
    var b = sparse.evaluate(value._denominator, point, budget, context=context)
    return _nf._divide[D](a, b, context, budget)

def evaluate[P: Copyable & Deinitable, D: NumericDomain, *Ts: Copyable & Deinitable](value: NumericFunction[P, D], point: Tuple[*Ts], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> D.Scalar:
    """Evaluate a numerical fraction; a denominator containing zero raises DOMAIN.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple in the stored numerical domain.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        The rounded or enclosing quotient in the stored numerical domain.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return evaluate(value, point, budget, context=context)

def enclose[P: Copyable & Deinitable, k: Int, *Ts: Copyable & Deinitable](value: Function[P, k], point: Tuple[*Ts], mut budget: Budget, *, context: BallContext) raises PolyError -> FractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]], False, True], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]]]:
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[_Scalar[k]] >= 0 and k != 3, "enclosure evaluation requires characteristic-zero coefficients"
    var a = sparse.enclose(numerator(value, budget), point, budget, context=context)
    var b = sparse.enclose(denominator(value, budget), point, budget, context=context)
    return rebind_var[FractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]], False, True], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]]]](_nf._enclose(a^, b^, context, budget))

def enclose[P: Copyable & Deinitable, k: Int, *Ts: Copyable & Deinitable](value: Function[P, k], point: Tuple[*Ts], *, context: BallContext, max_steps: UInt64 = UInt64.MAX) raises PolyError -> FractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]], False, True], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]]]:
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return enclose(value, point, budget, context=context)

def enclose[P: Copyable & Deinitable, D: NumericDomain, *Ts: Copyable & Deinitable](value: NumericFunction[P, D], point: Tuple[*Ts], mut budget: Budget, *, context: BallContext) raises PolyError -> FractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]], True, True], _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]]]:
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        context: The explicit output precision and rounding context.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var a = sparse.enclose(numerator(value, budget), point, budget, context=context)
    var b = sparse.enclose(denominator(value, budget), point, budget, context=context)
    return _nf._enclose(a^, b^, context, budget)

def enclose[P: Copyable & Deinitable, D: NumericDomain, *Ts: Copyable & Deinitable](value: NumericFunction[P, D], point: Tuple[*Ts], *, context: BallContext, max_steps: UInt64 = UInt64.MAX) raises PolyError -> FractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]], True, True], _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]]]:
    """Evaluate a fraction with independently checked numerator, denominator, and quotient bounds.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        context: The explicit output precision and rounding context.
        max_steps: The complete call's work allowance.

    Returns:
        A quotient enclosure retaining both polynomial evaluations and their evidence.

    Raises:
        PolyError: DOMAIN for incompatible input or a denominator containing zero, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return enclose(value, point, budget, context=context)


from . import _common

def derivative[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], mut budget: Budget, *, context: D.Context, axis: Int = 0) raises PolyError -> NumericFunction[P, D]:
    """Differentiate a numerical fraction using the quotient rule.

    Every arithmetic step uses context; no approximate cancellation is attempted.
    DOMAIN reports an uncertain zero denominator; LIMIT/CAPACITY preserve resource failures.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        value: The input numerical fraction.
        budget: The shared allowance for conversion, arithmetic, and result storage.
        context: The explicit arithmetic precision and rounding or enclosure settings.
        axis: The zero-based polynomial variable axis.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var p = sparse.subtract(sparse.multiply(sparse.derivative(value._numerator, budget, context=context, axis=axis), value._denominator, budget, context=context), sparse.multiply(value._numerator, sparse.derivative(value._denominator, budget, context=context, axis=axis), budget, context=context), budget, context=context)
    var q = sparse.multiply(value._denominator, value._denominator, budget, context=context)
    return _nf._make(p^, q^, budget)

def derivative[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], *, context: D.Context, axis: Int = 0, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Differentiate a numerical fraction using the quotient rule.

    Every arithmetic step uses context; no approximate cancellation is attempted.
    DOMAIN reports an uncertain zero denominator; LIMIT/CAPACITY preserve resource failures.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        value: The input numerical fraction.
        context: The explicit arithmetic precision and rounding or enclosure settings.
        axis: The zero-based polynomial variable axis.
        max_steps: The complete call’s work allowance.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return derivative(value, budget, context=context, axis=axis)

def power[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], exponent: Int, mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D]:
    """Raise a numerical fraction to an integer power.

    Every arithmetic step uses context; no approximate cancellation is attempted.
    DOMAIN reports an uncertain zero denominator; LIMIT/CAPACITY preserve resource failures.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        value: The input numerical fraction.
        exponent: The integer power; negative values swap numerator and denominator.
        budget: The shared allowance for conversion, arithmetic, and result storage.
        context: The explicit arithmetic precision and rounding or enclosure settings.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    if exponent == Int.MIN:
        raise PolyError(PolyErrorKind.CAPACITY, "fraction exponent is too large")
    var n = -exponent if exponent < 0 else exponent
    var p = sparse.power(value._numerator if exponent >= 0 else value._denominator, n, budget, context=context)
    var q = sparse.power(value._denominator if exponent >= 0 else value._numerator, n, budget, context=context)
    return _nf._make(p^, q^, budget)

def power[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], exponent: Int, *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Raise a numerical fraction to an integer power.

    Every arithmetic step uses context; no approximate cancellation is attempted.
    DOMAIN reports an uncertain zero denominator; LIMIT/CAPACITY preserve resource failures.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        value: The input numerical fraction.
        exponent: The integer power; negative values swap numerator and denominator.
        context: The explicit arithmetic precision and rounding or enclosure settings.
        max_steps: The complete call’s work allowance.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return power(value, exponent, budget, context=context)

def substitute[P: Copyable & Deinitable, D: NumericDomain, I: Iterable](value: NumericFunction[P, D], assignments: I, mut budget: Budget, *, context: D.Context) raises PolyError -> NumericFunction[P, D] where _Element[I] == Tuple[Int, sparse.NumericArray[P, D]]:
    """Substitute numerical polynomials simultaneously in a fraction.

    Every arithmetic step uses context; no approximate cancellation is attempted.
    DOMAIN reports an uncertain zero denominator; LIMIT/CAPACITY preserve resource failures.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.
        I: The replacement iterable type.

    Args:
        value: The input numerical fraction.
        assignments: Distinct axis and numerical-polynomial pairs, applied simultaneously.
        budget: The shared allowance for conversion, arithmetic, and result storage.
        context: The explicit arithmetic precision and rounding or enclosure settings.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    # Materialize once so generators are not consumed twice.
    var items = List[Tuple[Int, sparse.NumericArray[P, D]]]()
    var cursor = iter(assignments)
    while True:
        budget.charge(1)
        var item: Tuple[Int, sparse.NumericArray[P, D]]
        try:
            item = rebind_var[Tuple[Int, sparse.NumericArray[P, D]]](next(cursor))
        except StopIteration:
            break
        _grow(items, budget)
        items.append(item^)
    return _nf._make(sparse.substitute(value._numerator, items, budget, context=context), sparse.substitute(value._denominator, items, budget, context=context), budget)

def substitute[P: Copyable & Deinitable, D: NumericDomain, I: Iterable](value: NumericFunction[P, D], assignments: I, *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D] where _Element[I] == Tuple[Int, sparse.NumericArray[P, D]]:
    """Substitute numerical polynomials simultaneously in a fraction.

    Every arithmetic step uses context; no approximate cancellation is attempted.
    DOMAIN reports an uncertain zero denominator; LIMIT/CAPACITY preserve resource failures.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.
        I: The replacement iterable type.

    Args:
        value: The input numerical fraction.
        assignments: Distinct axis and numerical-polynomial pairs, applied simultaneously.
        context: The explicit arithmetic precision and rounding or enclosure settings.
        max_steps: The complete call’s work allowance.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return substitute(value, assignments, budget, context=context)


def _numeric_constant[P: Copyable & Deinitable, D: NumericDomain, T: Copyable & Deinitable](value: T, context: D.Context, mut budget: Budget) raises PolyError -> sparse.NumericArray[P, D]:
    from polynomials import NumericSparsePoly, NumericTerm, Monomial
    var scalar = _num._value[D](value, context, budget)
    _allocation[NumericTerm[D]](1, budget)
    var terms = List[NumericTerm[D]]()
    terms.append(NumericTerm[D](Monomial(), scalar^))
    return sparse.NumericArray[P, D](NumericSparsePoly[D].from_terms(terms^, sparse._order[P](budget), context, budget), budget)


def _numeric_compose_polynomial[P: Copyable & Deinitable, D: NumericDomain](p: sparse.NumericArray[P, D], inner: NumericFunction[P, D], axis: Int, context: D.Context, mut budget: Budget) raises PolyError -> NumericFunction[P, D]:
    var one = _numeric_constant[P, D](1, context, budget)
    var result = _nf._make(_numeric_constant[P, D](0, context, budget), one.copy(), budget)
    for term in p._core[].terms:
        var exponent = 0
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] != Int64(axis):
                raise PolyError(PolyErrorKind.DOMAIN, "composition requires univariate fractions on the selected axis")
            exponent = Int(pair[1])
        var coefficient = _nf._make(_numeric_constant[P, D](term.coefficient, context, budget), one.copy(), budget)
        var contribution = multiply(coefficient, power(inner, exponent, budget, context=context), budget, context=context)
        result = add(result, contribution, budget, context=context)
    return result^


def compose[P: Copyable & Deinitable, D: NumericDomain](outer: NumericFunction[P, D], inner: NumericFunction[P, D], mut budget: Budget, *, context: D.Context, axis: Int = 0) raises PolyError -> NumericFunction[P, D]:
    """Compose univariate numerical fractions on the selected axis.

    Uses the supplied arithmetic context throughout, without approximate GCD cancellation.
    DOMAIN reports another occurring axis or an uncertain zero denominator.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        outer: The outer univariate fraction.
        inner: The univariate replacement fraction.
        budget: The shared allowance for conversion, arithmetic, and result storage.
        context: The explicit arithmetic precision and rounding or enclosure settings.
        axis: The zero-based polynomial variable axis.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    sparse._axis[P](axis)
    for term in inner._numerator._core[].terms:
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] != Int64(axis):
                raise PolyError(PolyErrorKind.DOMAIN, "composition requires univariate fractions on the selected axis")
    for term in inner._denominator._core[].terms:
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] != Int64(axis):
                raise PolyError(PolyErrorKind.DOMAIN, "composition requires univariate fractions on the selected axis")
    var a = _numeric_compose_polynomial(outer._numerator, inner, axis, context, budget)
    var b = _numeric_compose_polynomial(outer._denominator, inner, axis, context, budget)
    return divide(a, b, budget, context=context)


def compose[P: Copyable & Deinitable, D: NumericDomain](outer: NumericFunction[P, D], inner: NumericFunction[P, D], *, context: D.Context, axis: Int = 0, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, D]:
    """Compose numerical fractions with explicit precision and a work limit.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        outer: The outer univariate fraction.
        inner: The univariate replacement fraction.
        context: The explicit arithmetic precision and rounding or enclosure settings.
        axis: The zero-based polynomial variable axis.
        max_steps: The complete call’s work allowance.

    Returns:
        A new numerical fraction with the same polynomial axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return compose(outer, inner, budget, context=context, axis=axis)


def to_core[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], mut budget: Budget, *, variables: List[Int64] = List[Int64]()) raises PolyError -> Tuple[_num.NumericSparsePoly[D], _num.NumericSparsePoly[D]]:
    """Copy numerical numerator and denominator into independent core storage.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        value: The input numerical fraction.
        budget: The shared allowance for conversion, arithmetic, and result storage.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.

    Returns:
        An independent pair of mutable numerical core polynomials.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    return sparse.to_core(value._numerator, budget, variables=variables), sparse.to_core(value._denominator, budget, variables=variables)


def to_core[P: Copyable & Deinitable, D: NumericDomain](value: NumericFunction[P, D], *, max_steps: UInt64 = UInt64.MAX, variables: List[Int64] = List[Int64]()) raises PolyError -> Tuple[_num.NumericSparsePoly[D], _num.NumericSparsePoly[D]]:
    """Copy numerical numerator and denominator with a work limit.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.

    Args:
        value: The input numerical fraction.
        max_steps: The complete call’s work allowance.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.

    Returns:
        An independent pair of mutable numerical core polynomials.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return to_core(value, budget, variables=variables)


def from_core[ndim: Int, D: NumericDomain](value: Tuple[_num.NumericSparsePoly[D], _num.NumericSparsePoly[D]], mut budget: Budget, *, variables: List[Int64] = List[Int64]()) raises PolyError -> NumericFunction[_Shape[ndim], D]:
    """Validate numerical core polynomials and retain the requested axis count.

    Parameters:
        ndim: The number of retained polynomial axes.
        D: The stored numerical coefficient domain.

    Args:
        value: The numerical sparse numerator and denominator to validate and copy.
        budget: The shared allowance for conversion, arithmetic, and result storage.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.

    Returns:
        A numerical fraction retaining exactly ndim axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    return _nf._make(sparse.from_core[D, ndim=ndim](value[0], budget, variables=variables), sparse.from_core[D, ndim=ndim](value[1], budget, variables=variables), budget)


def from_core[ndim: Int, D: NumericDomain](value: Tuple[_num.NumericSparsePoly[D], _num.NumericSparsePoly[D]], *, max_steps: UInt64 = UInt64.MAX, variables: List[Int64] = List[Int64]()) raises PolyError -> NumericFunction[_Shape[ndim], D]:
    """Import numerical core polynomials with an explicit shape and work limit.

    Parameters:
        ndim: The number of retained polynomial axes.
        D: The stored numerical coefficient domain.

    Args:
        value: The numerical sparse numerator and denominator to validate and copy.
        max_steps: The complete call’s work allowance.
        variables: Core variable IDs in array-axis order; empty uses consecutive IDs.

    Returns:
        A numerical fraction retaining exactly ndim axes.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return from_core[ndim, D](value, budget, variables=variables)


def asarray[P: Copyable & Deinitable, D: NumericDomain, Context: Copyable & Deinitable](value: NumericFunction[P, D], mut budget: Budget, *, context: Context) raises PolyError -> NumericFunction[P, _num._Domain[Context, D.Scalar]]:
    """Convert a numerical fraction into the domain selected by context.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.
        Context: The explicit target context type.

    Args:
        value: The input numerical fraction.
        budget: The shared allowance for conversion, arithmetic, and result storage.
        context: The explicit arithmetic precision and rounding or enclosure settings.

    Returns:
        A numerical fraction in the context-selected coefficient domain.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    comptime Target = _num._Domain[Context, D.Scalar]
    var c = rebind[Target.Context](context).copy()
    return _nf._make(_num._sparse[Target, sparse.NumericArray[P, D], P, D.Scalar](value._numerator, c, budget), _num._sparse[Target, sparse.NumericArray[P, D], P, D.Scalar](value._denominator, c, budget), budget)


def asarray[P: Copyable & Deinitable, D: NumericDomain, Context: Copyable & Deinitable](value: NumericFunction[P, D], *, context: Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericFunction[P, _num._Domain[Context, D.Scalar]]:
    """Convert a numerical fraction using an explicit context and work limit.

    Parameters:
        P: The inferred polynomial power shape.
        D: The stored numerical coefficient domain.
        Context: The explicit target context type.

    Args:
        value: The input numerical fraction.
        context: The explicit arithmetic precision and rounding or enclosure settings.
        max_steps: The complete call’s work allowance.

    Returns:
        A numerical fraction in the context-selected coefficient domain.

    Raises:
        PolyError: DOMAIN for invalid axes, incompatible inputs, or an uncertain zero denominator; LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return asarray(value, budget, context=context)

from ._numeric_fraction import AdaptiveEnclosure as AdaptiveFractionEnclosure
from polynomials import NumericOptions


def adaptive_evaluate[P: Copyable & Deinitable, k: Int, X: Copyable & Deinitable](value: Function[P, k], point: X, mut budget: Budget, *, options: NumericOptions) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], X], False, True], _EnclosureDomain[_Scalar[k], X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        options: Initial and maximum precision, plus absolute or relative width targets.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[_Scalar[k]] >= 0 and k != 3, "enclosure evaluation requires characteristic-zero coefficients"
    options._validate()
    var p = numerator(value, budget)
    var q = denominator(value, budget)
    var bits = options.initial_precision
    while True:
        budget.charge(1)
        var context = _ev._context(bits, options.max_precision)
        var a = sparse.enclose(p, point, budget, context=context)
        var b = sparse.enclose(q, point, budget, context=context)
        var result = _nf._adaptive_step(a^, b^, options, context, budget)
        if result.accurate() or bits == options.max_precision:
            return rebind_var[AdaptiveFractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], X], False, True], _EnclosureDomain[_Scalar[k], X]]](result^)
        bits = options.max_precision if bits > options.max_precision // 2 else bits * 2


def adaptive_evaluate[P: Copyable & Deinitable, k: Int, X: Copyable & Deinitable](value: Function[P, k], point: X, *, options: NumericOptions, max_steps: UInt64 = UInt64.MAX) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], X], False, True], _EnclosureDomain[_Scalar[k], X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        options: Initial and maximum precision, plus absolute or relative width targets.
        max_steps: The complete call's work allowance.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return adaptive_evaluate(value, point, budget, options=options)


def adaptive_evaluate[P: Copyable & Deinitable, D: NumericDomain, X: Copyable & Deinitable](value: NumericFunction[P, D], point: X, mut budget: Budget, *, options: NumericOptions) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, X], True, True], _EnclosureDomain[D.Scalar, X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        options: Initial and maximum precision, plus absolute or relative width targets.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    options._validate()
    var p = numerator(value, budget)
    var q = denominator(value, budget)
    var bits = options.initial_precision
    while True:
        budget.charge(1)
        var context = _ev._context(bits, options.max_precision)
        var a = sparse.enclose(p, point, budget, context=context)
        var b = sparse.enclose(q, point, budget, context=context)
        var result = _nf._adaptive_step(a^, b^, options, context, budget)
        if result.accurate() or bits == options.max_precision:
            return rebind_var[AdaptiveFractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, X], True, True], _EnclosureDomain[D.Scalar, X]]](result^)
        bits = options.max_precision if bits > options.max_precision // 2 else bits * 2


def adaptive_evaluate[P: Copyable & Deinitable, D: NumericDomain, X: Copyable & Deinitable](value: NumericFunction[P, D], point: X, *, options: NumericOptions, max_steps: UInt64 = UInt64.MAX) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, X], True, True], _EnclosureDomain[D.Scalar, X]] where (_Kind[X] >= 0 or _num._IsNumeric[X]):
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        X: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        options: Initial and maximum precision, plus absolute or relative width targets.
        max_steps: The complete call's work allowance.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return adaptive_evaluate(value, point, budget, options=options)


def adaptive_evaluate[P: Copyable & Deinitable, k: Int, *Ts: Copyable & Deinitable](value: Function[P, k], point: Tuple[*Ts], mut budget: Budget, *, options: NumericOptions) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]], False, True], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]]]:
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        options: Initial and maximum precision, plus absolute or relative width targets.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[_Scalar[k]] >= 0 and k != 3, "enclosure evaluation requires characteristic-zero coefficients"
    options._validate()
    var p = numerator(value, budget)
    var q = denominator(value, budget)
    var bits = options.initial_precision
    while True:
        budget.charge(1)
        var context = _ev._context(bits, options.max_precision)
        var a = sparse.enclose(p, point, budget, context=context)
        var b = sparse.enclose(q, point, budget, context=context)
        var result = _nf._adaptive_step(a^, b^, options, context, budget)
        if result.accurate() or bits == options.max_precision:
            return rebind_var[AdaptiveFractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]], False, True], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]]]](result^)
        bits = options.max_precision if bits > options.max_precision // 2 else bits * 2


def adaptive_evaluate[P: Copyable & Deinitable, k: Int, *Ts: Copyable & Deinitable](value: Function[P, k], point: Tuple[*Ts], *, options: NumericOptions, max_steps: UInt64 = UInt64.MAX) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[_Ring[k], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]], False, True], _EnclosureDomain[_Scalar[k], _num._TupleScalar[Tuple[*Ts]]]]:
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        k: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The exact fraction.
        point: A scalar or coordinate tuple.
        options: Initial and maximum precision, plus absolute or relative width targets.
        max_steps: The complete call's work allowance.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return adaptive_evaluate(value, point, budget, options=options)


def adaptive_evaluate[P: Copyable & Deinitable, D: NumericDomain, *Ts: Copyable & Deinitable](value: NumericFunction[P, D], point: Tuple[*Ts], mut budget: Budget, *, options: NumericOptions) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]], True, True], _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]]]:
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        budget: The shared work allowance.
        options: Initial and maximum precision, plus absolute or relative width targets.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    options._validate()
    var p = numerator(value, budget)
    var q = denominator(value, budget)
    var bits = options.initial_precision
    while True:
        budget.charge(1)
        var context = _ev._context(bits, options.max_precision)
        var a = sparse.enclose(p, point, budget, context=context)
        var b = sparse.enclose(q, point, budget, context=context)
        var result = _nf._adaptive_step(a^, b^, options, context, budget)
        if result.accurate() or bits == options.max_precision:
            return rebind_var[AdaptiveFractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]], True, True], _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]]]](result^)
        bits = options.max_precision if bits > options.max_precision // 2 else bits * 2


def adaptive_evaluate[P: Copyable & Deinitable, D: NumericDomain, *Ts: Copyable & Deinitable](value: NumericFunction[P, D], point: Tuple[*Ts], *, options: NumericOptions, max_steps: UInt64 = UInt64.MAX) raises PolyError -> AdaptiveFractionEnclosure[_ev.Enclosure[D, _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]], True, True], _EnclosureDomain[D.Scalar, _num._TupleScalar[Tuple[*Ts]]]]:
    """Refine a quotient until its requested width is certified or precision is exhausted.

    Parameters:
        P: Inferred from the input types.
        D: Inferred from the input types.
        Ts: Inferred from the input types.

    Args:
        value: The numerical fraction.
        point: A scalar or coordinate tuple.
        options: Initial and maximum precision, plus absolute or relative width targets.
        max_steps: The complete call's work allowance.

    Returns:
        Checked numerator and denominator evaluations, an accuracy outcome, and
        an optional finite quotient. has_value() is False when denominator
        separation remains unresolved; this does not by itself prove a pole.

    Raises:
        PolyError: DOMAIN for incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return adaptive_evaluate(value, point, budget, options=options)


def check[E: _nf._QuotientInput, D: NumericDomain](value: AdaptiveFractionEnclosure[E, D], mut budget: Budget) raises PolyError -> Bool:
    """Verify the quotient bounds and any retained accuracy outcome.

    Parameters:
        E: The retained polynomial evaluation type.
        D: The real or complex ball domain.

    Args:
        value: The adaptive quotient and its evidence.
        budget: The shared checking allowance.

    Returns:
        Whether the evidence supports the outcome, including an honest precision ceiling.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    return _nf._check_adaptive(value, budget)


def check[E: _nf._QuotientInput, D: NumericDomain](value: AdaptiveFractionEnclosure[E, D], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the quotient bounds and any retained accuracy outcome.

    Parameters:
        E: The retained polynomial evaluation type.
        D: The real or complex ball domain.

    Args:
        value: The adaptive quotient and its evidence.
        max_steps: The complete checking allowance.

    Returns:
        Whether the evidence supports the outcome, including an honest precision ceiling.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return _nf._check_adaptive(value, budget)
