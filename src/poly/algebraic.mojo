"""Exact real and complex algebraic values with functional arithmetic and evidence."""
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink, BallDomain, QQ
from ._common import _grow
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration, iter, next
from std.builtin.rebind import downcast
from polynomials import SparsePoly, AlgebraicSpecializationCertificate, algebraic_evaluate, check_algebraic_evaluation, algebraic_evaluate_sparse, check_algebraic_evaluate_sparse, algebraic_specialize, check_algebraic_specialization
from ._common import _Element, Values, _check_batch

from apn_mojo import Rational, ExactComplex
from polynomials import Budget, PolyError, PolyErrorKind, NoContext, DensePoly, RealAlgebraic
from polynomials import real_rational, real_root, check_real_algebraic, algebraic_compare, check_algebraic_arithmetic
from polynomials import algebraic_multiply, algebraic_inverse, refine as _refine, RootInterval
from ._common import _Kind, _convert, _allocation, _charge_size


from ._algebraic_value import Real, Complex, _ValueData, _HornerData, _SparseEvaluationData, _plain, _make, _arithmetic, _copy_witness, _same_witness


def asreal[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> Real:
    """Construct an exact real value from a scalar or core root witness.

    Parameters:
        T: The inferred input type.

    Args:
        value: An integer, Rational, existing Real, or mutable core RealAlgebraic witness.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime if T == Real:
        budget.charge(1)
        return rebind[Real](value).copy()
    elif T == RealAlgebraic:
        return _plain(_copy_witness(rebind[RealAlgebraic](value), budget), budget)
    else:
        comptime assert _Kind[T] == 0 or _Kind[T] == 1, "real algebraic values require exact integers or rationals"
        return _plain(real_rational(QQ(_convert[1](value, NoContext(), budget)), budget), budget)


def asreal[T: Copyable & Deinitable](value: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Real:
    """Construct an exact real value from a scalar or core root witness.

    Parameters:
        T: The inferred input type.

    Args:
        value: An integer, Rational, existing Real, or mutable core RealAlgebraic witness.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return asreal(value, budget)


def add[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, mut budget: Budget) raises PolyError -> _Arithmetic[A, B]:
    """Add exact real or complex algebraic values with retained component evidence.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        b: The right exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _ComplexInput[A] or _ComplexInput[B]:
        return rebind_var[_Arithmetic[A, B]](ascomplex(a, budget).add(ascomplex(b, budget), budget))
    else:
        return rebind_var[_Arithmetic[A, B]](_arithmetic(asreal(a, budget), Optional(asreal(b, budget)), 0, budget))


def add[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Arithmetic[A, B]:
    """Add exact real or complex algebraic values with retained component evidence.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        b: The right exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return add(a, b, budget)


def subtract[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, mut budget: Budget) raises PolyError -> _Arithmetic[A, B]:
    """Subtract exact real or complex algebraic values with retained component evidence.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        b: The right exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _ComplexInput[A] or _ComplexInput[B]:
        return rebind_var[_Arithmetic[A, B]](ascomplex(a, budget).sub(ascomplex(b, budget), budget))
    else:
        return rebind_var[_Arithmetic[A, B]](_arithmetic(asreal(a, budget), Optional(asreal(b, budget)), 1, budget))


def subtract[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Arithmetic[A, B]:
    """Subtract exact real or complex algebraic values with retained component evidence.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        b: The right exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return subtract(a, b, budget)


def multiply[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, mut budget: Budget) raises PolyError -> _Arithmetic[A, B]:
    """Multiply exact real or complex algebraic values with retained component evidence.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        b: The right exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _ComplexInput[A] or _ComplexInput[B]:
        return rebind_var[_Arithmetic[A, B]](ascomplex(a, budget).mul(ascomplex(b, budget), budget))
    else:
        return rebind_var[_Arithmetic[A, B]](_arithmetic(asreal(a, budget), Optional(asreal(b, budget)), 2, budget))


def multiply[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Arithmetic[A, B]:
    """Multiply exact real or complex algebraic values with retained component evidence.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        b: The right exact integer, rational, Real, ExactComplex, or algebraic.Complex.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return multiply(a, b, budget)


def negative[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> _Arithmetic[T, T]:
    """Negate an exact real or complex algebraic value.

    Parameters:
        T: The inferred input type.

    Args:
        value: An exact real, integer, or rational; inversion requires a nonzero value.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _ComplexInput[T]:
        return rebind_var[_Arithmetic[T, T]](ascomplex(value, budget).neg(budget))
    else:
        return rebind_var[_Arithmetic[T, T]](_arithmetic(asreal(value, budget), Optional[Real](), 3, budget))


def negative[T: Copyable & Deinitable](value: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Arithmetic[T, T]:
    """Negate an exact real or complex algebraic value.

    Parameters:
        T: The inferred input type.

    Args:
        value: An exact real, integer, or rational; inversion requires a nonzero value.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return negative(value, budget)


def reciprocal[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> _Arithmetic[T, T]:
    """Invert an exact real or complex algebraic value.

    Parameters:
        T: The inferred input type.

    Args:
        value: An exact real, integer, or rational; inversion requires a nonzero value.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _ComplexInput[T]:
        return rebind_var[_Arithmetic[T, T]](ascomplex(value, budget).inverse(budget))
    else:
        return rebind_var[_Arithmetic[T, T]](_arithmetic(asreal(value, budget), Optional[Real](), 4, budget))


def reciprocal[T: Copyable & Deinitable](value: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Arithmetic[T, T]:
    """Invert an exact real or complex algebraic value.

    Parameters:
        T: The inferred input type.

    Args:
        value: An exact real, integer, or rational; inversion requires a nonzero value.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return reciprocal(value, budget)


def divide[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, mut budget: Budget) raises PolyError -> _Arithmetic[A, B]:
    """Divide exact real or complex algebraic values.

    Real results retain both inverse and multiplication certificates. Complex
    results retain the latest arithmetic evidence for each real component.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The numerator.
        b: The nonzero denominator.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime if _ComplexInput[A] or _ComplexInput[B]:
        return rebind_var[_Arithmetic[A, B]](ascomplex(a, budget).mul(ascomplex(b, budget).inverse(budget), budget))
    else:
        var numerator = asreal(a, budget)
        var denominator = asreal(b, budget)
        var left = _copy_witness(numerator.witness(), budget)
        var divisor = _copy_witness(denominator.witness(), budget)
        var inverse = algebraic_inverse(divisor, budget)
        var right = _copy_witness(inverse.value, budget)
        var proof = algebraic_multiply(left, right, budget)
        var value = _copy_witness(proof.value, budget)
        _allocation[_ValueData](1, budget)
        return rebind_var[_Arithmetic[A, B]](_make(_ValueData(value^, Optional(proof^), Optional(left^), Optional(right^), Optional(divisor^), Optional(inverse^)), budget))


def divide[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Arithmetic[A, B]:
    """Divide exact real or complex algebraic values.

    Real results retain both inverse and multiplication certificates. Complex
    results retain the latest arithmetic evidence for each real component.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The numerator.
        b: The nonzero denominator.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return divide(a, b, budget)


def compare[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, mut budget: Budget) raises PolyError -> Int:
    """Compare two exact real values.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact value.
        b: The right exact value.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Minus one, zero, or one according to the exact ordering.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var left = asreal(a, budget)
    var right = asreal(b, budget)
    return algebraic_compare(left.witness(), right.witness(), budget).relation


def compare[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Int:
    """Compare two exact real values.

    Parameters:
        A: The inferred input type.
        B: The inferred input type.

    Args:
        a: The left exact value.
        b: The right exact value.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Minus one, zero, or one according to the exact ordering.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return compare(a, b, budget)


def check(value: Real, mut budget: Budget) raises PolyError -> Bool:
    """Verify the value and its retained arithmetic step, if present.

    Args:
        value: The exact real value to verify.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Whether the retained evidence proves the value or arithmetic step.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    if value._horner:
        ref data = value._horner.value()[]
        if not _same_witness(value.witness(), data.proof.value, budget) or not check_algebraic_evaluation(data.coefficients, data.point, data.proof, budget):
            return False
    if value._sparse_evaluation:
        ref data = value._sparse_evaluation.value()[]
        if not _same_witness(value.witness(), data.proof.value, budget) or not check_algebraic_evaluate_sparse(data.polynomial, data.point, data.proof, budget):
            return False
    if value._data[].division:
        if not value._data[].divisor or not value._data[].right:
            return False
        if not check_algebraic_arithmetic(value._data[].divisor.value(), Optional[RealAlgebraic](), value._data[].division.value(), budget) or not _same_witness(value._data[].right.value(), value._data[].division.value().value, budget):
            return False
    if value._data[].arithmetic:
        return check_real_algebraic(value.witness(), budget) and _same_witness(value.witness(), value._data[].arithmetic.value().value, budget) and check_algebraic_arithmetic(value._data[].left.value(), value._data[].right, value._data[].arithmetic.value(), budget)
    return check_real_algebraic(value.witness(), budget)


def check(value: Real, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify the value and its retained arithmetic step, if present.

    Args:
        value: The exact real value to verify.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the retained evidence proves the value or arithmetic step.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(value, budget)


def bounds(value: Real, mut budget: Budget) raises PolyError -> Tuple[Rational, Rational]:
    """Return rational bounds for an exact real value.

    Args:
        value: The exact real value.
        budget: The shared work allowance for preparation and computation.

    Returns:
        Equal endpoints for rationals, or an open isolating interval for an irrational root.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    ref core = value.witness()
    _charge_size(core.low.size(), budget)
    _charge_size(core.high.size(), budget)
    return core.low.value.copy(), core.high.value.copy()


def bounds(value: Real, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Tuple[Rational, Rational]:
    """Return rational bounds for an exact real value.

    Args:
        value: The exact real value.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Equal endpoints for rationals, or an open isolating interval for an irrational root.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return bounds(value, budget)


def refine[T: Copyable & Deinitable](value: Real, width: T, mut budget: Budget) raises PolyError -> Real:
    """Return the same exact root with a narrower rational enclosure.

    The refined result retains a root witness. Its earlier arithmetic derivation remains on the original value.

    Parameters:
        T: The inferred input type.

    Args:
        value: The exact value; the input remains unchanged.
        width: A positive exact bound on the returned interval width.
        budget: The shared work allowance for preparation and computation.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    comptime assert _Kind[T] == 0 or _Kind[T] == 1, "root width must be an exact integer or rational"
    var target = QQ(_convert[1](width, NoContext(), budget))
    if target.value.sign() <= 0:
        raise PolyError(PolyErrorKind.DOMAIN, "root width must be positive")
    ref core = value.witness()
    if core.rational:
        budget.charge(1)
        return value.copy()
    _charge_size(core.low.size(), budget)
    _charge_size(core.high.size(), budget)
    var interval = _refine(core.polynomial, RootInterval(core.low.copy(), core.high.copy()), target, budget)
    return _plain(real_root(core.polynomial, interval.low, interval.high, budget), budget)


def refine[T: Copyable & Deinitable](value: Real, width: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Real:
    """Return the same exact root with a narrower rational enclosure.

    The refined result retains a root witness. Its earlier arithmetic derivation remains on the original value.

    Parameters:
        T: The inferred input type.

    Args:
        value: The exact value; the input remains unchanged.
        width: A positive exact bound on the returned interval width.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The computed result, with its evidence retained.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return refine(value, width, budget)


def approximate[T: Copyable & Deinitable](value: Real, width: T, mut budget: Budget) raises PolyError -> Rational:
    """Return the rational midpoint of an explicitly narrowed enclosure.

    Parameters:
        T: The inferred input type.

    Args:
        value: The exact value to approximate.
        width: A positive interval-width target.
        budget: The shared work allowance for preparation and computation.

    Returns:
        A rational approximation with absolute error less than width/2; exact for rational inputs.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var narrowed = refine(value, width, budget)
    ref core = narrowed.witness()
    return core.low.add(core.high, budget).mul(QQ(Rational(2)).inverse(budget), budget).value


def approximate[T: Copyable & Deinitable](value: Real, width: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Rational:
    """Return the rational midpoint of an explicitly narrowed enclosure.

    Parameters:
        T: The inferred input type.

    Args:
        value: The exact value to approximate.
        width: A positive interval-width target.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A rational approximation with absolute error less than width/2; exact for rational inputs.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable operation, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return approximate(value, width, budget)


def evaluate[I: Iterable, X: Copyable & Deinitable](coefficients: I, point: X, mut budget: Budget) raises PolyError -> Real:
    """Evaluate ascending real coefficients exactly, retaining the complete Horner proof.

    Parameters:
        I: The inferred input iterable type.
        X: The evaluation point or seed type.

    Args:
        coefficients: Ascending real coefficients, including exact algebraic values.
        point: The evaluation scalar or ordered coordinate tuple.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        The evaluated scalar; algebraic results retain the evaluation proof.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    _check_batch(coefficients)
    var values = List[RealAlgebraic]()
    var cursor = iter(coefficients)
    while True:
        budget.charge(1)
        var value: downcast[_Element[I], Copyable & Deinitable]
        try:
            value = rebind_var[downcast[_Element[I], Copyable & Deinitable]](next(cursor))
        except StopIteration:
            break
        var real = asreal(value, budget)
        _grow(values, budget)
        values.append(_copy_witness(real.witness(), budget))
    var x = asreal(point, budget)
    var witness = _copy_witness(x.witness(), budget)
    var proof = algebraic_evaluate(values, witness, budget)
    var result = _plain(_copy_witness(proof.value, budget), budget)
    _allocation[_HornerData](1, budget)
    result._horner = ArcPointer(_HornerData(values^, witness^, proof^))
    return result^


def evaluate[I: Iterable, X: Copyable & Deinitable](coefficients: I, point: X, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Real:
    """Evaluate ascending real coefficients exactly with one work allowance.

    Parameters:
        I: The inferred input iterable type.
        X: The evaluation point or seed type.

    Args:
        coefficients: Ascending real coefficients, including exact algebraic values.
        point: The evaluation scalar or ordered coordinate tuple.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The evaluated scalar; algebraic results retain the evaluation proof.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return evaluate(coefficients, point, budget)


def _evaluate_sparse(var polynomial: SparsePoly[QQ], var point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> Real:
    var proof = algebraic_evaluate_sparse(polynomial, point, budget)
    var result = _plain(_copy_witness(proof.value, budget), budget)
    _allocation[_SparseEvaluationData](1, budget)
    result._sparse_evaluation = ArcPointer(_SparseEvaluationData(polynomial^, point^, proof^))
    return result^


@fieldwise_init
struct _SpecializationData(Copyable, Movable):
    var polynomial: SparsePoly[QQ]
    var point: List[Tuple[Int64, RealAlgebraic]]
    var proof: AlgebraicSpecializationCertificate


@fieldwise_init
struct Specialization[P: Copyable & Deinitable]( _Archive, Copyable, Movable, Iterable, Sized):
    """Sparse algebraic coefficients after a partial assignment; original axes are retained.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.
    """
    var _data: ArcPointer[_SpecializationData]
    var _values: Values[Tuple[Self.P, Real], 1]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _SpecializationIterator[Self.P, ImmOrigin(iterable_origin)]
    """Iterate over remaining powers and exact real coefficients.

    Parameters:
        iterable_mut: Whether the borrowing source origin allows mutation.
        iterable_origin: The lifetime of the borrowed source.

    Returns:
        The borrowing iterator type.
    """

    def __len__(self) -> Int:
        """Return the number of retained coefficient groups, including exact zeros.

        Returns:
            The number of logical entries.
        """
        return len(self._values)

    def __getitem__(self, index: Int) -> Tuple[Self.P, Real]:
        """Copy a power and algebraic value pair at a valid index.

        Args:
            index: A valid zero-based index into the logical entries.

        Returns:
            An independent value copy; polynomial and evidence buffers may remain shared.
        """
        return self._values[index]

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over powers and exact coefficients with shared witnesses.

        Returns:
            An iterator borrowing the immutable source storage.
        """
        return _SpecializationIterator[Self.P, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] AlgebraicSpecializationCertificate:
        """Borrow the immutable proofs for every specialized coefficient.

        Returns:
            An immutable reference to the retained core evidence.
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
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_SpecializationData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)
        writer.line(String(len(self._values)), sink, budget)
        for field736 in self._values:
            _result_io._write_node[Self.P, _ArchiveRing, _ArchiveDomain, _ArchiveProof](field736[0], writer, sink, budget, depth + 1)
            _result_io._write_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](field736[1], writer, sink, budget, depth + 1)
        writer.line(String(self._values.modulus()), sink, budget)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = QQ
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field738 = _result_io._read_node[_SpecializationData, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_SpecializationData](1, budget)
        var field737 = ArcPointer(field738^)
        var field741 = reader.count(budget)
        var field740 = List[Tuple[Self.P, Real]]()
        for _ in range(field741):
            var field743 = _result_io._read_node[Self.P, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            var field744 = _result_io._read_node[Real, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
            var field742 = (field743^, field744^)
            _grow(field740, budget)
            field740.append(field742^)
        var field739 = _common.Values[Tuple[Self.P, Real], 1](field740^, reader.number(UInt64.MAX, budget), budget)
        return Self(field737^, field739^)

@fieldwise_init
struct _SpecializationIterator[P: Copyable & Deinitable, origin: ImmOrigin](Iterator):
    comptime Element = Tuple[Self.P, Real]
    var source: Pointer[Specialization[Self.P], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.Element:
        if self.position == len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position]
        self.position += 1
        return result^


def _specialize[P: Copyable & Deinitable](var polynomial: SparsePoly[QQ], var point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> Specialization[P]:
    from .sparse import _powers
    var proof = algebraic_specialize(polynomial, point, budget)
    _allocation[Tuple[P, Real]](len(proof.monomials), budget)
    var values = List[Tuple[P, Real]](capacity=len(proof.monomials))
    for i in range(len(proof.monomials)):
        budget.charge(1)
        values.append((_powers[P](proof.monomials[i]), _plain(_copy_witness(proof.evaluations[i].value, budget), budget)))
    _allocation[_SpecializationData](1, budget)
    return Specialization[P](ArcPointer(_SpecializationData(polynomial^, point^, proof^)), Values[Tuple[P, Real], 1](values^, 0, budget))


def check[P: Copyable & Deinitable](result: Specialization[P], mut budget: Budget) raises PolyError -> Bool:
    """Check every specialized coefficient and bind the proof to the returned terms.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Whether the evidence binds the specialized coefficients to the original
        polynomial and assignments.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    from .sparse import _monomial
    ref data = result._data[]
    if len(result) != len(data.proof.monomials) or not check_algebraic_specialization(data.polynomial, data.point, data.proof, budget):
        return False
    for i in range(len(result)):
        var term = result[i]
        if _monomial(term[0], budget) != data.proof.monomials[i] or not _same_witness(term[1].witness(), data.proof.evaluations[i].value, budget):
            return False
    return True


def check[P: Copyable & Deinitable](result: Specialization[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify a specialization against its retained input and assignments.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the evidence binds the specialized coefficients to the original
        polynomial and assignments.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


comptime _AnyReal[previous: Bool, T: Movable] = previous or T == Real
comptime _HasReal[T: type_of(Tuple)] = T.Ts.reduce[False, _AnyReal]


def coefficients[P: Copyable & Deinitable](value: Specialization[P], mut budget: Budget, *, axis: Int = 0) raises PolyError -> Values[Real, 4]:
    """Extract ascending algebraic coefficients when only one unassigned axis remains.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        value: The value to convert or inspect.
        budget: The shared work allowance, including preparation, computation, and result storage.
        axis: The zero-based polynomial axis; defaults to zero.

    Returns:
        Ascending exact algebraic coefficients of the remaining univariate polynomial.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    from .sparse import _width
    from ._common import _sum
    if axis < 0 or axis >= _width[P]():
        raise PolyError(PolyErrorKind.DOMAIN, "coefficient axis is outside the power shape")
    var count = 0
    for monomial in value.evidence().monomials:
        for pair in monomial.pairs:
            budget.charge(1)
            if pair[0] != Int64(axis):
                raise PolyError(PolyErrorKind.DOMAIN, "specialized polynomial still has another free axis")
        count = max(count, _sum(Int(UInt64(monomial.exponent(Int64(axis), budget))), 1))
    _allocation[Real](count, budget)
    var zero = asreal(0, budget)
    var result = List[Real](length=count, fill=zero)
    for i in range(len(value)):
        budget.charge(1)
        var exponent = Int(UInt64(value.evidence().monomials[i].exponent(Int64(axis), budget)))
        result[exponent] = value[i][1].copy()
    return Values[Real, 4](result^, 0, budget)


def coefficients[P: Copyable & Deinitable](value: Specialization[P], *, axis: Int = 0, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Values[Real, 4]:
    """Materialize the remaining univariate algebraic coefficient sequence.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        value: The value to convert or inspect.
        axis: The zero-based polynomial axis; defaults to zero.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Ascending exact algebraic coefficients of the remaining univariate polynomial.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return coefficients(value, budget, axis=axis)


from polynomials import AlgebraicComparisonCertificate, check_algebraic_comparison


@fieldwise_init
struct _ComparisonData(Copyable, Movable):
    var left: RealAlgebraic
    var right: RealAlgebraic
    var proof: AlgebraicComparisonCertificate


@fieldwise_init
struct Comparison(Copyable, Movable):
    """An exact ordering relation retaining both values and independent evidence."""
    var _data: ArcPointer[_ComparisonData]

    def relation(self) -> Int:
        """Return the claimed exact ordering.

        Returns:
            Minus one, zero or one; use algebraic.check to verify the evidence.
        """
        return self._data[].proof.relation

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] AlgebraicComparisonCertificate:
        """Borrow the independent comparison witness.

        Returns:
            Immutable evidence bound to the original values.
        """
        return self._data[].proof


def compare_checked[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, mut budget: Budget) raises PolyError -> Comparison:
    """Compare exact real values and retain the ordering proof.

    Parameters:
        A: Inferred from the input types unless explicitly supplied.
        B: Inferred from the input types unless explicitly supplied.

    Args:
        a: The first polynomial.
        b: The second polynomial.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        A comparison retaining the original exact values and its claimed relation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var left = asreal(a, budget)
    var right = asreal(b, budget)
    var proof = algebraic_compare(left.witness(), right.witness(), budget)
    _allocation[_ComparisonData](1, budget)
    return Comparison(ArcPointer(_ComparisonData(_copy_witness(left.witness(), budget), _copy_witness(right.witness(), budget), proof^)))


def compare_checked[A: Copyable & Deinitable, B: Copyable & Deinitable](a: A, b: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Comparison:
    """Compare exact real values and retain the ordering proof.

    Parameters:
        A: Inferred from the input types unless explicitly supplied.
        B: Inferred from the input types unless explicitly supplied.

    Args:
        a: The first polynomial.
        b: The second polynomial.
        max_steps: The complete call’s work allowance.

    Returns:
        A comparison retaining the original exact values and its claimed relation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return compare_checked(a, b, budget)


def sign_checked[A: Copyable & Deinitable](value: A, mut budget: Budget) raises PolyError -> Comparison:
    """Determine the exact sign and retain its comparison with zero.

    Parameters:
        A: Inferred from the input types unless explicitly supplied.

    Args:
        value: The supplied value.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        A checked comparison against exact zero.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    return compare_checked(value, 0, budget)


def sign_checked[A: Copyable & Deinitable](value: A, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Comparison:
    """Determine the exact sign and retain its comparison with zero.

    Parameters:
        A: Inferred from the input types unless explicitly supplied.

    Args:
        value: The supplied value.
        max_steps: The complete call’s work allowance.

    Returns:
        A checked comparison against exact zero.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return sign_checked(value, budget)


def check(result: Comparison, mut budget: Budget) raises PolyError -> Bool:
    """Verify an exact algebraic comparison without trusting its relation flag.

    Args:
        result: The result retaining its original inputs and evidence.
        budget: The shared work allowance, including preparation and output storage.

    Returns:
        Whether the evidence establishes the retained ordering relation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    return check_algebraic_comparison(result._data[].left, result._data[].right, result._data[].proof, budget)


def check(result: Comparison, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify an exact algebraic comparison without trusting its relation flag.

    Args:
        result: The result retaining its original inputs and evidence.
        max_steps: The complete call’s work allowance.

    Returns:
        Whether the evidence establishes the retained ordering relation.

    Raises:
        PolyError: DOMAIN for incompatible input, UNSUPPORTED for an unavailable coefficient domain, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)


def ascomplex[A: Copyable & Deinitable, B: Copyable & Deinitable](real: A, imag: B, mut budget: Budget) raises PolyError -> Complex:
    """Build an exact complex algebraic value from its real and imaginary parts.

    Parameters:
        A: The real component's input type.
        B: The imaginary component's input type.

    Args:
        real: An integer, rational, or algebraic.Real.
        imag: An integer, rational, or algebraic.Real.
        budget: The work allowance shared by both conversions.

    Returns:
        The exact complex value, sharing existing real witnesses.

    Raises:
        PolyError: LIMIT or CAPACITY when conversion exceeds available resources.
    """
    return Complex(asreal(real, budget), asreal(imag, budget))


def ascomplex[A: Copyable & Deinitable, B: Copyable & Deinitable](real: A, imag: B, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Complex:
    """Build an exact complex value from two real components with one work limit.

    Parameters:
        A: The real component's input type.
        B: The imaginary component's input type.

    Args:
        real: An integer, rational, or algebraic.Real.
        imag: An integer, rational, or algebraic.Real.
        max_steps: The complete call's work allowance.

    Returns:
        The exact complex algebraic value.

    Raises:
        PolyError: LIMIT or CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return ascomplex(real, imag, budget)


def ascomplex[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> Complex:
    """Embed an exact scalar in the complex algebraic field.

    Parameters:
        T: The scalar's input type.

    Args:
        value: An integer, Rational, ExactComplex, Real, or Complex.
        budget: The shared work allowance.

    Returns:
        The promoted exact complex value.

    Raises:
        PolyError: LIMIT or CAPACITY on resource limits.
    """
    return _convert[5](value, NoContext(), budget)


def ascomplex[T: Copyable & Deinitable](value: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Complex:
    """Embed an exact scalar in the complex algebraic field with one work limit.

    Parameters:
        T: The scalar's input type.

    Args:
        value: An integer, Rational, ExactComplex, Real, or Complex.
        max_steps: The complete call's work allowance.

    Returns:
        The promoted exact complex value.

    Raises:
        PolyError: LIMIT or CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return ascomplex(value, budget)


def conjugate(value: Complex, mut budget: Budget) raises PolyError -> Complex:
    """Conjugate an exact complex value by negating its imaginary part.

    Args:
        value: The exact complex value.
        budget: The shared work allowance.

    Returns:
        The exact complex conjugate.

    Raises:
        PolyError: LIMIT or CAPACITY on resource limits.
    """
    budget.charge(1)
    return Complex(value.real(), value.imag().neg(budget))


def conjugate(value: Complex, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Complex:
    """Conjugate an exact complex value with one work limit.

    Args:
        value: The exact complex value.
        max_steps: The complete call's work allowance.

    Returns:
        The exact complex conjugate.

    Raises:
        PolyError: LIMIT or CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return conjugate(value, budget)


def check(value: Complex, mut budget: Budget) raises PolyError -> Bool:
    """Check both component witnesses and their retained arithmetic steps.

    Args:
        value: The exact complex value to verify.
        budget: The shared checking allowance.

    Returns:
        Whether both real component checks pass.

    Raises:
        PolyError: DOMAIN for malformed evidence, or LIMIT/CAPACITY on resource limits.
    """
    return check(value.real(), budget) and check(value.imag(), budget)


def check(value: Complex, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Check both exact complex components with one work limit.

    Args:
        value: The exact complex value to verify.
        max_steps: The complete check's work allowance.

    Returns:
        Whether both real component checks pass.

    Raises:
        PolyError: DOMAIN for malformed evidence, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(value, budget)


comptime _ComplexInput[T: AnyType] = T == Complex or T == ExactComplex
comptime _Arithmetic[A: AnyType, B: AnyType] = Complex if _ComplexInput[A] or _ComplexInput[B] else Real


from ._common import _Scalar, _Ring, _values_from_core

def coefficients[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](value: I, mut budget: Budget, *, axis: Int = 0) raises PolyError -> Values[_Scalar[_Kind[C]], _Kind[C]] where _Element[I] == Tuple[P, C]:
    """Extract one remaining coefficient axis from algebraic sparse terms.

    Parameters:
        I: The sparse input iterable type.
        P: The integer or tuple power shape.
        C: The real or complex algebraic coefficient type.

    Args:
        value: Algebraic sparse terms or a prepared sparse array.
        budget: The shared work allowance, including dense output allocation.
        axis: The only axis permitted to have nonzero powers.

    Returns:
        Ascending coefficients in AA or QQbar, including any missing zeros.

    Raises:
        PolyError: DOMAIN for an invalid axis or another free variable, or LIMIT/CAPACITY on resource limits.
    """
    from .sparse import _prepare, _axis
    comptime k = _Kind[C]
    comptime assert k == 4 or k == 5, "algebraic coefficient extraction requires Real or algebraic.Complex values"
    _axis[P](axis)
    var p = _prepare[I, P, C](value, budget)
    var core = DensePoly[_Ring[k]].from_sparse(p._core[], budget, variable=Optional(Int64(axis)))
    var values = List[_Ring[k]]()
    # from_sparse allocates one contiguous list at its known output length.
    debug_assert(core.coefficients.is_contiguous(), "dense conversion must be contiguous")
    swap(values, core.coefficients._flat)
    return _values_from_core[k](values^, core.ring, budget)


def coefficients[I: Iterable, P: Copyable & Deinitable, C: Copyable & Deinitable](value: I, *, axis: Int = 0, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Values[_Scalar[_Kind[C]], _Kind[C]] where _Element[I] == Tuple[P, C]:
    """Extract one algebraic coefficient axis with one work allowance.

    Parameters:
        I: The sparse input iterable type.
        P: The integer or tuple power shape.
        C: The real or complex algebraic coefficient type.

    Args:
        value: Algebraic sparse terms or a prepared sparse array.
        axis: The only axis permitted to have nonzero powers.
        max_steps: The complete call's work allowance.

    Returns:
        Ascending coefficients in AA or QQbar, retaining missing zeros.

    Raises:
        PolyError: DOMAIN for an invalid axis or another free variable, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return coefficients[I, P, C](value, budget, axis=axis)


from . import _common

from . import _specialization as _specialized
from ._specialization import AlgebraicSpecialization


def check[P: Copyable & Deinitable](result: AlgebraicSpecialization[P], mut budget: Budget) raises PolyError -> Bool:
    """Check every specialized coefficient and bind the proof to the returned terms.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        budget: The shared work allowance, including preparation, computation, and result storage.

    Returns:
        Whether the evidence binds the specialized coefficients to the original
        polynomial and assignments.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    return _specialized._check(result, budget)


def check[P: Copyable & Deinitable](result: AlgebraicSpecialization[P], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Bool:
    """Verify a specialization against its retained input and assignments.

    Parameters:
        P: The sparse power shape or retained source array type, inferred from the input.

    Args:
        result: The result whose retained evidence or values are being inspected.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Whether the evidence binds the specialized coefficients to the original
        polynomial and assignments.

    Raises:
        PolyError: DOMAIN for invalid or incompatible input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return check(result, budget)
