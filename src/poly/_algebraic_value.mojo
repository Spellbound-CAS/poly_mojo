"""Immutable real and complex algebraic coefficient fields behind the array API."""
from std.memory import ArcPointer
from std.format import Writer
from apn_mojo import Rational
from polynomials import Budget, PolyError, PolyErrorKind, QQ, NoContext, DensePoly, SparsePoly, CoefficientField, isolate_real_roots
from polynomials import RealAlgebraic, AlgebraicArithmeticCertificate, AlgebraicEvaluationCertificate, AlgebraicPolynomialEvaluation
from polynomials import real_root, evaluate_dense, algebraic_add, algebraic_subtract, algebraic_multiply, algebraic_negate, algebraic_inverse
from std.sys import size_of
from polynomials import check_real_algebraic, check_algebraic_arithmetic, algebraic_compare, check_algebraic_comparison

def _charge_size(n: Int, mut budget: Budget) raises PolyError:
    if n < 0:
        raise PolyError(PolyErrorKind.CAPACITY, "algebraic witness size overflow")
    budget.charge(UInt64(n))

def _sum(a: Int, b: Int) raises PolyError -> Int:
    if a < 0 or b < 0 or a > Int.MAX - b:
        raise PolyError(PolyErrorKind.CAPACITY, "algebraic witness size overflow")
    return a + b

def _allocation[T: AnyType](count: Int, mut budget: Budget) raises PolyError:
    if count < 0 or count > Int.MAX // size_of[T]():
        raise PolyError(PolyErrorKind.CAPACITY, "algebraic allocation size overflow")
    _charge_size(count, budget)

@fieldwise_init
struct _Canonical(Copyable, Movable):
    var polynomial: DensePoly[QQ]
    var index: Int
    var zero: Bool
    var bytes: Int

def _root_index(value: RealAlgebraic, mut budget: Budget) raises PolyError -> Int:
    # The number of real roots of the minimal polynomial below the value: its
    # isolating intervals are disjoint, and narrowing the value by sign until
    # each one holds it or lies clear of it leaves exactly one holding it.
    var isolation = isolate_real_roots(value.polynomial, budget)
    var low = value.low.copy()
    var high = value.high.copy()
    var low_sign = evaluate_dense(value.polynomial, low, budget).value.sign()
    var half = QQ(Rational(2)).inverse(budget)
    while True:
        budget.charge(1)
        var below = 0
        var settled = True
        for q in isolation.exact_roots:
            if q.sub(low, budget).value.sign() < 0:
                below += 1
        for interval in isolation.intervals:
            if interval.high.sub(low, budget).value.sign() <= 0:
                below += 1
            elif interval.low.sub(high, budget).value.sign() < 0 and not (interval.low.sub(low, budget).value.sign() <= 0 and high.sub(interval.high, budget).value.sign() <= 0):
                settled = False
        if settled:
            return below
        var middle = low.add(high, budget).mul(half, budget)
        if evaluate_dense(value.polynomial, middle, budget).value.sign() == low_sign:
            low = middle^
        else:
            high = middle^


def _make(var data: _ValueData, mut budget: Budget) raises PolyError -> Real:
    # Normalize at construction so fixed-signature equality never runs a solver.
    if not check_real_algebraic(data.plain, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic coefficient witness")
    var value = real_root(data.plain.polynomial, data.plain.low, data.plain.high, budget)
    if not check_real_algebraic(value, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic normalization produced an invalid witness")
    if not _same_witness(value, data.plain, budget):
        var comparison = algebraic_compare(value, data.plain, budget)
        if comparison.relation != 0 or not check_algebraic_comparison(value, data.plain, comparison, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "algebraic normalization changed the coefficient")
    var index = 0 if value.rational else _root_index(value, budget)
    return _finish(value^, index, data^, budget)


def _normal(var value: RealAlgebraic, index: Int, mut budget: Budget) raises PolyError -> Real:
    """A Real from a witness already on its monic irreducible polynomial, the index-th real root of it.

    Library producers that isolated every real root of that polynomial in
    order know the index and skip renormalization.
    """
    _allocation[_ValueData](1, budget)
    var data = _ValueData(value.copy(), None, None, None, None, None)
    return _finish(value^, 0 if value.rational else index, data^, budget)


def _finish(var value: RealAlgebraic, index: Int, var data: _ValueData, mut budget: Budget) raises PolyError -> Real:
    var size = _sum(value.low.size(), value.high.size())
    for coefficient_index in range(len(value.polynomial.coefficients)):
        ref c = value.polynomial.coefficients[coefficient_index]
        _charge_size(c.size(), budget)
        size = _sum(size, c.size())
    _allocation[_Canonical](1, budget)
    _allocation[_ValueData](1, budget)
    var zero = value.rational and value.low.is_zero()
    var polynomial = DensePoly[QQ](0, NoContext())
    swap(polynomial, value.polynomial)
    return Real(ArcPointer(_Canonical(polynomial^, index, zero, size)), ArcPointer(data^), None, None)

def _constant(value: Int) -> Real:
    # Only zero and one call this fixed-size protocol constructor.
    var p = DensePoly[QQ](0, NoContext())
    var values: List[QQ] = [QQ(Rational(-value)), QQ(Rational(1))]
    p.coefficients = values^
    var witness = RealAlgebraic(p^, QQ(Rational(value)), QQ(Rational(value)), True)
    return Real(ArcPointer(_Canonical(witness.polynomial.copy(), 0, value == 0, 8)), ArcPointer(_ValueData(witness^, None, None, None, None, None)), None, None)

def _copy_witness(value: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    _allocation[QQ](len(value.polynomial.coefficients), budget)
    for coefficient_index in range(len(value.polynomial.coefficients)):
        ref c = value.polynomial.coefficients[coefficient_index]
        _charge_size(c.size(), budget)
    _charge_size(value.low.size(), budget)
    _charge_size(value.high.size(), budget)
    return value.copy()


def _same_witness(a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    return a.rational == b.rational and a.low.sub(b.low, budget).is_zero() and a.high.sub(b.high, budget).is_zero() and a.polynomial.sub(b.polynomial, budget).is_zero()


@fieldwise_init
struct _HornerData(Copyable, Movable):
    var coefficients: List[RealAlgebraic]
    var point: RealAlgebraic
    var proof: AlgebraicEvaluationCertificate


@fieldwise_init
struct _SparseEvaluationData(Copyable, Movable):
    var polynomial: SparsePoly[QQ]
    var point: List[Tuple[Int64, RealAlgebraic]]
    var proof: AlgebraicPolynomialEvaluation


@fieldwise_init
struct _ValueData(Copyable, Movable):
    var plain: RealAlgebraic
    var arithmetic: Optional[AlgebraicArithmeticCertificate]
    var left: Optional[RealAlgebraic]
    var right: Optional[RealAlgebraic]
    var divisor: Optional[RealAlgebraic]
    var division: Optional[AlgebraicArithmeticCertificate]


@fieldwise_init
struct Real(CoefficientField):
    """An exact real number sharing its root witness and any arithmetic proof.

    Use asreal or a polynomial's real_roots function to construct a value.
    Arithmetic stays exact; bounds and approximate make approximation explicit.
    """
    var _canonical: ArcPointer[_Canonical]
    var _data: ArcPointer[_ValueData]
    var _horner: Optional[ArcPointer[_HornerData]]
    var _sparse_evaluation: Optional[ArcPointer[_SparseEvaluationData]]

    def witness(self) -> ref[ImmOrigin(origin_of(self._data[].plain))] RealAlgebraic:
        """Borrow the exact value's defining polynomial and root witness.

        Returns:
            Immutable core evidence tied to this value's shared storage.
        """
        return self._data[].plain


    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].arithmetic))] Optional[AlgebraicArithmeticCertificate]:
        """Borrow the arithmetic proof, when this value came from arithmetic.

        Returns:
            The optional immutable arithmetic certificate.
        """
        return self._data[].arithmetic

    def division_evidence(self) -> ref[ImmOrigin(origin_of(self._data[].division))] Optional[AlgebraicArithmeticCertificate]:
        """Borrow the inverse proof retained by divide.

        Returns:
            An inverse certificate when this result came from divide.
        """
        return self._data[].division



    comptime Context = NoContext
    """The context shared by all real algebraic values."""

    @staticmethod
    def zero(ctx: NoContext) -> Self:
        """Construct exact zero in the real algebraic field.

        Args:
            ctx: The field context.

        Returns:
            The exact immutable value in this coefficient field.
        """
        return _constant(0)

    @staticmethod
    def one(ctx: NoContext) -> Self:
        """Construct exact one in the real algebraic field.

        Args:
            ctx: The field context.

        Returns:
            The exact immutable value in this coefficient field.
        """
        return _constant(1)

    def context(self) -> NoContext:
        """Return the real algebraic field's unique context.

        Returns:
            The unique field context.
        """
        return NoContext()

    def is_zero(self) -> Bool:
        """Test zero using the normalized rational representation.

        Returns:
            Whether the exact condition holds.
        """
        return self._canonical[].zero

    def size(self) -> Int:
        """Return the cached witness size used for work accounting.

        Returns:
            The stored size used for precharging work.
        """
        return self._canonical[].bytes

    def __eq__(self, other: Self) -> Bool:
        """Compare minimal polynomials and real-root indices exactly.

        Args:
            other: The other coefficient in this field.

        Returns:
            Whether the exact condition holds.
        """
        return self._canonical[].index == other._canonical[].index and self._canonical[].polynomial == other._canonical[].polynomial

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Add real algebraic coefficients with retained arithmetic evidence.

        Args:
            other: The other coefficient in this field.
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return _field_arithmetic(self, Optional(other.copy()), 0, budget)

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Subtract real algebraic coefficients exactly.

        Args:
            other: The other coefficient in this field.
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return _field_arithmetic(self, Optional(other.copy()), 1, budget)

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Multiply real algebraic coefficients exactly.

        Args:
            other: The other coefficient in this field.
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return _field_arithmetic(self, Optional(other.copy()), 2, budget)

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Negate a real algebraic coefficient exactly.

        Args:
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return _field_arithmetic(self, Optional[Self](), 3, budget)

    def inverse(self, mut budget: Budget) raises PolyError -> Self:
        """Invert a nonzero coefficient; zero raises DOMAIN.

        Args:
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return _field_arithmetic(self, Optional[Self](), 4, budget)

    def write_to[W: Writer](self, mut writer: W):
        """Write an exact rational or a defining polynomial and isolating interval.

        Parameters:
            W: The output writer type.

        Args:
            writer: The destination for the exact representation.
        """
        if self.witness().rational:
            writer.write(self.witness().low)
        else:
            writer.write("Real(", self.witness().polynomial, ", ", self.witness().low, ", ", self.witness().high, ")")


def _plain(var value: RealAlgebraic, mut budget: Budget) raises PolyError -> Real:
    _allocation[_ValueData](1, budget)
    return _make(_ValueData(value^, None, None, None, None, None), budget)


def _arithmetic(a: Real, b: Optional[Real], operation: Int, mut budget: Budget) raises PolyError -> Real:
    var left = _copy_witness(a.witness(), budget)
    var right = Optional[RealAlgebraic]()
    if b:
        right = _copy_witness(b.value().witness(), budget)
    var proof: AlgebraicArithmeticCertificate
    if operation == 0:
        proof = algebraic_add(left, right.value(), budget)
    elif operation == 1:
        proof = algebraic_subtract(left, right.value(), budget)
    elif operation == 2:
        proof = algebraic_multiply(left, right.value(), budget)
    elif operation == 3:
        proof = algebraic_negate(left, budget)
    else:
        proof = algebraic_inverse(left, budget)
    _allocation[_ValueData](1, budget)
    return _make(_ValueData(_copy_witness(proof.value, budget), Optional(proof^), Optional(left^), right^, None, None), budget)


@fieldwise_init
struct Complex(CoefficientField):
    """An exact complex algebraic number represented by two real algebraic parts.

    Construct values with algebraic.ascomplex(real, imag). No precision context
    or chosen primitive element is required.
    """
    var _real: Real
    var _imag: Real
    comptime Context = NoContext
    """The context shared by all complex algebraic values."""

    @staticmethod
    def zero(ctx: NoContext) -> Self:
        """Construct exact complex zero.

        Args:
            ctx: The field context.

        Returns:
            The exact immutable value in this coefficient field.
        """
        return Self(Real.zero(ctx), Real.zero(ctx))

    @staticmethod
    def one(ctx: NoContext) -> Self:
        """Construct exact complex one.

        Args:
            ctx: The field context.

        Returns:
            The exact immutable value in this coefficient field.
        """
        return Self(Real.one(ctx), Real.zero(ctx))

    def context(self) -> NoContext:
        """Return the complex algebraic field's unique context.

        Returns:
            The unique field context.
        """
        return NoContext()

    def real(self) -> Real:
        """Return the immutable exact real part.

        Returns:
            The exact immutable value in this coefficient field.
        """
        return self._real.copy()

    def imag(self) -> Real:
        """Return the immutable exact imaginary part.

        Returns:
            The exact immutable value in this coefficient field.
        """
        return self._imag.copy()

    def is_zero(self) -> Bool:
        """Test whether both real parts are exactly zero.

        Returns:
            Whether the exact condition holds.
        """
        return self._real.is_zero() and self._imag.is_zero()

    def size(self) -> Int:
        """Return the larger component witness size; copies share both components.

        Returns:
            The stored size used for precharging work.
        """
        return max(self._real.size(), self._imag.size())

    def __eq__(self, other: Self) -> Bool:
        """Compare both algebraic components exactly.

        Args:
            other: The other coefficient in this field.

        Returns:
            Whether the exact condition holds.
        """
        return self._real == other._real and self._imag == other._imag

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Add exact complex algebraic coefficients.

        Args:
            other: The other coefficient in this field.
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return Self(self._real.add(other._real, budget), self._imag.add(other._imag, budget))

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Subtract exact complex algebraic coefficients.

        Args:
            other: The other coefficient in this field.
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return Self(self._real.sub(other._real, budget), self._imag.sub(other._imag, budget))

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Multiply exact complex algebraic coefficients.

        Args:
            other: The other coefficient in this field.
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        var re = self._real.mul(other._real, budget).sub(self._imag.mul(other._imag, budget), budget)
        var im = self._real.mul(other._imag, budget).add(self._imag.mul(other._real, budget), budget)
        return Self(re^, im^)

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Negate both exact algebraic components.

        Args:
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        return Self(self._real.neg(budget), self._imag.neg(budget))

    def inverse(self, mut budget: Budget) raises PolyError -> Self:
        """Invert using the squared modulus; zero raises DOMAIN.

        Args:
            budget: The caller-owned work allowance.

        Returns:
            The exact immutable value in this coefficient field.

        Raises:
            PolyError: DOMAIN for division by zero, or LIMIT/CAPACITY on resource limits.
        """
        var norm = self._real.mul(self._real, budget).add(self._imag.mul(self._imag, budget), budget)
        var inverse = norm.inverse(budget)
        return Self(self._real.mul(inverse, budget), self._imag.neg(budget).mul(inverse, budget))

    def write_to[W: Writer](self, mut writer: W):
        """Write both exact algebraic components.

        Parameters:
            W: The output writer type.

        Args:
            writer: The destination for the exact representation.
        """
        writer.write("Complex(", self._real, ", ", self._imag, ")")


def _field_arithmetic(a: Real, b: Optional[Real], operation: Int, mut budget: Budget) raises PolyError -> Real:
    # Generic polynomial checkers call coefficient arithmetic too. Verify the
    # scalar certificate here so their identities cannot accept an unchecked
    # elimination result as a primitive field operation.
    var result = _arithmetic(a, b, operation, budget)
    ref data = result._data[]
    if data.arithmetic.value().operation != operation or not check_algebraic_arithmetic(data.left.value(), data.right, data.arithmetic.value(), budget):
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic coefficient arithmetic failed verification")
    return result^
