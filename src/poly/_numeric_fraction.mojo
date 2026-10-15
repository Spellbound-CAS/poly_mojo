"""Numerical fractions and independently checked quotient enclosures."""
from std.memory import ArcPointer
from apn_mojo import Rational, BallContext, Float, Complex, Ball, ComplexBall, divide
from polynomials import Budget, PolyError, PolyErrorKind, NumericDomain, NumericBox, QQ, numeric_bounds, numeric_value, FloatDomain, ComplexDomain, BallDomain
from . import sparse, evaluation
from ._io_meta import _Record, _Fraction
from ._archive_meta import _Archive, _QuotientInput
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink
from . import _numeric as _num
from ._common import _allocation, _charge_size, _sum


@fieldwise_init
struct NumericFunction[P: Copyable & Deinitable, D: NumericDomain](Copyable, Movable, _Record):
    """A numerical numerator and denominator with explicit arithmetic contexts.

    No approximate GCD cancellation is performed. Ball coefficients describe
    independent families; evaluation requires a denominator enclosure excluding zero.

    Parameters:
        P: The retained polynomial power shape.
        D: The inferred numerical coefficient domain.
    """
    comptime _Format = _Fraction
    comptime _Item = Self.D.Scalar
    comptime _Power = Self.P
    var _numerator: sparse.NumericArray[Self.P, Self.D]
    var _denominator: sparse.NumericArray[Self.P, Self.D]

    def ndim(self) -> Int:
        """Return the number of polynomial axes.

        Returns:
            The retained axis count, including unused axes.
        """
        return sparse._width[Self.P]()


def _make[P: Copyable & Deinitable, D: NumericDomain](var p: sparse.NumericArray[P, D], var q: sparse.NumericArray[P, D], mut budget: Budget) raises PolyError -> NumericFunction[P, D]:
    # A coefficient excluding zero ensures every represented denominator is nonzero.
    var nonzero = False
    for term in q._core[].terms:
        if not _zero(numeric_bounds[D](term.coefficient, budget)):
            nonzero = True
    if not nonzero:
        raise PolyError(PolyErrorKind.DOMAIN, "denominator may be the zero polynomial")
    return NumericFunction[P, D](p^, q^)


def _zero(box: NumericBox) -> Bool:
    return box.real_low.value.sign() <= 0 and box.real_high.value.sign() >= 0 and box.imag_low.value.sign() <= 0 and box.imag_high.value.sign() >= 0


def _divide[D: NumericDomain](a: D.Scalar, b: D.Scalar, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    _num._context[D](context, budget)
    var left = numeric_bounds[D](a, budget)
    var right = numeric_bounds[D](b, budget)
    if _zero(right):
        raise PolyError(PolyErrorKind.DOMAIN, "denominator enclosure contains zero")
    _charge_size(_sum(_sum(left.real_low.size(), left.imag_low.size()), _sum(right.real_low.size(), right.imag_low.size())), budget)
    var result: D.Scalar
    try:
        comptime if D == FloatDomain:
            result = rebind_var[D.Scalar](divide(rebind[Float](a), rebind[Float](b), context=rebind[_num.ArithmeticContext](context)))
        elif D == ComplexDomain:
            result = rebind_var[D.Scalar](divide(rebind[Complex](a), rebind[Complex](b), context=rebind[_num.ComplexContext](context)))
        elif D == BallDomain:
            result = rebind_var[D.Scalar](divide(rebind[Ball](a), rebind[Ball](b), context=rebind[BallContext](context)))
        else:
            result = rebind_var[D.Scalar](divide(rebind[ComplexBall](a), rebind[ComplexBall](b), context=rebind[BallContext](context)))
    except e:
        raise PolyError(PolyErrorKind.DOMAIN, String(e))

    _ = numeric_bounds[D](result, budget)
    return result^


def _product(a: QQ, b: QQ, c: QQ, d: QQ, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    var low = a.mul(c, budget)
    var high = low.copy()
    for v in [a.mul(d, budget), b.mul(c, budget), b.mul(d, budget)]:
        _charge_size(v.size(), budget)
        if v.sub(low, budget).value.sign() < 0: low = v.copy()
        if v.sub(high, budget).value.sign() > 0: high = v.copy()
    return low^, high^


def _square(a: QQ, b: QQ, mut budget: Budget) raises PolyError -> Tuple[QQ, QQ]:
    var aa = a.mul(a, budget)
    var bb = b.mul(b, budget)
    var low = aa.copy()
    var high = bb.copy()
    if aa.sub(bb, budget).value.sign() > 0:
        swap(low, high)
    if a.value.sign() <= 0 and b.value.sign() >= 0:
        low = QQ(Rational(0))
    return low^, high^


def _quotient(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    if _zero(b):
        raise PolyError(PolyErrorKind.DOMAIN, "denominator enclosure contains zero")
    var re = _square(b.real_low, b.real_high, budget)
    var im = _square(b.imag_low, b.imag_high, budget)
    var low = re[0].add(im[0], budget)
    var high = re[1].add(im[1], budget)
    var invlow = high.inverse(budget)
    var invhigh = low.inverse(budget)
    var ac = _product(a.real_low, a.real_high, b.real_low, b.real_high, budget)
    var bd = _product(a.imag_low, a.imag_high, b.imag_low, b.imag_high, budget)
    var bc = _product(a.imag_low, a.imag_high, b.real_low, b.real_high, budget)
    var ad = _product(a.real_low, a.real_high, b.imag_low, b.imag_high, budget)
    var real = _product(ac[0].add(bd[0], budget), ac[1].add(bd[1], budget), invlow, invhigh, budget)
    var imag = _product(bc[0].sub(ad[1], budget), bc[1].sub(ad[0], budget), invlow, invhigh, budget)
    return NumericBox(real[0].copy(), real[1].copy(), imag[0].copy(), imag[1].copy())


@fieldwise_init
struct Enclosure[E: _QuotientInput, D: NumericDomain](_Archive, Copyable, Movable):
    """A quotient enclosure retaining checked numerator and denominator evaluations.

    Parameters:
        E: The retained polynomial evaluation type.
        D: The real or complex ball domain.
    """
    var _numerator: Self.E
    var _denominator: Self.E
    var _value: ArcPointer[Self.D.Scalar]

    def value(self, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Copy the certified enclosure using the shared work allowance.

        Args:
            budget: The allowance for copying the enclosing scalar.

        Returns:
            The real or complex ball quotient.

        Raises:
            PolyError: LIMIT or CAPACITY on resource limits.
        """
        _ = numeric_bounds[Self.D](self._value[], budget)
        return self._value[].copy()

    def value(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Self.D.Scalar:
        """Copy the enclosing value with a bounded work allowance.

        Args:
            max_steps: The copying allowance.

        Returns:
            The real or complex ball quotient.

        Raises:
            PolyError: LIMIT or CAPACITY on resource limits.
        """
        var budget = Budget(max_steps)
        return self.value(budget)

    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return _check(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io as codec
        codec._depth(depth, budget)
        codec._write_node[Self.E, QQ, Self.D, Int](self._numerator, writer, sink, budget, depth + 1)
        codec._write_node[Self.E, QQ, Self.D, Int](self._denominator, writer, sink, budget, depth + 1)
        codec._write_node[Self.D.Scalar, QQ, Self.D, Int](self._value[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io as codec
        codec._depth(depth, budget)
        var a = codec._read_node[Self.E, QQ, Self.D, Int](reader, budget, depth + 1)
        var b = codec._read_node[Self.E, QQ, Self.D, Int](reader, budget, depth + 1)
        var value = codec._read_node[Self.D.Scalar, QQ, Self.D, Int](reader, budget, depth + 1)
        _allocation[Self.D.Scalar](1, budget)
        return Self(a^, b^, ArcPointer(value^))


def _enclose[C: AnyType, D: NumericDomain, numerical: Bool](var a: evaluation.Enclosure[C, D, numerical, True], var b: evaluation.Enclosure[C, D, numerical, True], context: BallContext, mut budget: Budget) raises PolyError -> Enclosure[evaluation.Enclosure[C, D, numerical, True], D]:
    var bounds = _quotient(numeric_bounds[D](a.value(budget), budget), numeric_bounds[D](b.value(budget), budget), budget)
    var value = numeric_value[D](bounds, rebind[D.Context](context), budget)
    _allocation[D.Scalar](1, budget)
    return Enclosure[evaluation.Enclosure[C, D, numerical, True], D](a^, b^, ArcPointer(value^))


def _check[E: _QuotientInput, D: NumericDomain](value: Enclosure[E, D], mut budget: Budget) raises PolyError -> Bool:
    if not value._numerator._archive_check(budget) or not value._denominator._archive_check(budget):
        return False
    if not value._numerator._same_point(value._denominator, budget): return False
    var a = value._numerator._quotient_bounds(budget)
    var b = value._denominator._quotient_bounds(budget)
    if _zero(b): return False
    return numeric_bounds[D](value._value[], budget).contains(_quotient(a, b, budget), budget)


from polynomials import NumericOptions, NumericOutcome, NoContext
from apn_mojo import Integer


@fieldwise_init
struct AdaptiveEnclosure[E: _QuotientInput, D: NumericDomain](_Archive, Copyable, Movable):
    """A checked quotient with an accuracy outcome and optional finite enclosure.

    Parameters:
        E: The retained polynomial evaluation type.
        D: The real or complex ball domain.
    """
    var _numerator: Self.E
    var _denominator: Self.E
    var _value: Optional[ArcPointer[Self.D.Scalar]]
    var _options: NumericOptions
    var _precision: Int
    var _outcome: NumericOutcome

    def has_value(self) -> Bool:
        """Report whether the denominator was certified away from zero.

        Returns:
            Whether value() can return a finite enclosure.
        """
        return Bool(self._value)

    def accurate(self) -> Bool:
        """Report whether the requested quotient width was established.

        Returns:
            True for SUCCESS; verify the evidence with rational.check.
        """
        return self._outcome == NumericOutcome.SUCCESS

    def precision(self) -> Int:
        """Return the last working precision in bits.

        Returns:
            The precision used for both polynomial evaluations and the quotient.
        """
        return self._precision

    def options(self) -> NumericOptions:
        """Return the original accuracy request.

        Returns:
            The retained precision limits and width targets.
        """
        return self._options

    def outcome(self) -> NumericOutcome:
        """Report SUCCESS or PRECISION_CEILING.

        Returns:
            PRECISION_CEILING when accuracy or denominator separation remains unresolved.
        """
        return self._outcome

    def value(self, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Copy the finite quotient enclosure when one was established.

        Args:
            budget: The copying allowance.

        Returns:
            An enclosing ball, which may be wider than requested at the precision ceiling.

        Raises:
            PolyError: DOMAIN if the denominator still contains zero, or LIMIT/CAPACITY.
        """
        budget.charge(1)
        if not self._value:
            raise PolyError(PolyErrorKind.DOMAIN, "denominator separation remains unresolved")
        _ = numeric_bounds[Self.D](self._value.value()[], budget)
        return self._value.value()[].copy()

    def value(self, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Self.D.Scalar:
        """Copy the finite quotient enclosure under a work limit.

        Args:
            max_steps: The copying allowance.

        Returns:
            The enclosing real or complex ball.

        Raises:
            PolyError: DOMAIN without a finite enclosure, or LIMIT/CAPACITY.
        """
        var budget = Budget(max_steps)
        return self.value(budget)

    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return _check_adaptive(self, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io as codec
        codec._depth(depth, budget)
        codec._write_node[Self.E, QQ, Self.D, Int](self._numerator, writer, sink, budget, depth + 1)
        codec._write_node[Self.E, QQ, Self.D, Int](self._denominator, writer, sink, budget, depth + 1)
        codec._write_node[Bool, QQ, Self.D, Int](Bool(self._value), writer, sink, budget, depth + 1)
        if self._value:
            codec._write_node[Self.D.Scalar, QQ, Self.D, Int](self._value.value()[], writer, sink, budget, depth + 1)
        codec._write_node[NumericOptions, QQ, Self.D, Int](self._options, writer, sink, budget, depth + 1)
        codec._write_node[Int, QQ, Self.D, Int](self._precision, writer, sink, budget, depth + 1)
        codec._write_node[NumericOutcome, QQ, Self.D, Int](self._outcome, writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io as codec
        codec._depth(depth, budget)
        var a = codec._read_node[Self.E, QQ, Self.D, Int](reader, budget, depth + 1)
        var b = codec._read_node[Self.E, QQ, Self.D, Int](reader, budget, depth + 1)
        var value = Optional[ArcPointer[Self.D.Scalar]]()
        if codec._read_node[Bool, QQ, Self.D, Int](reader, budget, depth + 1):
            var scalar = codec._read_node[Self.D.Scalar, QQ, Self.D, Int](reader, budget, depth + 1)
            _allocation[Self.D.Scalar](1, budget)
            value = ArcPointer(scalar^)
        var options = codec._read_node[NumericOptions, QQ, Self.D, Int](reader, budget, depth + 1)
        var precision = codec._read_node[Int, QQ, Self.D, Int](reader, budget, depth + 1)
        var outcome = codec._read_node[NumericOutcome, QQ, Self.D, Int](reader, budget, depth + 1)
        return Self(a^, b^, value^, options, precision, outcome)


def _adaptive_step[C: AnyType, D: NumericDomain, numerical: Bool](var a: evaluation.Enclosure[C, D, numerical, True], var b: evaluation.Enclosure[C, D, numerical, True], options: NumericOptions, context: BallContext, mut budget: Budget) raises PolyError -> AdaptiveEnclosure[evaluation.Enclosure[C, D, numerical, True], D]:
    var value = Optional[ArcPointer[D.Scalar]]()
    var outcome = NumericOutcome.PRECISION_CEILING
    var denominator = b._quotient_bounds(budget)
    if not _zero(denominator):
        var bounds = _quotient(a._quotient_bounds(budget), denominator, budget)
        var scalar = numeric_value[D](bounds, rebind[D.Context](context), budget)
        if _accurate(numeric_bounds[D](scalar, budget), options, budget):
            outcome = NumericOutcome.SUCCESS
        _allocation[D.Scalar](1, budget)
        value = ArcPointer(scalar^)
    return AdaptiveEnclosure[evaluation.Enclosure[C, D, numerical, True], D](a^, b^, value^, options, context.precision().value(), outcome)


def _check_adaptive[E: _QuotientInput, D: NumericDomain](value: AdaptiveEnclosure[E, D], mut budget: Budget) raises PolyError -> Bool:
    value._options._validate()
    if not value._numerator._archive_check(budget) or not value._denominator._archive_check(budget):
        return False
    if not value._numerator._same_point(value._denominator, budget):
        return False
    var bits = value._options.initial_precision
    while bits < value._precision and bits < value._options.max_precision:
        budget.charge(1)
        bits = value._options.max_precision if bits > value._options.max_precision // 2 else bits * 2
    if bits != value._precision or value._numerator._quotient_precision() != bits or value._denominator._quotient_precision() != bits:
        return False
    var denominator = value._denominator._quotient_bounds(budget)
    var accurate = False
    if value._value:
        if _zero(denominator):
            return False
        var box = numeric_bounds[D](value._value.value()[], budget)
        if not box.contains(_quotient(value._numerator._quotient_bounds(budget), denominator, budget), budget):
            return False
        accurate = _accurate(box, value._options, budget)
    elif not _zero(denominator):
        return False
    if value._outcome == NumericOutcome.SUCCESS:
        return accurate
    return value._outcome == NumericOutcome.PRECISION_CEILING and not accurate and bits == value._options.max_precision


def _dyadic(bits: Int, mut budget: Budget) raises PolyError -> QQ:
    _charge_size(bits, budget)
    try:
        return QQ(Rational(Integer(1), Integer(1) << bits))
    except error:
        raise PolyError(PolyErrorKind.CAPACITY, String(error))


def _accurate(box: NumericBox, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    # The relative scale is a certified lower bound on the infinity norm.
    # A rectangle crossing both axes has scale zero, not its midpoint's norm.
    var scale = QQ.zero(NoContext())
    for bounds in [(box.real_low.copy(), box.real_high.copy()), (box.imag_low.copy(), box.imag_high.copy())]:
        var lower = QQ.zero(NoContext())
        if bounds[0].value.sign() > 0:
            lower = bounds[0].copy()
        elif bounds[1].value.sign() < 0:
            lower = bounds[1].neg(budget)
        if lower.sub(scale, budget).value.sign() > 0:
            scale = lower^
    var target = _dyadic(options.absolute_bits, budget) if options.absolute_bits else QQ.zero(NoContext())
    if options.relative_bits:
        target = target.add(_dyadic(options.relative_bits, budget).mul(scale, budget), budget)
    return box.real_high.sub(box.real_low, budget).sub(target, budget).value.sign() <= 0 and box.imag_high.sub(box.imag_low, budget).sub(target, budget).value.sign() <= 0
