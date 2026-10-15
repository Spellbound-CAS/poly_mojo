"""Explicit apn_mojo numerical domains, separate from exact coefficient rings."""
from apn_mojo import Integer, Float, Complex, Ball, ComplexBall, ArithmeticContext, ComplexContext, BallContext, add, subtract, multiply, divide, contains_zero
from .errors import Budget, PolyError, PolyErrorKind, _sum, _product, _charge_size


trait NumericDomain:
    """Associate a numerical scalar with its explicit arithmetic context."""
    comptime Scalar: Copyable & Writable & Deinitable
    """The underlying apn_mojo scalar."""
    comptime Context: Copyable & Writable & Deinitable
    """The explicit output context."""


struct FloatDomain(NumericDomain):
    """Rounded arbitrary-precision real arithmetic."""
    comptime Scalar = Float
    """The apn_mojo point value."""
    comptime Context = ArithmeticContext
    """Precision, rounding mode and traps."""


struct ComplexDomain(NumericDomain):
    """Rounded arbitrary-precision complex arithmetic."""
    comptime Scalar = Complex
    """The apn_mojo complex point value."""
    comptime Context = ComplexContext
    """Independent real and imaginary output settings."""


struct BallDomain(NumericDomain):
    """Real interval arithmetic with rigorous enclosures."""
    comptime Scalar = Ball
    """The apn_mojo real enclosure."""
    comptime Context = BallContext
    """Explicit midpoint and working precision."""


struct ComplexBallDomain(NumericDomain):
    """Complex rectangular interval arithmetic with rigorous enclosures."""
    comptime Scalar = ComplexBall
    """The apn_mojo rectangular enclosure."""
    comptime Context = BallContext
    """Explicit midpoint and working precision."""


def _float_cost(value: Float) raises PolyError -> Int:
    if not value.is_finite():
        raise PolyError(PolyErrorKind.DOMAIN, "numerical coefficients must be finite")
    var exponent: Int
    try:
        exponent = value.exponent() if not value.is_zero() else 0
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    if exponent == Int.MIN:
        raise PolyError(PolyErrorKind.CAPACITY, "exponent magnitude exceeds capacity")
    return _sum(value.precision(), abs(exponent))


def _numeric_cost[D: NumericDomain](value: D.Scalar) raises PolyError -> Int:
    comptime if D == FloatDomain:
        return _float_cost(rebind[Float](value))
    elif D == ComplexDomain:
        var z = rebind[Complex](value)
        return _sum(_float_cost(z.real()), _float_cost(z.imag()))
    elif D == BallDomain:
        var b = rebind[Ball](value)
        if not b.is_finite():
            raise PolyError(PolyErrorKind.DOMAIN, "numerical balls must be bounded")
        var radius: Float
        try:
            radius = b.radius()
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        return _sum(_float_cost(b.midpoint()), _float_cost(radius))
    elif D == ComplexBallDomain:
        var z = rebind[ComplexBall](value)
        return _sum(_numeric_cost[BallDomain](z.real()), _numeric_cost[BallDomain](z.imag()))
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "unknown numerical domain")


def _numeric_zero[D: NumericDomain](value: D.Scalar) -> Bool:
    comptime if D == FloatDomain:
        return rebind[Float](value).is_zero()
    elif D == ComplexDomain:
        return rebind[Complex](value).is_zero()
    elif D == BallDomain:
        var b = rebind[Ball](value)
        return b.is_exact() and b.midpoint().is_zero()
    elif D == ComplexBallDomain:
        var z = rebind[ComplexBall](value)
        return _numeric_zero[BallDomain](z.real()) and _numeric_zero[BallDomain](z.imag())
    else:
        return False


def _numeric_nonzero[D: NumericDomain](value: D.Scalar) raises PolyError -> Bool:
    _ = _numeric_cost[D](value)
    try:
        comptime if D == BallDomain:
            return not contains_zero(rebind[Ball](value))
        elif D == ComplexBallDomain:
            var z = rebind[ComplexBall](value)
            return not contains_zero(z.real()) or not contains_zero(z.imag())
        else:
            return not _numeric_zero[D](value)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _context_cost[D: NumericDomain](context: D.Context) raises PolyError -> Int:
    comptime if D == FloatDomain:
        return rebind[ArithmeticContext](context).format().precision()
    elif D == ComplexDomain:
        var c = rebind[ComplexContext](context)
        return _sum(c.real().format().precision(), c.imag().format().precision())
    elif D == BallDomain or D == ComplexBallDomain:
        var bits = rebind[BallContext](context).precision()
        if not bits:
            raise PolyError(PolyErrorKind.DOMAIN, "numerical polynomial operations require explicit ball precision")
        return bits.value()
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "unknown numerical domain")


def _numeric_integer[D: NumericDomain](value: Int, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    _charge_size(_context_cost[D](context), budget)
    try:
        comptime if D == FloatDomain:
            return rebind_var[D.Scalar](Float(Integer(value), context=rebind[ArithmeticContext](context)))
        elif D == ComplexDomain:
            return rebind_var[D.Scalar](Complex(Integer(value), context=rebind[ComplexContext](context)))
        elif D == BallDomain:
            return rebind_var[D.Scalar](Ball(Integer(value), precision=rebind[BallContext](context).precision().value()))
        elif D == ComplexBallDomain:
            return rebind_var[D.Scalar](ComplexBall(Ball(Integer(value), precision=rebind[BallContext](context).precision().value()), Ball(0, precision=rebind[BallContext](context).precision().value())))
        else:
            raise Error("unknown numerical domain")
    except e:
        raise PolyError(PolyErrorKind.DOMAIN, String(e))


def _numeric_apply[D: NumericDomain](a: D.Scalar, b: D.Scalar, operation: Int, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    var size = _sum(_sum(_numeric_cost[D](a), _numeric_cost[D](b)), _context_cost[D](context))
    _charge_size(_product(size, size) if operation >= 2 else size, budget)
    if operation == 3 and not _numeric_nonzero[D](b):
        raise PolyError(PolyErrorKind.DOMAIN, "scalar divisor is not certified nonzero")
    var result: D.Scalar
    try:
        comptime if D == FloatDomain:
            var x = rebind[Float](a)
            var y = rebind[Float](b)
            var c = rebind[ArithmeticContext](context)
            if operation == 0:
                result = rebind_var[D.Scalar](add(x, y, context=c))
            elif operation == 1:
                result = rebind_var[D.Scalar](subtract(x, y, context=c))
            elif operation == 2:
                result = rebind_var[D.Scalar](multiply(x, y, context=c))
            else:
                result = rebind_var[D.Scalar](divide(x, y, context=c))
        elif D == ComplexDomain:
            var x = rebind[Complex](a)
            var y = rebind[Complex](b)
            var c = rebind[ComplexContext](context)
            if operation == 0:
                result = rebind_var[D.Scalar](add(x, y, context=c))
            elif operation == 1:
                result = rebind_var[D.Scalar](subtract(x, y, context=c))
            elif operation == 2:
                result = rebind_var[D.Scalar](multiply(x, y, context=c))
            else:
                result = rebind_var[D.Scalar](divide(x, y, context=c))
        elif D == BallDomain:
            var x = rebind[Ball](a)
            var y = rebind[Ball](b)
            var c = rebind[BallContext](context)
            if operation == 0:
                result = rebind_var[D.Scalar](add(x, y, context=c))
            elif operation == 1:
                result = rebind_var[D.Scalar](subtract(x, y, context=c))
            elif operation == 2:
                result = rebind_var[D.Scalar](multiply(x, y, context=c))
            else:
                result = rebind_var[D.Scalar](divide(x, y, context=c))
        elif D == ComplexBallDomain:
            var x = rebind[ComplexBall](a)
            var y = rebind[ComplexBall](b)
            var c = rebind[BallContext](context)
            if operation == 0:
                result = rebind_var[D.Scalar](add(x, y, context=c))
            elif operation == 1:
                result = rebind_var[D.Scalar](subtract(x, y, context=c))
            elif operation == 2:
                result = rebind_var[D.Scalar](multiply(x, y, context=c))
            else:
                result = rebind_var[D.Scalar](divide(x, y, context=c))
        else:
            raise Error("unknown numerical domain")
    except e:
        raise PolyError(PolyErrorKind.DOMAIN, String(e))
    _ = _numeric_cost[D](result)
    return result^


def _numeric_power[D: NumericDomain](value: D.Scalar, exponent: UInt64, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    var result = _numeric_integer[D](1, context, budget)
    var base = value.copy()
    var n = exponent
    while n:
        budget.charge(1)
        if n & 1:
            result = _numeric_apply[D](result, base, 2, context, budget)
        n >>= 1
        if n:
            base = _numeric_apply[D](base, base, 2, context, budget)
    return result^
