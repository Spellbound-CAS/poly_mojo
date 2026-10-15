"""Budgeted canonical scalar text and bounded unsigned metadata parsing."""
from apn_mojo import Integer, Rational, ExactComplex, ConversionLimits
from .coefficients import CoefficientRing, ZZ, QQ, ModInt, Modulus
from .exact_complex import QQi
from .errors import Budget, PolyError, PolyErrorKind, _charge_size, _sum, _product
from .primes import _residue


def _text_work(length: Int, mut budget: Budget) raises PolyError -> Int:
    var size = _sum(length, 1)
    var work = _product(size, size)
    _charge_size(work, budget)
    return _product(work, 64)


def _unsigned(text: String, maximum: UInt64, mut budget: Budget) raises PolyError -> UInt64:
    _charge_size(text.byte_length(), budget)
    if not text.byte_length():
        raise PolyError(PolyErrorKind.ENCODING, "empty unsigned integer")
    var result = UInt64(0)
    for byte in text.as_bytes():
        if byte < 48 or byte > 57:
            raise PolyError(PolyErrorKind.ENCODING, "invalid unsigned decimal integer")
        var digit = UInt64(byte - 48)
        if digit > maximum or result > (maximum - digit) // 10:
            raise PolyError(PolyErrorKind.CAPACITY, "unsigned integer exceeds allowed range")
        result = result * 10 + digit
    return result


def _integer_text(text: String, mut budget: Budget) raises PolyError -> Integer:
    var capacity = _text_work(text.byte_length(), budget)
    try:
        return Integer.parse(text, limits=ConversionLimits(max_input_bytes=text.byte_length(), max_digits=text.byte_length(), max_allocated_bytes=capacity))
    except e:
        raise PolyError(PolyErrorKind.ENCODING, String(e))


def _rational_text(text: String, mut budget: Budget) raises PolyError -> Rational:
    var capacity = _text_work(text.byte_length(), budget)
    for byte in text.as_bytes():
        if (byte < 48 or byte > 57) and byte != 45 and byte != 47:
            raise PolyError(PolyErrorKind.ENCODING, "expected exact integer or rational text")
    try:
        return Rational.parse(text, limits=ConversionLimits(max_input_bytes=text.byte_length(), max_digits=text.byte_length(), max_allocated_bytes=capacity))
    except e:
        raise PolyError(PolyErrorKind.ENCODING, String(e))


def _literal[C: CoefficientRing](text: String, ring: C.Context, mut budget: Budget) raises PolyError -> C:
    var integer = _integer_text(text, budget)
    comptime if C == ZZ:
        return rebind_var[C](ZZ(integer^))
    elif C == QQ:
        return rebind_var[C](QQ(Rational(integer)))
    elif C == QQi:
        return rebind_var[C](QQi(ExactComplex(Rational(integer), Rational(0))))
    elif C == ModInt:
        return rebind_var[C](ModInt(_residue(integer, rebind[Modulus](ring).value, budget), rebind[Modulus](ring)))
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "text literals require a built-in exact coefficient domain")


def _scalar_text[C: CoefficientRing](value: C, mut budget: Budget) raises PolyError -> String:
    _ = _text_work(value.size(), budget)
    try:
        comptime if C == ZZ:
            return rebind[ZZ](value).value.to_string()
        elif C == QQ:
            return rebind[QQ](value).value.to_string()
        elif C == QQi:
            return String(rebind[QQi](value).value.real().to_string(), ",", rebind[QQi](value).value.imag().to_string())
        elif C == ModInt:
            return String(rebind[ModInt](value).value)
        else:
            raise Error("scalar text requires a built-in exact coefficient domain")
    except e:
        raise PolyError(PolyErrorKind.ENCODING, String(e))


def _scalar_parse[C: CoefficientRing](text: String, ring: C.Context, mut budget: Budget) raises PolyError -> C:
    comptime if C == ZZ:
        return rebind_var[C](ZZ(_integer_text(text, budget)))
    elif C == QQ:
        return rebind_var[C](QQ(_rational_text(text, budget)))
    elif C == QQi:
        _charge_size(text.byte_length(), budget)
        var parts = text.split(",")
        if len(parts) != 2:
            raise PolyError(PolyErrorKind.ENCODING, "Gaussian rational text needs two rational components")
        return rebind_var[C](QQi(ExactComplex(_rational_text(String(parts[0]), budget), _rational_text(String(parts[1]), budget))))
    elif C == ModInt:
        return rebind_var[C](ModInt(_unsigned(text, rebind[Modulus](ring).value - 1, budget), rebind[Modulus](ring)))
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "scalar text requires a built-in exact coefficient domain")


def _domain[C: CoefficientRing](ring: C.Context) raises PolyError -> String:
    comptime if C == ZZ:
        return "ZZ"
    elif C == QQ:
        return "QQ"
    elif C == QQi:
        return "QQi"
    elif C == ModInt:
        return String("Fp:", rebind[Modulus](ring).value)
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "interchange requires a built-in exact coefficient domain")
