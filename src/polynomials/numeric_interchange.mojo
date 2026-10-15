"""Lossless numerical polynomial interchange using public apn_mojo scalar JSON."""
from apn_mojo import Float, Complex, Ball, ComplexBall, ConversionLimits
from .numeric_domains import NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, _numeric_cost, _numeric_zero
from .numeric_dense import NumericDensePoly
from .numeric_sparse import NumericSparsePoly, NumericTerm
from .monomial import MonomialOrder
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from ._text_ops import _text_work
from ._interchange import _read_header, _write_header


def _domain[D: NumericDomain]() raises PolyError -> String:
    comptime if D == FloatDomain:
        return "Float"
    elif D == ComplexDomain:
        return "Complex"
    elif D == BallDomain:
        return "Ball"
    elif D == ComplexBallDomain:
        return "ComplexBall"
    else:
        raise PolyError(PolyErrorKind.UNSUPPORTED, "unknown numerical interchange domain")


def _scalar_text[D: NumericDomain](value: D.Scalar, mut budget: Budget) raises PolyError -> String:
    _ = _text_work(_numeric_cost[D](value), budget)
    try:
        comptime if D == FloatDomain:
            return rebind[Float](value).to_json()
        elif D == ComplexDomain:
            return rebind[Complex](value).to_json()
        elif D == BallDomain:
            return rebind[Ball](value).to_json()
        elif D == ComplexBallDomain:
            return rebind[ComplexBall](value).to_json()
        else:
            raise Error("unknown numerical interchange domain")
    except e:
        raise PolyError(PolyErrorKind.ENCODING, String(e))


def _scalar_parse[D: NumericDomain](text: String, mut budget: Budget) raises PolyError -> D.Scalar:
    var capacity = _text_work(text.byte_length(), budget)
    try:
        var limits = ConversionLimits(max_input_bytes=text.byte_length(), max_digits=text.byte_length(), max_allocated_bytes=capacity)
        comptime if D == FloatDomain:
            return rebind_var[D.Scalar](Float.from_json(text, limits=limits))
        elif D == ComplexDomain:
            return rebind_var[D.Scalar](Complex.from_json(text, limits=limits))
        elif D == BallDomain:
            return rebind_var[D.Scalar](Ball.from_json(text, limits=limits))
        elif D == ComplexBallDomain:
            return rebind_var[D.Scalar](ComplexBall.from_json(text, limits=limits))
        else:
            raise Error("unknown numerical interchange domain")
    except e:
        raise PolyError(PolyErrorKind.ENCODING, String(e))


@always_inline
def _coefficient[D: NumericDomain](text: String, mut budget: Budget) raises PolyError -> D.Scalar:
    var value = _scalar_parse[D](text, budget)
    _charge_size(_numeric_cost[D](value), budget)
    var canonical = _scalar_text[D](value, budget)
    _charge_size(_sum(text.byte_length(), canonical.byte_length()), budget)
    if canonical != text:
        raise PolyError(PolyErrorKind.ENCODING, "noncanonical numerical coefficient")
    return value^


def serialize_numeric_dense[D: NumericDomain](poly: NumericDensePoly[D], mut budget: Budget) raises PolyError -> String:
    """Serialize numerical dense coefficients without rounding them.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The canonical numerical polynomial.
        budget: The shared work allowance.

    Returns:
        Version-one text preserving scalar formats and ball enclosures exactly.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    if len(poly.coefficients) and _numeric_zero[D](poly.coefficients[len(poly.coefficients) - 1]):
        raise PolyError(PolyErrorKind.DOMAIN, "noncanonical numerical dense polynomial")
    var writer = _write_header(True, _domain[D](), poly.variable, MonomialOrder(), len(poly.coefficients), budget)
    for value in poly.coefficients:
        writer.add(_scalar_text[D](value, budget), budget)
    return writer.finish(budget)


def deserialize_numeric_dense[D: NumericDomain](text: String, mut budget: Budget) raises PolyError -> NumericDensePoly[D]:
    """Read numerical dense coefficients exactly from versioned interchange.

    Parameters:
        D: The expected numerical domain.

    Args:
        text: Canonical version-one polynomial text.
        budget: The shared work allowance.

    Returns:
        The original scalar values, formats and retained polynomial variable.

    Raises:
        PolyError: ENCODING for malformed text, DOMAIN for nonfinite coefficients, LIMIT or CAPACITY on work limits.
    """
    var record = _read_header(text, budget)
    _charge_size(record.domain.byte_length(), budget)
    if not record.dense or record.domain != _domain[D]():
        raise PolyError(PolyErrorKind.ENCODING, "expected numerical dense interchange in requested domain")
    _allocation[D.Scalar](record.count, budget)
    var coefficients = List[D.Scalar](capacity=record.count)
    for _ in range(record.count):
        coefficients.append(_coefficient[D](record.reader.line(budget), budget))
    record.finish()
    if len(coefficients) and _numeric_zero[D](coefficients[len(coefficients) - 1]):
        raise PolyError(PolyErrorKind.ENCODING, "trailing certain-zero numerical coefficient")
    return NumericDensePoly[D].from_coefficients(record.variable, coefficients^, budget)


def serialize_numeric_sparse[D: NumericDomain](poly: NumericSparsePoly[D], mut budget: Budget) raises PolyError -> String:
    """Serialize sparse numerical coefficients and monomial-order metadata.

    Parameters:
        D: The numerical domain.

    Args:
        poly: The canonical numerical sparse polynomial.
        budget: The shared work allowance.

    Returns:
        Version-one text preserving every scalar's format and enclosure.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    var writer = _write_header(False, _domain[D](), 0, poly.order, len(poly.terms), budget)
    for i in range(len(poly.terms)):
        ref term = poly.terms[i]
        term.monomial.check_invariants(budget)
        poly.order.validate(term.monomial, budget)
        if i and poly.order.compare(poly.terms[i - 1].monomial, term.monomial, budget) >= 0:
            raise PolyError(PolyErrorKind.DOMAIN, "noncanonical numerical term order")
        if _numeric_zero[D](term.coefficient):
            raise PolyError(PolyErrorKind.DOMAIN, "zero numerical sparse term")
        writer.monomial(term.monomial, budget)
        writer.add(_scalar_text[D](term.coefficient, budget), budget)
    return writer.finish(budget)


def deserialize_numeric_sparse[D: NumericDomain](text: String, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    """Read canonical numerical sparse interchange without arithmetic or rounding.

    Parameters:
        D: The expected numerical domain.

    Args:
        text: Canonical version-one polynomial text.
        budget: The shared work allowance.

    Returns:
        The exact stored scalar formats, enclosures and monomial order.

    Raises:
        PolyError: ENCODING for malformed text, DOMAIN for nonfinite coefficients, LIMIT or CAPACITY on work limits.
    """
    var record = _read_header(text, budget)
    _charge_size(record.domain.byte_length(), budget)
    if record.dense or record.domain != _domain[D]():
        raise PolyError(PolyErrorKind.ENCODING, "expected numerical sparse interchange in requested domain")
    _allocation[NumericTerm[D]](record.count, budget)
    _charge_size(len(record.order.variables), budget)
    var result = NumericSparsePoly[D](record.order.copy())
    result.terms = List[NumericTerm[D]](capacity=record.count)
    for i in range(record.count):
        var monomial = record.monomial(budget)
        if i:
            record.increasing(result.terms[i - 1].monomial, monomial, budget)
        var coefficient = _coefficient[D](record.reader.line(budget), budget)
        if _numeric_zero[D](coefficient):
            raise PolyError(PolyErrorKind.ENCODING, "zero sparse numerical coefficient")
        result.terms.append(NumericTerm[D](monomial^, coefficient^))
    record.finish()
    return result^

from .text_stream import PolynomialTextSource, PolynomialTextSink, _chunk_header, _flush


def serialize_numeric_dense_to[D: NumericDomain, S: PolynomialTextSink](poly: NumericDensePoly[D], mut sink: S, mut budget: Budget, *, chunk_bytes: Int = 4096) raises PolyError:
    """Write canonical polynomial text in bounded batches.

    Parameters:
        D: The coefficient domain.
        S: The text destination type.

    Args:
        poly: The canonical polynomial to serialize.
        sink: The destination receiving ordered chunks.
        budget: The shared work and transport allowance.
        chunk_bytes: Positive target size; a header or single term can exceed it.

    Raises:
        PolyError: DOMAIN for invalid inputs, ENCODING for transport failure, LIMIT or CAPACITY on work limits; a failed write can leave an incomplete destination.
    """
    if chunk_bytes < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "chunk size must be positive")
    if len(poly.coefficients) and _numeric_zero[D](poly.coefficients[len(poly.coefficients) - 1]):
        raise PolyError(PolyErrorKind.DOMAIN, "noncanonical numerical polynomial")
    var writer = _write_header(True, _domain[D](), poly.variable, MonomialOrder(), len(poly.coefficients), budget)
    _flush(writer, sink, chunk_bytes, budget)
    for i in range(len(poly.coefficients)):
        ref value = poly.coefficients[i]
        writer.add(_scalar_text[D](value, budget), budget)
        _flush(writer, sink, chunk_bytes, budget)
    _flush(writer, sink, chunk_bytes, budget, final=True)


def deserialize_numeric_dense_from[D: NumericDomain, S: PolynomialTextSource](var source: S, mut budget: Budget, *, max_field_bytes: Int = 16777216) raises PolyError -> NumericDensePoly[D]:
    """Decode one canonical polynomial from arbitrary text chunks.

    Parameters:
        D: The expected coefficient domain.
        S: The owned text source type.

    Args:
        source: The source consumed through end of input.
        budget: The shared parsing and transport allowance.
        max_field_bytes: Positive limit for any individual metadata or scalar field.

    Returns:
        The complete polynomial, only after canonical validation and end of input.

    Raises:
        PolyError: ENCODING for malformed input, DOMAIN for invalid limits or nonfinite scalars, LIMIT or CAPACITY on work limits.
    """
    var record = _chunk_header(source^, True, _domain[D](), max_field_bytes, budget)
    var values = List[D.Scalar]()
    for _ in range(record.count):
        var value = _coefficient[D](record.reader.line(budget), budget)
        _allocation[D.Scalar](1, budget)
        values.append(value^)
    record.reader.finish(budget)
    if len(values) and _numeric_zero[D](values[len(values) - 1]):
        raise PolyError(PolyErrorKind.ENCODING, "trailing zero dense coefficient")
    return NumericDensePoly[D].from_coefficients(record.variable, values^, budget)


def serialize_numeric_sparse_to[D: NumericDomain, S: PolynomialTextSink](poly: NumericSparsePoly[D], mut sink: S, mut budget: Budget, *, chunk_bytes: Int = 4096) raises PolyError:
    """Write canonical polynomial text in bounded batches.

    Parameters:
        D: The coefficient domain.
        S: The text destination type.

    Args:
        poly: The canonical polynomial to serialize.
        sink: The destination receiving ordered chunks.
        budget: The shared work and transport allowance.
        chunk_bytes: Positive target size; a header or single term can exceed it.

    Raises:
        PolyError: DOMAIN for invalid inputs, ENCODING for transport failure, LIMIT or CAPACITY on work limits; a failed write can leave an incomplete destination.
    """
    if chunk_bytes < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "chunk size must be positive")
    var writer = _write_header(False, _domain[D](), 0, poly.order, len(poly.terms), budget)
    _flush(writer, sink, chunk_bytes, budget)
    for i in range(len(poly.terms)):
        ref value = poly.terms[i]
        value.monomial.check_invariants(budget)
        poly.order.validate(value.monomial, budget)
        if _numeric_zero[D](value.coefficient) or (i and poly.order.compare(poly.terms[i - 1].monomial, value.monomial, budget) >= 0):
            raise PolyError(PolyErrorKind.DOMAIN, "noncanonical numerical sparse term")
        writer.monomial(value.monomial, budget)
        writer.add(_scalar_text[D](value.coefficient, budget), budget)
        _flush(writer, sink, chunk_bytes, budget)
    _flush(writer, sink, chunk_bytes, budget, final=True)


def deserialize_numeric_sparse_from[D: NumericDomain, S: PolynomialTextSource](var source: S, mut budget: Budget, *, max_field_bytes: Int = 16777216) raises PolyError -> NumericSparsePoly[D]:
    """Decode one canonical polynomial from arbitrary text chunks.

    Parameters:
        D: The expected coefficient domain.
        S: The owned text source type.

    Args:
        source: The source consumed through end of input.
        budget: The shared parsing and transport allowance.
        max_field_bytes: Positive limit for any individual metadata or scalar field.

    Returns:
        The complete polynomial, only after canonical validation and end of input.

    Raises:
        PolyError: ENCODING for malformed input, DOMAIN for invalid limits or nonfinite scalars, LIMIT or CAPACITY on work limits.
    """
    var record = _chunk_header(source^, False, _domain[D](), max_field_bytes, budget)
    var result = NumericSparsePoly[D](record.order.copy())
    for i in range(record.count):
        var monomial = record.monomial(budget)
        if i and record.order.compare(result.terms[i - 1].monomial, monomial, budget) >= 0:
            raise PolyError(PolyErrorKind.ENCODING, "nonincreasing serialized terms")
        var value = _coefficient[D](record.reader.line(budget), budget)
        if _numeric_zero[D](value):
            raise PolyError(PolyErrorKind.ENCODING, "zero sparse coefficient")
        _allocation[NumericTerm[D]](1, budget)
        result.terms.append(NumericTerm[D](monomial^, value^))
    record.reader.finish(budget)
    return result^
