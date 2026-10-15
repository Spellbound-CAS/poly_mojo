"""Canonical versioned text interchange for built-in exact coefficient domains."""
from .coefficients import CoefficientRing
from .dense import DensePoly
from .sparse import SparsePoly, Term
from .monomial import Monomial, MonomialOrder
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from ._text_ops import _scalar_text, _scalar_parse, _domain
from ._interchange import _read_header, _write_header, _TextTerm


def serialize_dense[C: CoefficientRing](poly: DensePoly[C], mut budget: Budget) raises PolyError -> String:
    """Serialize exact dense coefficients and their variable and domain.

    Parameters:
        C: ZZ, QQ, QQi or ModInt.

    Args:
        poly: The canonical dense polynomial.
        budget: The shared work allowance.

    Returns:
        Canonical version-one polynomial text, including a final newline.

    Raises:
        PolyError: DOMAIN for invalid input, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    var writer = _write_header(True, _domain[C](poly.ring), poly.variable, MonomialOrder(), len(poly.coefficients), budget)
    for coefficient_index in range(len(poly.coefficients)):
        ref value = poly.coefficients[coefficient_index]
        writer.add(_scalar_text(value, budget), budget)
    return writer.finish(budget)


@always_inline
def _coefficient[C: CoefficientRing](text: String, ring: C.Context, mut budget: Budget) raises PolyError -> C:
    var value = _scalar_parse[C](text, ring, budget)
    var canonical = _scalar_text(value, budget)
    _charge_size(_sum(text.byte_length(), canonical.byte_length()), budget)
    if text != canonical:
        raise PolyError(PolyErrorKind.ENCODING, "noncanonical coefficient text")
    return value^


def deserialize_dense[C: CoefficientRing](text: String, ring: C.Context, mut budget: Budget) raises PolyError -> DensePoly[C]:
    """Decode canonical exact dense interchange without repairing malformed data.

    Parameters:
        C: The expected built-in exact coefficient domain.

    Args:
        text: Version-one polynomial text.
        ring: The expected context, including the prime modulus.
        budget: The shared work allowance.

    Returns:
        The exact dense polynomial with its original variable.

    Raises:
        PolyError: ENCODING for malformed or incompatible text, LIMIT or CAPACITY on work limits.
    """
    var record = _read_header(text, budget)
    _charge_size(record.domain.byte_length(), budget)
    if not record.dense or record.domain != _domain[C](ring):
        raise PolyError(PolyErrorKind.ENCODING, "expected dense polynomial in requested ring")
    _allocation[C](record.count, budget)
    var result = DensePoly[C](record.variable, ring.copy())
    result.coefficients = List[C](capacity=record.count)
    for _ in range(record.count):
        result.coefficients.append(_coefficient[C](record.reader.line(budget), ring, budget))
    record.finish()
    if len(result.coefficients):
        _charge_size(result.coefficients[len(result.coefficients) - 1].size(), budget)
        if result.coefficients[len(result.coefficients) - 1].is_zero():
            raise PolyError(PolyErrorKind.ENCODING, "trailing zero dense coefficient")
    return result^


def serialize_sparse[C: CoefficientRing](poly: SparsePoly[C], mut budget: Budget) raises PolyError -> String:
    """Serialize exact sparse terms, their order and coefficient domain.

    Parameters:
        C: ZZ, QQ, QQi or ModInt.

    Args:
        poly: The canonical sparse polynomial.
        budget: The shared work allowance.

    Returns:
        Canonical version-one polynomial text.

    Raises:
        PolyError: DOMAIN for invalid input, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
    """
    poly.check_invariants(budget)
    var writer = _write_header(False, _domain[C](poly.ring), 0, poly.order, len(poly.terms), budget)
    for term in poly.terms:
        writer.monomial(term.monomial, budget)
        writer.add(_scalar_text(term.coefficient, budget), budget)
    return writer.finish(budget)


def deserialize_sparse[C: CoefficientRing](text: String, ring: C.Context, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Decode exact sparse interchange while enforcing canonical term order.

    Parameters:
        C: The expected built-in exact coefficient domain.

    Args:
        text: Version-one polynomial text.
        ring: The expected coefficient context.
        budget: The shared work allowance.

    Returns:
        The exact sparse polynomial, retaining its monomial order.

    Raises:
        PolyError: ENCODING for malformed or incompatible text, LIMIT or CAPACITY on work limits.
    """
    var record = _read_header(text, budget)
    _charge_size(record.domain.byte_length(), budget)
    if record.dense or record.domain != _domain[C](ring):
        raise PolyError(PolyErrorKind.ENCODING, "expected sparse polynomial in requested ring")
    _allocation[Term[C]](record.count, budget)
    _charge_size(len(record.order.variables), budget)
    var result = SparsePoly[C](record.order.copy(), ring.copy())
    result.terms = List[Term[C]](capacity=record.count)
    var position = 0
    while position < record.count:
        # Keep metadata scanning and scalar conversion in short contiguous
        # passes, without retaining an input-sized intermediate record.
        var count = min(64, record.count - position)
        _allocation[_TextTerm](count, budget)
        var block = List[_TextTerm](capacity=count)
        for i in range(count):
            var monomial = record.monomial(budget)
            if i:
                record.increasing(block[i - 1].monomial, monomial, budget)
            elif position:
                record.increasing(result.terms[position - 1].monomial, monomial, budget)
            block.append(_TextTerm(monomial^, record.reader.line(budget)))
        for i in range(count):
            var coefficient = _coefficient[C](block[i].coefficient, ring, budget)
            _charge_size(coefficient.size(), budget)
            if coefficient.is_zero():
                raise PolyError(PolyErrorKind.ENCODING, "zero sparse coefficient")
            var monomial = Monomial()
            swap(monomial, block[i].monomial)
            result.terms.append(Term[C](monomial^, coefficient^))
        position += count
    record.finish()
    return result^

from .text_stream import PolynomialTextSource, PolynomialTextSink, _chunk_header, _flush


def serialize_dense_to[C: CoefficientRing, S: PolynomialTextSink](poly: DensePoly[C], mut sink: S, mut budget: Budget, *, chunk_bytes: Int = 4096) raises PolyError:
    """Write canonical polynomial text in bounded batches.

    Parameters:
        C: The coefficient domain.
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
    poly.check_invariants(budget)
    var writer = _write_header(True, _domain[C](poly.ring), poly.variable, MonomialOrder(), len(poly.coefficients), budget)
    _flush(writer, sink, chunk_bytes, budget)
    for i in range(len(poly.coefficients)):
        ref value = poly.coefficients[i]
        writer.add(_scalar_text(value, budget), budget)
        _flush(writer, sink, chunk_bytes, budget)
    _flush(writer, sink, chunk_bytes, budget, final=True)


def deserialize_dense_from[C: CoefficientRing, S: PolynomialTextSource](var source: S, ring: C.Context, mut budget: Budget, *, max_field_bytes: Int = 16777216) raises PolyError -> DensePoly[C]:
    """Decode one canonical polynomial from arbitrary text chunks.

    Parameters:
        C: The expected coefficient domain.
        S: The owned text source type.

    Args:
        source: The source consumed through end of input.
        ring: The expected exact coefficient context.
        budget: The shared parsing and transport allowance.
        max_field_bytes: Positive limit for any individual metadata or scalar field.

    Returns:
        The complete polynomial, only after canonical validation and end of input.

    Raises:
        PolyError: ENCODING for malformed input, DOMAIN for invalid limits or nonfinite scalars, LIMIT or CAPACITY on work limits.
    """
    var record = _chunk_header(source^, True, _domain[C](ring), max_field_bytes, budget)
    var values = List[C]()
    for _ in range(record.count):
        var value = _coefficient[C](record.reader.line(budget), ring, budget)
        _allocation[C](1, budget)
        values.append(value^)
    record.reader.finish(budget)
    if len(values) and values[len(values) - 1].is_zero():
        raise PolyError(PolyErrorKind.ENCODING, "trailing zero dense coefficient")
    return DensePoly[C].from_coefficients(record.variable, ring, values^, budget)


def serialize_sparse_to[C: CoefficientRing, S: PolynomialTextSink](poly: SparsePoly[C], mut sink: S, mut budget: Budget, *, chunk_bytes: Int = 4096) raises PolyError:
    """Write canonical polynomial text in bounded batches.

    Parameters:
        C: The coefficient domain.
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
    poly.check_invariants(budget)
    var writer = _write_header(False, _domain[C](poly.ring), 0, poly.order, len(poly.terms), budget)
    _flush(writer, sink, chunk_bytes, budget)
    for i in range(len(poly.terms)):
        ref value = poly.terms[i]
        writer.monomial(value.monomial, budget)
        writer.add(_scalar_text(value.coefficient, budget), budget)
        _flush(writer, sink, chunk_bytes, budget)
    _flush(writer, sink, chunk_bytes, budget, final=True)


def deserialize_sparse_from[C: CoefficientRing, S: PolynomialTextSource](var source: S, ring: C.Context, mut budget: Budget, *, max_field_bytes: Int = 16777216) raises PolyError -> SparsePoly[C]:
    """Decode one canonical polynomial from arbitrary text chunks.

    Parameters:
        C: The expected coefficient domain.
        S: The owned text source type.

    Args:
        source: The source consumed through end of input.
        ring: The expected exact coefficient context.
        budget: The shared parsing and transport allowance.
        max_field_bytes: Positive limit for any individual metadata or scalar field.

    Returns:
        The complete polynomial, only after canonical validation and end of input.

    Raises:
        PolyError: ENCODING for malformed input, DOMAIN for invalid limits or nonfinite scalars, LIMIT or CAPACITY on work limits.
    """
    var record = _chunk_header(source^, False, _domain[C](ring), max_field_bytes, budget)
    var result = SparsePoly[C](record.order.copy(), ring.copy())
    for i in range(record.count):
        var monomial = record.monomial(budget)
        if i and record.order.compare(result.terms[i - 1].monomial, monomial, budget) >= 0:
            raise PolyError(PolyErrorKind.ENCODING, "nonincreasing serialized terms")
        var value = _coefficient[C](record.reader.line(budget), ring, budget)
        if value.is_zero():
            raise PolyError(PolyErrorKind.ENCODING, "zero sparse coefficient")
        _allocation[Term[C]](1, budget)
        result.terms.append(Term[C](monomial^, value^))
    record.reader.finish(budget)
    return result^
