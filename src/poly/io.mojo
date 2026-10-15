"""Save and restore complete polynomial values without losing their metadata.

Use dumps/loads for strings and dump/load for incremental text transports.
The expected result type is explicit: io.loads[type_of(value)](text). Records
preserve exact witnesses and stored numerical formats; subsequent arithmetic
still takes an explicit precision context.
"""
from polynomials import Budget, PolyError, PolyErrorKind, PolynomialTextSource, PolynomialTextSink
from ._algebraic_io import _Reader, _Writer, _Source, _Sink, _bad
from ._common import _charge_size
from ._value_io import _write, _read

def dump[T: Copyable & Deinitable, S: PolynomialTextSink](value: T, mut sink: S, mut budget: Budget) raises PolyError:
    """Stream a complete value, including its representation metadata.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.
        S: The text transport implementation.

    Args:
        value: A prepared polynomial, series, basis, rational function, supported scalar, or batch.
        sink: The text sink receiving fields in order; a failure can leave partial output.
        budget: The shared work allowance, including decoding and validation.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    var writer = _Writer()
    writer.line("POLYNOMIALS-VALUE 1", sink, budget)
    _write(value, writer, sink, budget)
    writer.flush(sink, budget)


def dump[T: Copyable & Deinitable, S: PolynomialTextSink](value: T, mut sink: S, *, max_steps: UInt64 = UInt64.MAX) raises PolyError:
    """Stream a complete value, including its representation metadata.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.
        S: The text transport implementation.

    Args:
        value: A prepared polynomial, series, basis, rational function, supported scalar, or batch.
        sink: The text sink receiving fields in order; a failure can leave partial output.
        max_steps: The complete call's work allowance.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    dump(value, sink, budget)


def dumps[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> String:
    """Serialize a complete value to versioned text.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.

    Args:
        value: A prepared polynomial, series, basis, rational function, supported scalar, or batch.
        budget: The shared work allowance, including decoding and validation.

    Returns:
        Version-one text with a final newline.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    var sink = _Sink(String())
    dump(value, sink, budget)
    var result = String()
    swap(result, sink.text)
    return result^


def dumps[T: Copyable & Deinitable](value: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> String:
    """Serialize a complete value to versioned text.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.

    Args:
        value: A prepared polynomial, series, basis, rational function, supported scalar, or batch.
        max_steps: The complete call's work allowance.

    Returns:
        Version-one text with a final newline.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return dumps(value, budget)


def load[T: Copyable & Deinitable, S: PolynomialTextSource](var source: S, mut budget: Budget, *, max_field_bytes: Int = 16777216) raises PolyError -> T:
    """Validate and restore a value from arbitrary text chunks.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.
        S: The text transport implementation.

    Args:
        source: An owned source ending immediately after this record.
        budget: The shared work allowance, including decoding and validation.
        max_field_bytes: Maximum bytes per input field; defaults to 16 MiB.

    Returns:
        An independent value with its domain, metadata, and stored precision preserved.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    if max_field_bytes < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "field limit must be positive")
    var reader = _Reader[S](source^, String(), 0, max_field_bytes)
    try:
        if reader.line(budget) != "POLYNOMIALS-VALUE 1":
            raise _bad("unsupported polynomial value format")
        var result = _read[T](reader, budget)
        reader.finish(budget)
        return result^
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            raise _bad(e.message)
        raise e^


def load[T: Copyable & Deinitable, S: PolynomialTextSource](var source: S, *, max_field_bytes: Int = 16777216, max_steps: UInt64 = UInt64.MAX) raises PolyError -> T:
    """Validate and restore a value from arbitrary text chunks.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.
        S: The text transport implementation.

    Args:
        source: An owned source ending immediately after this record.
        max_field_bytes: Maximum bytes per input field; defaults to 16 MiB.
        max_steps: The complete call's work allowance.

    Returns:
        An independent value with its domain, metadata, and stored precision preserved.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return load[T](source^, budget, max_field_bytes=max_field_bytes)


def loads[T: Copyable & Deinitable](text: String, mut budget: Budget, *, max_field_bytes: Int = 16777216) raises PolyError -> T:
    """Validate and restore a value from versioned text.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.

    Args:
        text: One complete record produced by this namespace.
        budget: The shared work allowance, including decoding and validation.
        max_field_bytes: Maximum bytes per input field; defaults to 16 MiB.

    Returns:
        An independent value with its domain, metadata, and stored precision preserved.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    _charge_size(text.byte_length(), budget)
    return load[T](_Source(Optional(text.copy())), budget, max_field_bytes=max_field_bytes)


def loads[T: Copyable & Deinitable](text: String, *, max_field_bytes: Int = 16777216, max_steps: UInt64 = UInt64.MAX) raises PolyError -> T:
    """Validate and restore a value from versioned text.

    Parameters:
        T: The concrete value type, usually inferred on write or supplied as type_of(value) on read.

    Args:
        text: One complete record produced by this namespace.
        max_field_bytes: Maximum bytes per input field; defaults to 16 MiB.
        max_steps: The complete call's work allowance.

    Returns:
        An independent value with its domain, metadata, and stored precision preserved.

    Raises:
        PolyError: ENCODING for malformed or incompatible records, DOMAIN for invalid settings or inputs, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return loads[T](text, budget, max_field_bytes=max_field_bytes)


