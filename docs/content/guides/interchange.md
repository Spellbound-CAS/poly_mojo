# Save, format, and stream polynomial data

Choose an interchange format according to its reader. Named expressions are
convenient for people; canonical text preserves supported polynomial data;
the kernel codec serves a specific proof-system boundary. These formats are
not interchangeable.

## Text for people and programs

Use a `PolynomialRing` to parse and format expressions with stable variable
names. Its grammar covers polynomial expressions with explicit multiplication,
nonnegative integer powers, and supported scalar fractions. It is not a
symbolic expression parser for functions such as `sin` or `exp`.

Use `serialize_dense` and `serialize_sparse`, with their matching deserializers,
for canonical exact polynomial text. Numerical containers have separate text
entry points. The format records the domain and representation needed by the
reader; malformed or noncanonical input is rejected. Diagnostic `Writable`
output is not a serialization contract.

When text is large, the `*_to` writers emit chunks to a `PolynomialTextSink`.
The corresponding reader APIs consume bytes from a `PolynomialTextSource`.
A sink can write to storage or process each chunk without collecting a complete
string. The [toolkit reference](../reference/toolkit.md#incremental-text-and-cad-transport)
defines the protocols and the byte-for-byte compatibility with materialized
text.

## The kernel codec

`to_kernel_bytes` and `from_kernel_bytes` use the kernel's `FORMAL_POLYNOMIAL` encoding
for exact rational polynomials. This is a bounded binary format with canonical
rational coefficients and fixed count/exponent encodings. Zero and one use
dedicated kernel constructors and are rejected by the polynomial codec.

The codec does not encode `QQi`, numerical point or ball polynomials, or general
certificate objects. See [Contract P](../reference/contract-p.md) for the
integration and compatibility boundary.

## Consume proof data incrementally

<!-- example: docs/examples/streaming.mojo -->

The first part counts canonical text bytes as they arrive. The remaining
parts send real-root, complex-root, and CAD evidence directly into independent
stream checkers. Each producer returns a receipt; `finish` verifies that the
received evidence establishes the complete claim.

The lifecycle is the same for each proof stream:

1. Prepare a checker from the original inputs and options.
2. Pass it as the producer's sink.
3. Let it validate the ordered events as they arrive.
4. Call `finish` with the returned receipt before accepting the result.

For enclosure evaluation, the incremental APIs similarly allow witness steps
to be consumed without retaining the completed materialized trace. See the
[numerical reference](../spec/numerical.md#incremental-evaluation-witnesses).

Streaming reduces retained intermediate payloads; it does not promise
constant memory. Producers and checkers still need input polynomials, graphs,
live values, search state, and the coverage information required by the claim.
The [streaming architecture](../architecture/streaming.md) explains those
lifetimes and the tradeoffs with materialized results.
