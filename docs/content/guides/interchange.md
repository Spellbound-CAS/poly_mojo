# Save, format, and stream polynomial data

Choose an interchange format according to its reader. Formatted text is
convenient for people; value records preserve complete polynomial values;
the kernel codec serves a specific proof-system boundary. These formats are
not interchangeable.

## Text for people and programs

Use `format` to display a polynomial with the variable names you choose. It
prints every term with its exact coefficient, but omits domain, modulus, and
shape metadata, so it is not a way to save data.

Use `io.dumps` and `io.loads` for complete values. A record keeps the
coefficient domain and the representation's metadata: series order, basis
nodes, numerator and denominator, algebraic root witnesses, numerical
precision, and batch shape. `dense.dumps` and `sparse.dumps` exchange canonical
coefficient records with the matching `loads`. Each reader checks that a record
matches the requested type; malformed or noncanonical input raises `ENCODING`.
Diagnostic `Writable` output is not a serialization contract.

When text is large, `io.dump` and the `dump` functions of `poly.dense` and
`poly.sparse` emit chunks to a `PolynomialTextSink`. The corresponding `load`
functions consume bytes from a `PolynomialTextSource`. A sink can write to
storage or process each chunk without collecting a complete string. The
[serialization tutorial](../tutorial/array-serialization.md) shows each format
and its transport limits.

## The kernel codec

`dense.to_kernel_bytes` and `sparse.to_kernel_bytes` use the kernel's
`FORMAL_POLYNOMIAL` encoding for exact rational polynomials. This is a bounded
binary format with canonical rational coefficients and fixed count/exponent
encodings. Zero and one use dedicated kernel constructors and are rejected by
the polynomial codec.

The codec does not encode `QQi`, numerical point or ball polynomials, or general
certificate objects. See [kernel literal bytes](../tutorial/array-serialization.md#exchange-kernel-literal-bytes)
for the variable map that travels with each payload.

## Consume proof data incrementally

<!-- example: docs/examples/array_solver_streams.mojo -->

The program sends real-root and CAD evidence directly into independent stream
checkers. Each producer returns a receipt; `finish` verifies that the received
evidence establishes the complete claim.

The lifecycle is the same for each proof stream:

1. Prepare a checker from the original inputs and options.
2. Pass it as the producer's sink.
3. Let it validate the ordered events as they arrive.
4. Call `finish` with the returned receipt before accepting the result.

For enclosure evaluation, the incremental APIs similarly allow witness steps
to be consumed without retaining the completed materialized trace. See
[check one witness step at a time](../tutorial/array-certified.md#check-one-witness-step-at-a-time).

Streaming reduces retained intermediate payloads; it does not promise
constant memory. Producers and checkers still need input polynomials, graphs,
live values, search state, and the coverage information required by the claim.
The [streaming architecture](../architecture/streaming.md) explains those
lifetimes and the tradeoffs with materialized results.
