# Value serialization

`poly.io` saves complete values: dense and sparse polynomials, formal series,
Newton/Chebyshev/Bernstein representations, rational functions, supported
coefficient scalars, checked results, and shaped batches of these values. Nested batches work
as well. Records retain exact domains, representation metadata, algebraic
witnesses, and stored numerical formats.

Use `io.loads[type_of(value)](text)` when restoring an existing type, or name
the concrete type directly, such as `series.Array[1]` for a rational series.
The reader validates that the record matches that type. Construction and
arithmetic keep their usual domain inference and explicit precision contexts.

All calls accept a shared `Budget` or `max_steps=`. Incremental readers also
accept `max_field_bytes=`. See [saving and restoring values](../tutorial/array-serialization.md)
for examples, transport behavior, and the difference between value records
and the existing coefficient-array records.

<!-- api: src/poly/io.mojo dumps -->

<!-- api: src/poly/io.mojo loads -->

<!-- api: src/poly/io.mojo dump -->

<!-- api: src/poly/io.mojo load -->

`dump` and `load` and the coefficient-record `dump` and `load` functions in
`poly.dense` and `poly.sparse` exchange text through these two protocols. Import
them from `polynomials` and implement them over your file, network, or other
storage.

<!-- api: src/polynomials/text_stream.mojo PolynomialTextSink -->

<!-- api: src/polynomials/text_stream.mojo PolynomialTextSource -->

Checked archives support `exact.Result`, exact real and complex root sets,
interval counts, reusable Sturm plans and their counts, accepted root covers,
multiplicities, diagnostics, successful
polishing results, specialization results, CAD results, rational partial
fractions and Laurent expansions, and polynomial or fraction enclosure results.
This includes algebraic covers with conversion witnesses, `AA` specialization
transcripts, and adaptive quotient outcomes with or without a finite enclosure.
The reader restores their evidence and runs the corresponding independent
checker. It does not rerun the producing solver. Unsupported result types fail
at compilation; records with rejected evidence fail with `ENCODING`.

Certificate archives have their own versioned `CHECKED-*` type tag inside the
value framing. CAD decompositions and full coverings now write `CHECKED-3`,
with compact final data or an explicitly tagged legacy payload. Their readers
also accept the old `CHECKED-1` records and independently check the retained
evidence without rerunning construction. Older readers reject the new tag.
Selective solutions remain `CHECKED-2`. Keep the producing package version alongside long-lived archives;
the tag is based on the Mojo type name plus an explicit axis count and is not a cross-version ABI promise.

CAD decomposition, solution, and covering archives retain any `fixed=`
coordinate bindings with their exact witnesses. Loading verifies the geometry
on that slice and rejects evidence for a different conjugate. Existing CAD
archives without user-fixed coordinates remain readable; older library
versions reject the new coordinate-scope records.
