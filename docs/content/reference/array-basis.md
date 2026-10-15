# Alternative basis arrays

Construct and convert Newton, first-kind Chebyshev, and Bernstein representations. Coefficient extent and basis metadata are retained.

Every computational call accepts a shared `Budget` or `max_steps=`. See the
[tutorial](../tutorial/array-bases.md) for examples, domain restrictions,
and the meaning of retained metadata.

Save and restore these values with [poly.io](array-io.md).

<!-- api: src/poly/basis.mojo Array -->

<!-- api: src/poly/basis.mojo asarray -->

<!-- api: src/poly/basis.mojo from_polynomial -->

<!-- api: src/poly/basis.mojo to_polynomial -->

<!-- api: src/poly/basis.mojo coefficients -->

<!-- api: src/poly/basis.mojo nodes -->

<!-- api: src/poly/basis.mojo evaluate -->

<!-- api: src/poly/basis.mojo subdivide -->

<!-- api: src/poly/basis.mojo add -->

<!-- api: src/poly/basis.mojo subtract -->

<!-- api: src/poly/basis.mojo multiply -->

<!-- api: src/poly/basis.mojo negative -->

<!-- api: src/poly/basis.mojo derivative -->

<!-- api: src/poly/basis.mojo integral -->

<!-- api: src/poly/basis.mojo to_core -->

<!-- api: src/poly/basis.mojo from_core -->

Numerical construction uses `context=`; arithmetic keeps precision explicit on every call.
See the [numerical tutorial](../tutorial/array-numeric.md#numerical-series-and-alternative-bases).

<!-- api: src/poly/basis.mojo NumericArray -->
