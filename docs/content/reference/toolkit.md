# Polynomial toolkit

The toolkit builds on the existing exact and numerical polynomial types.
Every computational call takes a shared `Budget`; numerical arithmetic also
takes an explicit apn_mojo context. The calculus, interpolation,
rational-function, basis, and series groups each start with a focused example.

<!-- example: docs/examples/toolkit.mojo -->

## Choosing the operation

| Need | Entry points and guarantees |
|---|---|
| Calculus and transforms | `integrate_dense`, `definite_integral`, `compose_dense`, `affine_transform`, `taylor_shift`, `jet_evaluate`, `substitute_many`; corresponding numerical functions use explicit contexts |
| Repeated evaluation or interpolation | `MultipointPlan.prepare` stores a balanced product tree; `.evaluate` permits repeated nodes, `.interpolate` requires distinct nodes |
| Derivative constraints | `hermite_interpolate` takes raw derivatives by default; `derivatives=False` takes Taylor/Hasse coefficients and avoids factorial assumptions |
| Named expressions | `PolynomialRing.create` maps names to stable IDs; `.parse` accepts a bounded polynomial grammar and `.format` round trips exact coefficients |
| Persistence | `serialize_dense`/`serialize_sparse` and numerical counterparts retain domains and representation metadata; decoding is strict |
| Rational functions | Evaluation, differentiation, simultaneous polynomial substitution, univariate rational composition, partial fractions, pole factors, Laurent coefficients and residues |
| Alternative bases | `BasisPolynomial` and `NumericBasisPolynomial` provide Newton, first-kind Chebyshev and Bernstein representations; Bernstein retains its degree |
| Formal series | `TruncatedSeries` and `NumericSeries` carry a known order modulo x^n; binary operations require equal orders, and `.truncate` only reduces precision |
| Root refinement | `refine_real_roots_*` and `refine_complex_roots_*` check earlier evidence, retain satisfactory regions and apply a new subdivision allowance |
| Repeated exact real roots | `real_root_multiplicities` adds a squarefree decomposition and the vanishing block of each root to a complete Descartes result |
| Root diagnostics | Residual, absolute coefficient backward-error and local derivative-sensitivity bounds are rational; sensitivity is a root condition bound only at a root |
| Newton polishing | `polish_root_*` treats iterates as candidates and proves a one-root rectangle independently before reporting success |
| Numerical division | Point domains return rounded q/r and full residual bounds; ball domains also prove uniform enclosure of exact Euclidean q/r |

Basis and series types use exact coefficient fields; `rational_dense` promotes
integer inputs for those operations. Exact antiderivatives require the
coefficient divisions to exist; use
`integrate_integer_dense` or `integrate_integer_sparse` to promote ZZ to QQ.
Chebyshev and Bernstein conversions can reject finite-characteristic fields
when their conversion factors are not invertible. Series composition requires
an exactly zero inner constant, logarithm requires constant one, exponential
requires constant zero, and reversion requires a nonzero linear coefficient.
A formal truncation supplies no analytic tail bound.

Numerical point operations round; ball operations enclose every represented
coefficient selection. Numerical interpolation requires separated node
intervals. In family division, the exact quotient and remainder belonging to
an input pair lie in the reported coefficient boxes; arbitrary independent
selections from those output boxes need not satisfy the division identity.
The separate residual bounds cover `f-q*g-r` even for independent selections.
Newton polishing uses the midpoint representative only to generate candidates
for uncertain coefficient families; its contour witness checks the original
family. A Newton-correction bound is a step-size bound, not a root-distance
bound.

Core partial-fraction checkers prove properness and recomposition without asserting
irreducibility. Pole factors represent all their complex zeros. Laurent
expansions and residues currently take rational centers. General algebraic
coefficient fields are outside these core rational utilities. The
[array rational API](array-rational.md) supports QQ, QQi, AA, and QQbar and adds
complete-normalization and irreducibility checks.

## Calculus and transformations

For $p(x)=3x^2+2x$, `integrate_dense` returns $x^3+x^2$, and the integral
from $0$ to $1$ is exactly $2$. The example also substitutes $x+1$ with
`compose_dense`:

<!-- example: docs/examples/tutorial_calculus.mojo -->

<!-- api: src/polynomials/calculus.mojo derivative_dense -->

<!-- api: src/polynomials/calculus.mojo rational_dense -->

<!-- api: src/polynomials/calculus.mojo rational_sparse -->

<!-- api: src/polynomials/calculus.mojo integrate_dense -->

<!-- api: src/polynomials/calculus.mojo integrate_integer_dense -->

<!-- api: src/polynomials/calculus.mojo definite_integral -->

<!-- api: src/polynomials/calculus.mojo compose_dense -->

<!-- api: src/polynomials/calculus.mojo affine_transform -->

<!-- api: src/polynomials/calculus.mojo taylor_shift -->

<!-- api: src/polynomials/calculus.mojo jet_evaluate -->

<!-- api: src/polynomials/calculus.mojo integrate_sparse -->

<!-- api: src/polynomials/calculus.mojo integrate_integer_sparse -->

<!-- api: src/polynomials/calculus.mojo substitute_many -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_integrate_dense -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_definite_integral -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_compose_dense -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_affine_transform -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_jet_evaluate -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_taylor_shift -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_integrate_sparse -->

<!-- api: src/polynomials/numeric_calculus.mojo numeric_substitute_many -->

## Interpolation

Use `MultipointPlan.prepare` for evaluation, including repeated nodes. Use
`MultipointPlan.prepare_interpolation` when interpolating multiple value lists
at the same distinct nodes: it retains the inverse derivative weights as well
as the product tree. Both plans are immutable; preparation consumes its own
budget, and each interpolation call has a predictable work charge. The ordinary
plan and the one-shot `interpolate` function still work without a weight cache.

The values $1,4,9$ at $0,1,2$ determine $(x+1)^2$. The example builds a plan,
interpolates, checks the result with `check_interpolation`, and reuses the
plan for evaluation:

<!-- example: docs/examples/tutorial_interpolation.mojo -->
For arbitrary-precision exact fields, calls with at least 32 points release coefficient buffers
after their last tree consumer. The reusable product tree stays in the plan.

<!-- api: src/polynomials/interpolation.mojo MultipointPlan -->

<!-- api: src/polynomials/interpolation.mojo interpolate -->

<!-- api: src/polynomials/interpolation.mojo hermite_interpolate -->

<!-- api: src/polynomials/interpolation.mojo check_interpolation -->

<!-- api: src/polynomials/interpolation.mojo check_hermite_interpolation -->

<!-- api: src/polynomials/numeric_interpolation.mojo numeric_multipoint_evaluate -->

<!-- api: src/polynomials/numeric_interpolation.mojo numeric_interpolate -->

<!-- api: src/polynomials/numeric_interpolation.mojo numeric_hermite_interpolate -->

## Construction and interchange

Version-one serializers build one output buffer; decoders validate terms as
they read them into the final polynomial, using at most 64 text rows of scratch
storage for exact sparse input. They preserve canonical bytes and
reject malformed metadata, noncanonical scalars and trailing data. These are
internal memory improvements; the text format and public APIs are unchanged.

<!-- api: src/polynomials/builders.mojo PolynomialRing -->

<!-- api: src/polynomials/interchange.mojo serialize_dense -->

<!-- api: src/polynomials/interchange.mojo deserialize_dense -->

<!-- api: src/polynomials/interchange.mojo serialize_sparse -->

<!-- api: src/polynomials/interchange.mojo deserialize_sparse -->

<!-- api: src/polynomials/numeric_interchange.mojo serialize_numeric_dense -->

<!-- api: src/polynomials/numeric_interchange.mojo deserialize_numeric_dense -->

<!-- api: src/polynomials/numeric_interchange.mojo serialize_numeric_sparse -->

<!-- api: src/polynomials/numeric_interchange.mojo deserialize_numeric_sparse -->

## Rational functions

For $f(x)=x^2/(x-1)$, `partial_fractions` gives $x+1+1/(x-1)$. Its value at
$2$ is $4$, and its residue at the pole $1$ is $1$:

<!-- example: docs/examples/tutorial_rational.mojo -->

<!-- api: src/polynomials/rational_tools.mojo rational_evaluate -->

<!-- api: src/polynomials/rational_tools.mojo rational_derivative -->

<!-- api: src/polynomials/rational_tools.mojo rational_substitute -->

<!-- api: src/polynomials/rational_tools.mojo rational_compose -->

<!-- api: src/polynomials/rational_tools.mojo PartialFractionTerm -->

<!-- api: src/polynomials/rational_tools.mojo PartialFractions -->

<!-- api: src/polynomials/rational_tools.mojo partial_fractions -->

<!-- api: src/polynomials/rational_tools.mojo check_partial_fractions -->

<!-- api: src/polynomials/rational_tools.mojo rational_poles -->

<!-- api: src/polynomials/rational_tools.mojo LaurentExpansion -->

<!-- api: src/polynomials/rational_tools.mojo laurent_expansion -->

<!-- api: src/polynomials/rational_tools.mojo check_laurent_expansion -->

<!-- api: src/polynomials/rational_tools.mojo rational_residue -->

## Alternative bases

The polynomial $x^2$ has Bernstein coefficients $0,0,1$ at degree two and
equals $\tfrac12 T_0(x)+\tfrac12 T_2(x)$ in the Chebyshev basis. The example
also evaluates the Bernstein form at $1/2$:

<!-- example: docs/examples/tutorial_bases.mojo -->

<!-- api: src/polynomials/bases.mojo BasisPolynomial -->

<!-- api: src/polynomials/numeric_bases.mojo NumericBasisPolynomial -->

## Formal series

Modulo $x^4$, the reciprocal of $1-x$ is $1+x+x^2+x^3$:

<!-- example: docs/examples/tutorial_series.mojo -->

<!-- api: src/polynomials/series.mojo TruncatedSeries -->

<!-- api: src/polynomials/numeric_series.mojo NumericSeries -->

## Numerical workflows

<!-- api: src/polynomials/numeric_refine.mojo refine_real_roots_exact -->

<!-- api: src/polynomials/numeric_refine.mojo refine_real_roots_numeric -->

<!-- api: src/polynomials/numeric_refine.mojo refine_complex_roots_exact -->

<!-- api: src/polynomials/numeric_refine.mojo refine_complex_roots_numeric -->

<!-- api: src/polynomials/numeric_multiplicity.mojo RootMultiplicityCertificate -->

<!-- api: src/polynomials/numeric_multiplicity.mojo real_root_multiplicities -->

<!-- api: src/polynomials/numeric_multiplicity.mojo check_real_root_multiplicities -->

<!-- api: src/polynomials/numeric_diagnostics.mojo RootDiagnostics -->

<!-- api: src/polynomials/numeric_diagnostics.mojo PolishedRootResult -->

<!-- api: src/polynomials/numeric_diagnostics.mojo root_diagnostics_exact -->

<!-- api: src/polynomials/numeric_diagnostics.mojo check_root_diagnostics_exact -->

<!-- api: src/polynomials/numeric_diagnostics.mojo polish_root_exact -->

<!-- api: src/polynomials/numeric_diagnostics.mojo check_polished_root_exact -->

<!-- api: src/polynomials/numeric_diagnostics.mojo root_diagnostics_numeric -->

<!-- api: src/polynomials/numeric_diagnostics.mojo check_root_diagnostics_numeric -->

<!-- api: src/polynomials/numeric_diagnostics.mojo polish_root_numeric -->

<!-- api: src/polynomials/numeric_diagnostics.mojo check_polished_root_numeric -->

<!-- api: src/polynomials/numeric_division.mojo NumericDivisionStep -->

<!-- api: src/polynomials/numeric_division.mojo NumericDivisionResult -->

<!-- api: src/polynomials/numeric_division.mojo numeric_divmod -->

<!-- api: src/polynomials/numeric_division.mojo check_numeric_division -->

## Incremental text and CAD transport

The `serialize_*_to` functions write the same version-one bytes as the
String-returning serializers through `PolynomialTextSink.write`. The
`deserialize_*_from` functions consume an owned `PolynomialTextSource` through
end of input. Sources may return arbitrary chunks, including empty chunks;
fields may cross chunk boundaries. The default individual-field limit is
16 MiB and can be raised explicitly. Untrusted record counts do not cause
up-front allocation. The returned polynomial is complete and canonical.

The output chunk target defaults to 4 KiB. A single header or term may exceed
that target; scalar codecs still construct one scalar string. A source controls
its own chunk allocation. A failed operation can leave a partial destination,
which callers must discard rather than treat as a complete record.

`cad_to` calls `CADSink.begin` with shared projection evidence, then writes
`CADStreamEvent` values in depth-first order. Internal events contain a complete
stack, whose `parent` is `level - 1`; terminal events contain evaluations and
signs. `CADStreamChecker.prepare(inputs, budget)` can be passed directly as the
sink. Call its `finish(receipt, budget)` after production. A receipt alone does
not certify a decomposition. Existing `cad` retains its breadth-first result
and original certificate layout.

The checker retains the projection basis and active sample stacks. It checks
and discards completed projection, lifting and sign proof payloads. Empty,
truncated, duplicated and reordered traversals do not certify completion.
After a rejected event or failed checking computation, construct a new checker.

<!-- api: src/polynomials/cad_stream.mojo CADSink -->

<!-- api: src/polynomials/cad_stream.mojo CADStreamChecker -->

<!-- api: src/polynomials/cad_stream.mojo CADStreamEvent -->

<!-- api: src/polynomials/cad_stream.mojo CADStreamReceipt -->

<!-- api: src/polynomials/text_stream.mojo PolynomialTextSink -->

<!-- api: src/polynomials/text_stream.mojo PolynomialTextSource -->

<!-- api: src/polynomials/cad_stream.mojo cad_to -->

<!-- api: src/polynomials/interchange.mojo deserialize_dense_from -->

<!-- api: src/polynomials/interchange.mojo deserialize_sparse_from -->

<!-- api: src/polynomials/interchange.mojo serialize_dense_to -->

<!-- api: src/polynomials/interchange.mojo serialize_sparse_to -->

<!-- api: src/polynomials/numeric_interchange.mojo deserialize_numeric_dense_from -->

<!-- api: src/polynomials/numeric_interchange.mojo deserialize_numeric_sparse_from -->

<!-- api: src/polynomials/numeric_interchange.mojo serialize_numeric_dense_to -->

<!-- api: src/polynomials/numeric_interchange.mojo serialize_numeric_sparse_to -->

### Complete streaming example

<!-- example: docs/examples/streaming.mojo -->
