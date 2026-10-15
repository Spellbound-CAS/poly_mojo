# Arbitrary-precision numerical polynomials

These APIs provide rounded arithmetic, rigorous evaluation, and validated
real and complex root covers. Start with [numerical evaluation](../tutorial/numerical.md)
or [root refinement](../tutorial/numerical-roots.md) for a worked example.
Numerical declarations have a separate compatibility snapshot and do not
extend the exact kernel codec.

## Domains and operations

| Adapter | Scalar | Explicit context | Result semantics |
|---|---|---|---|
| `FloatDomain` | `apn_mojo.Float` | `ArithmeticContext` | Rounded real point arithmetic |
| `ComplexDomain` | `apn_mojo.Complex` | `ComplexContext` | Rounded complex point arithmetic |
| `BallDomain` | `apn_mojo.Ball` | `BallContext` with precision | Encloses all real inputs |
| `ComplexBallDomain` | `apn_mojo.ComplexBall` | `BallContext` with precision | Encloses all complex inputs |

`NumericDensePoly[D]` stores ascending univariate coefficients.
`NumericSparsePoly[D]` uses the existing monomials and orders. Both implement
addition, subtraction, negation, multiplication, nonnegative powers, scalar
multiplication/division, differentiation and evaluation. Dense evaluation also
provides a simultaneous value/derivative operation. Sparse evaluation caches
variable powers within the call. Scalar division requires a divisor that
certainly excludes zero. `numeric_divmod` provides polynomial division with residual or family-enclosure
certificates; see the [toolkit reference](../reference/toolkit.md#numerical-workflows).

Only finite coefficients and bounded enclosures are accepted. Constructors
remove point zeros or singleton-zero balls, never coefficients that merely
contain zero. `degree_bound` and `certified_degree` distinguish possible degree
from a degree shared by the entire coefficient family. Public fields support
certificate transport and inspection; callers must preserve container invariants.

Precision is never ambient. Scalar operations follow the supplied apn_mojo
contexts, including rounding and traps. Operands retain their stored values.
Dense convolution visits increasing left then right indices; Horner descends
through coefficients; powers use binary powering. Sparse duplicate coefficients
are accumulated stably in input order after monomial sorting and removal of
certain-zero input terms. An intermediate zero is retained until its duplicate
group ends, so subsequent additions still receive the prescribed rounding. A coefficient
that needs no sparse arithmetic retains its original format. Integer derivative
multipliers are converted under the context before multiplication. These are
sequences of rounded operations, not a promise of one correctly rounded result
for the entire polynomial calculation. No algebraic reassociation is applied.

Explicit conversion functions support all four domains. Exact ZZ/QQ/QQi conversion
constructs balls directly from rational bounds. Modular polynomials, including
zero polynomials, have no implicit numerical embedding. Conversion to a point
uses apn_mojo’s public certainty predicate, preserving its signed-zero and
rounding-status distinctions; otherwise it raises `DOMAIN`. Exact singleton
rationals are rounded directly. This never silently substitutes
a midpoint. Numerical domains do not implement the exact coefficient traits.

## Enclosure certificates and accuracy

`enclose_dense` and `enclose_sparse` evaluate existing ZZ/QQ/QQi polynomials.
`enclose_numeric` and `enclose_numeric_sparse` evaluate numerical families.
The output domain is Ball or ComplexBall. Each result carries a deterministic
trace of exact rational rectangles. Its `value(budget)` method returns an owned scalar.
The checker reconstructs the evaluation graph from the original polynomial and
assignments, checks every local rational bound, and checks the reported value.
It performs no floating-point arithmetic or adaptive search.

Adaptive evaluation starts at `initial_precision`, doubles it, and caps the
last attempt at `max_precision`. It rebuilds coefficient enclosures from the
original exact coefficients and preserves argument uncertainty. `SUCCESS` means
both coordinate widths are at most

`2^-absolute_bits + 2^-relative_bits * lower_norm(region)`.

A zero exponent disables that term. `lower_norm` is a lower bound on the
infinity norm throughout the rectangle, so relative accuracy is uniform over
its values. Relative tolerance alone cannot establish a positive width at zero.
A precision ceiling returns the last valid enclosure with an inconclusive
outcome. `check_numeric_accuracy` checks its accuracy claim; the corresponding
evaluation checker must also be called.

## Real and complex roots

Real root entry points accept exact ZZ/QQ and real numerical coefficients.
Finite Float coefficients denote their exact stored binary rationals. Point
inputs use one shared Sturm certificate and a bounded subdivision worklist;
counts are **distinct real roots**. Repeated roots do not defeat this path.
Real ball families use uniform image exclusion, endpoint signs, derivative
bounds and interval Newton predicates. Their output includes excluded cells,
forming an ordered cover of a checked global Cauchy interval. If two explicit
endpoint coefficient choices have different exact global Sturm counts, the
request returns `INPUT_UNCERTAINTY` immediately; equal counts do not establish a
uniform claim. All counted roots
are strictly inside their cells. There is no inferred multiplicity for a ball
family and no automatic multiplicity attachment for exact real inputs.

Complex entry points accept exact ZZ/QQ/QQi and every numerical domain. Each output
rectangle has a zero-free boundary trace and an exact rational winding count,
**including multiplicity**. The implementation translates coefficients around
each boundary segment’s midpoint using rational interval arithmetic before evaluating its
boundary. This reduces cancellation around repeated roots while enclosing the
original family. Boundary segments are oriented counterclockwise and cover all
four sides. The checker computes the winding of an exact representative's image
polygon, validates every enclosure, and checks disjoint interiors and a count
sum equal to the uniform degree. Shared rectangle edges are permitted because
both boundaries are certified root-free.

Search is deterministic and bounded. Complex search bisects boundary segments
and tries asymmetric spatial cuts at fractions 2/7, 3/7 and 4/7 along the longer
axis. Each trial receives at most `32 + 8*degree` boundary subdivisions, within
the shared `max_subdivisions` allowance, leaving room for another contour.
The initial global contour can use the remaining allowance. Precision retries
occur when exact interval bounds show that additional bits may help. Failed
contours never discard a previously counted region.

`complete`, `accurate` and `isolated` are separate claims. A complete cover can
remain too wide or contain a multiple-root cluster. `require_isolation=False`
permits counted clusters. Otherwise an inseparable repeated complex root returns
an inconclusive result, retaining its checked cluster cover when available.
The `refine_real_roots_*` and `refine_complex_roots_*` entry points check a
previous cover against its original options, retain satisfactory regions, and
continue with the new options. The [refinement tutorial](../tutorial/numerical-roots.md)
shows how to keep the result and options together.

A nonzero constant has a successful empty root set. The zero polynomial raises
`DOMAIN`. An uncertain leading coefficient returns `UNCERTAIN_DEGREE`, never a
successful empty set. Exhausted subdivisions or unresolved contours return
explicit outcomes. `complete_only=True` raises unless **all requested claims**
succeed, including accuracy and isolation. It maps inherent unsupported input
to `UNSUPPORTED` and search ceilings to `LIMIT`. Work-budget exhaustion always
raises `LIMIT`; no partial result is returned for work exhaustion.

Root checkers verify complete covers, including valid complete cluster covers
with an inconclusive accuracy/isolation outcome. They return false for incomplete
analyses. Checkers validate counts and coverage independently of search; they do
not certify that a resource ceiling was reached at precisely the reported step.
All numerical checkers share the caller's work budget and preserve borrowed
inputs on failure. Unrepresentable sizes raise `CAPACITY`.

The [theory page](../theory/numerical.md) states the proof obligations.

## Example

<!-- example: docs/examples/numerical.mojo -->

## API

<!-- api: src/polynomials/numeric_domains.mojo NumericDomain -->

<!-- api: src/polynomials/numeric_domains.mojo FloatDomain -->

<!-- api: src/polynomials/numeric_domains.mojo ComplexDomain -->

<!-- api: src/polynomials/numeric_domains.mojo BallDomain -->

<!-- api: src/polynomials/numeric_domains.mojo ComplexBallDomain -->

<!-- api: src/polynomials/numeric_dense.mojo NumericDensePoly -->

<!-- api: src/polynomials/numeric_sparse.mojo NumericTerm -->

<!-- api: src/polynomials/numeric_sparse.mojo NumericSparsePoly -->

<!-- api: src/polynomials/numeric_bounds.mojo NumericBox -->

<!-- api: src/polynomials/numeric_bounds.mojo numeric_bounds -->

<!-- api: src/polynomials/numeric_bounds.mojo numeric_value -->

<!-- api: src/polynomials/numeric_convert.mojo convert_numeric_dense -->

<!-- api: src/polynomials/numeric_convert.mojo convert_numeric_sparse -->

<!-- api: src/polynomials/numeric_convert.mojo numeric_from_dense -->

<!-- api: src/polynomials/numeric_convert.mojo numeric_from_sparse -->

<!-- api: src/polynomials/numeric_options.mojo NumericOutcome -->

<!-- api: src/polynomials/numeric_options.mojo NumericOptions -->

<!-- api: src/polynomials/numeric_evaluate.mojo EnclosureEvaluation -->

<!-- api: src/polynomials/numeric_evaluate.mojo AdaptiveEvaluation -->

### Incremental evaluation witnesses

`EnclosureStream` and `EnclosureStreamChecker` are additive alternatives for
transporting the existing evaluation witness one `NumericBox` at a time. Their
`from_dense`, `from_sparse`, `from_numeric` and `from_numeric_sparse` factories
own the prepared graph, so later edits to the original inputs cannot affect
an active producer or checker. Exact inputs are ZZ/QQ/QQi; the four numerical
coefficient families are accepted by the numerical factories. Output domains
remain Ball and ComplexBall, with an explicit producer context.

Call `next(budget)` until it returns None. Only then may `finish(budget)` return
a `StreamedEnclosure` receipt. Every successful `next` returns the corresponding
box from the existing materialized certificate, in the same instruction order.
An abandoned prefix is incomplete evidence. A failed computation does not
consume a step, allowing retry with a fresh allowance.

The independent checker receives each box through `accept(box, budget)` and
uses exact rational arithmetic to check that instruction. It retains bounds
with future uses and recycles other slots. A malformed or extra box permanently
invalidates the checker. Its `finish(receipt, budget)` returns True only after
the complete stream, the reported step count and precision, and the final
scalar enclosure have all been checked. A receipt alone is not a proof.

Producer and checker can be connected in one loop or separated by a caller's
transport/storage layer. Neither requires storing all witness payloads.
Memory still includes the input graph and live scalar/bound slots; dense input
coefficients are loaded before their Horner consumers, so this is not a
constant-memory algorithm in the polynomial size. Existing APIs keep their
full `EnclosureEvaluation.steps` list and unchanged certificate format.

<!-- api: src/polynomials/numeric_stream.mojo StreamedEnclosure -->

<!-- api: src/polynomials/numeric_stream.mojo EnclosureStream -->

<!-- api: src/polynomials/numeric_stream.mojo EnclosureStreamChecker -->

<!-- api: src/polynomials/numeric_evaluate.mojo enclose_dense -->

<!-- api: src/polynomials/numeric_evaluate.mojo check_enclose_dense -->

<!-- api: src/polynomials/numeric_evaluate.mojo enclose_numeric -->

<!-- api: src/polynomials/numeric_evaluate.mojo check_enclose_numeric -->

<!-- api: src/polynomials/numeric_evaluate.mojo enclose_sparse -->

<!-- api: src/polynomials/numeric_evaluate.mojo check_enclose_sparse -->

<!-- api: src/polynomials/numeric_evaluate.mojo enclose_numeric_sparse -->

<!-- api: src/polynomials/numeric_evaluate.mojo check_enclose_numeric_sparse -->

<!-- api: src/polynomials/numeric_evaluate.mojo adaptive_enclose_dense -->

<!-- api: src/polynomials/numeric_evaluate.mojo adaptive_enclose_sparse -->

<!-- api: src/polynomials/numeric_evaluate.mojo check_numeric_accuracy -->

<!-- api: src/polynomials/numeric_real_roots.mojo RealRootCell -->

<!-- api: src/polynomials/numeric_real_roots.mojo RealRootResult -->

<!-- api: src/polynomials/numeric_real_roots.mojo real_roots_exact -->

<!-- api: src/polynomials/numeric_real_roots.mojo real_roots_numeric -->

<!-- api: src/polynomials/numeric_real_roots.mojo check_real_roots_exact -->

<!-- api: src/polynomials/numeric_real_roots.mojo check_real_roots_numeric -->

<!-- api: src/polynomials/numeric_complex_roots.mojo ComplexBoundarySegment -->

<!-- api: src/polynomials/numeric_complex_roots.mojo ComplexRootRegion -->

<!-- api: src/polynomials/numeric_complex_roots.mojo ComplexRootResult -->

<!-- api: src/polynomials/numeric_complex_roots.mojo complex_roots_exact -->

<!-- api: src/polynomials/numeric_complex_roots.mojo complex_roots_numeric -->

<!-- api: src/polynomials/numeric_complex_roots.mojo check_complex_roots_exact -->

<!-- api: src/polynomials/numeric_complex_roots.mojo check_complex_roots_numeric -->

## General toolkit workflows

The toolkit adds numerical calculus, interpolation, basis conversion, formal series,
lossless interchange, resumable real/complex root refinement, exact real-root
multiplicity certificates, diagnostic bounds, certified Newton polishing and
numerical division. These additive APIs are documented in the
[toolkit reference](../reference/toolkit.md). Existing result and certificate
formats remain unchanged.

## Incremental root witnesses

`real_roots_exact_to`, `real_roots_numeric_to`, `complex_roots_exact_to` and
`complex_roots_numeric_to` send completed local witnesses directly to a sink.
The materialized APIs and stream producers share their refinement workers.
Real streams start with borrowed shared Sturm evidence, absent for uncertain
coefficient families, and emit an ordered cover including excluded cells.
Complex streams emit positive-count rectangles in the original traversal order.

Prepare `RealRootStreamChecker` or `ComplexRootStreamChecker` from the original
polynomial and options, use it as the sink, then call `finish(receipt, budget)`.
For transported evidence, call `begin` for the real header, `accept` for each
local witness, and `finish` for the receipt. Header/witness rejection or an error
during `begin` or `accept` invalidates the checker permanently. `finish` checks
the receipt without changing state. The checker never reruns root isolation.

`RootStreamReceipt` reports outcome, completeness, accuracy, isolation and the
number of transported pieces. It requires the complete accepted witness stream;
its fields alone establish no mathematical claim. The independent checkers
accept only certified local counts and complete covers. Inconclusive producers
can return an incomplete receipt to a general sink; such a receipt does not
pass `finish`. A checker used directly as a sink raises on unresolved local
counts. `complete_only=True` can raise after transporting a prefix, which must
not be treated as a completed proof.

Real checking stores coefficient bounds, shared Sturm evidence and a current
endpoint. Complex checking stores rectangle geometry for pairwise disjointness;
completed contour traces are discarded. Producers retain their active search
frontier. One cell or contour witness is still materialized at a time.

<!-- api: src/polynomials/numeric_complex_roots.mojo ComplexRootSink -->

<!-- api: src/polynomials/root_stream_check.mojo ComplexRootStreamChecker -->

<!-- api: src/polynomials/numeric_real_roots.mojo RealRootSink -->

<!-- api: src/polynomials/root_stream_check.mojo RealRootStreamChecker -->

<!-- api: src/polynomials/root_stream.mojo RootStreamReceipt -->

<!-- api: src/polynomials/numeric_complex_roots.mojo complex_roots_exact_to -->

<!-- api: src/polynomials/numeric_complex_roots.mojo complex_roots_numeric_to -->

<!-- api: src/polynomials/numeric_real_roots.mojo real_roots_exact_to -->

<!-- api: src/polynomials/numeric_real_roots.mojo real_roots_numeric_to -->
