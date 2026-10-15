# Check an evaluation and request an accuracy

Use `dense.enclose` or `sparse.enclose` when you need an evaluation result that
another calculation can check independently. These functions retain the original
polynomial, argument enclosures, and a trace of rational bounds. `evaluation.check`
verifies that trace without rerunning floating-point evaluation.

<!-- example: docs/examples/array_extension_certified.mojo -->

The polynomial $1+2x+3x^2$ equals 2 at $x=1/3$. The first call encloses that value
at a fixed precision. The second starts at 16 bits and increases the working
precision until it establishes the requested width or reaches 128 bits.
`absolute_bits=64` requests a width no greater than $2^{-64}$. A nonzero
`relative_bits` adds a scale-dependent relative allowance.

## Read the result's claim

`value()` returns an independent `Ball` or `ComplexBall`. `precision()` reports
the last working midpoint precision; `evidence()` borrows its immutable trace.
For adaptive calls, `options()` retains the original request and `outcome()` is
`SUCCESS` or `PRECISION_CEILING`. Fixed-precision calls have no accuracy outcome.

A successful `evaluation.check(result)` proves the enclosure and validates any
reported accuracy outcome. Check `result.accurate()` as well when your application
needs the requested width. The final example evaluates $x$ on $[-1,1]$. Its
uncertainty persists at every precision, so it returns a valid enclosure with
`PRECISION_CEILING`. The checker accepts that truthful outcome.

## Choose inputs and precision

Both workflows accept integer, rational, Gaussian-rational, `AA`, and `QQbar`
coefficient sequences, or prepared numerical arrays. Modular coefficients are
rejected.
Fixed evaluation requires a `BallContext` with an explicit precision. Adaptive
evaluation requires `NumericOptions`; it preserves the original coefficients
through all refinement attempts. Scalar arguments are converted outward to a
ball. Adaptive calls perform that conversion at the configured maximum precision
and retain the resulting argument enclosure throughout refinement.

Sparse inputs use a coordinate tuple in the original axis order; univariate
inputs also accept a scalar argument. Complex coefficients or coordinates select a `ComplexBall` result. Ball coefficients
describe an entire polynomial family, and the returned enclosure covers every
member of that family at every allowed argument. Raising working precision
cannot remove uncertainty supplied with the inputs.

For algebraic coefficients, each attempt refines their exact root bounds at
the current working precision. The result retains those witnesses and the
numerical coefficient enclosures. Its checker verifies the original and refined
root evidence, containment, and the arithmetic trace without repeating root
refinement or floating-point evaluation. Algebraic coordinates also retain
conversion witnesses. Keep the complete result: `evidence()` exposes the core
arithmetic trace, not the additional algebraic conversion evidence.

Every call accepts `max_steps=` or a shared `Budget`. Exhaustion raises
`PolyError` with kind `LIMIT`; it is distinct from completing an adaptive request
at its precision ceiling. See the [evaluation reference](../reference/array-evaluation.md)
for signatures and the [numerical design](../architecture/numerical.md) for the
underlying proof schedule.

## Check one witness step at a time

Use `enclose_stream` and `enclose_checker` when retaining a complete trace would
cost too much memory. They accept the same exact or numerical inputs and
explicit `BallContext` as fixed `enclose`.

<!-- example: docs/examples/array_parity_stream.mojo -->

Forward every event in order. `accept` checks the next event; `finish` checks
completeness and the final receipt. A receipt alone is not a proof. Calling the
producer's `finish` before consuming its steps raises `DOMAIN`. A rejected event
permanently invalidates that checker. If a call raises `LIMIT`, give it more
budget and retry the same event; its position has not advanced.

### Refine algebraic inputs automatically

Use the same calls with `AA` or `QQbar` coefficients or algebraic evaluation
points. Each real and imaginary component is refined at the requested
`BallContext` precision before it enters the arithmetic trace. You don't need
to call `algebraic.refine` first.

<!-- example: docs/examples/array_algebraic_stream.mojo -->

Here the coefficient and evaluation point both represent $\sqrt{2}$. Their
product is $2$, which lies inside the returned enclosure. The checker verifies
the root identities and narrower intervals before checking the multiplication. It does not repeat the
producer's root search. The original coefficients and point are unchanged.

For algebraic inputs, `next` returns an `evaluation.EnclosureStep`. The stream
sends coefficient refinements in canonical coefficient order, algebraic point
refinements in axis order, then arithmetic bounds. Your forwarding loop stays
the same. `step_count()` and the receipt's `steps` include all these events.
Streams without algebraic inputs continue to emit the core `NumericBox` type.
Code that explicitly annotated an algebraic stream as a core `EnclosureStream`
should let Mojo infer the facade stream type instead.

Higher precision narrows the algebraic input intervals and reduces arithmetic
rounding. It does not promise a particular final width: cancellation, repeated
use of the same uncertain value, or numerical input uncertainty can leave a
wider result. Use retained `adaptive_evaluate` when you need an explicit
absolute or relative accuracy target. Streaming evaluation performs one pass
at the supplied precision.

The producer refines one input when you request its event. Neither producer nor
checker keeps a history of refinement witnesses or arithmetic events. They do
retain the input polynomial, prepared instruction graph, and live intermediate
values, so memory depends on the input and calculation. Preparation, refinement,
checking, and completion all use your shared work budget.
