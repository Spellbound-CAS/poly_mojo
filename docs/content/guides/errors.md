# Handle errors and work limits

Polynomial algorithms can encounter invalid input, grow beyond representable
sizes, or need more work than an application can afford. Computational entry
points raise `PolyError` and share a caller-owned `Budget`, so you can keep
these cases separate from mathematical results.

<!-- example: docs/examples/tutorial_errors.mojo -->

The first error is deliberate work exhaustion. The second is division by the
zero polynomial. Neither produces a polynomial that could be mistaken for a
successful result.

## Set a budget for a workflow

Pass the same `Budget` through related calls to limit their combined work.
`Budget()` uses the maximum `UInt64` allowance; choose a finite allowance when
your application needs one. A rejected charge leaves the remaining allowance
unchanged. Algorithms charge before the corresponding work and preserve
borrowed inputs on failure.

### What a step costs

A budget unit is one step of the library's work model. Each call charges its
steps before doing the corresponding work:

| Work | Steps charged |
|---|---|
| Visiting a term, coefficient, or pending subproblem | 1 |
| Allocating storage | 1 per element |
| Adding or subtracting two integers | the sum of their sizes |
| Multiplying two integers; adding, subtracting, or multiplying two rationals | the square of the sum of their sizes |
| An operation on prime-field residues | 1 |

An integer's size is its magnitude in bytes, plus one; a rational's size is
the sum of its numerator's and denominator's sizes. Numerical coefficients are
charged the same way, by the storage their precision needs. Computing
$(1+x)^{100}$ with `poly.dense.power`, for example, takes about 330,000 steps;
its coefficients have up to 30 digits.

The same call on the same input charges the same steps on every machine,
whatever its load or speed. Reusing prepared state from an earlier call can
lower the charge without changing the result. To measure a call, pass it a
fresh `Budget()` and subtract its `remaining()` from `UInt64.MAX` afterwards,
as [the first tutorial](../tutorial/first-polynomial.md#bound-a-calculation)
shows.

Steps do not convert to time or memory. On one laptop core, different
operations ran at anywhere from about 100 million to over 10 billion steps
per second: the table charges a multiplication of large numbers by the square
of their sizes, while the processor multiplies eight bytes at a time. An
allocation is charged per element, not per byte, so steps are not a hard
memory limit either. Step counts can change between releases, so measure
representative inputs when choosing an application limit, and check it again
after upgrading.

## Respond to the error kind

| Kind | Typical cause | Useful response |
|---|---|---|
| `DOMAIN` | Zero divisor, incompatible domains or moduli, invalid nodes, or an invalid axis | Correct the input or choose an appropriate domain |
| `LIMIT` | Exhausted work allowance or an unmet strict search limit | Increase the relevant limit if more work is acceptable |
| `CAPACITY` | A size, exponent, or encoding exceeds its representable bound | Reduce or restructure the problem |
| `UNSUPPORTED` | The requested operation cannot handle this case | Choose a supported formulation or algorithm |
| `ENCODING` | Malformed or noncanonical interchange data | Reject or correct the data |

Do not interpret an error as zero, nonmembership, absence of roots, or a false
constraint. Those are mathematical answers with their own successful return
values. Catch `PolyError` where the application can make a recovery decision;
let it propagate when it cannot.

## Numerical outcomes are separate

A numerical producer can return useful certified information without meeting
all requested targets. Its `NumericOutcome` explains the remaining obstacle:

| Outcome | What to inspect |
|---|---|
| `SUCCESS` | The result and its checker, plus the claims your application needs |
| `UNCERTAIN_DEGREE` | Whether the leading coefficient can be zero |
| `INPUT_UNCERTAINTY` | Whether different members of the input family require different answers |
| `PRECISION_CEILING` | `initial_precision`, `max_precision`, and the requested width |
| `SUBDIVISION_LIMIT` | `max_subdivisions` and unresolved or overly wide regions |
| `UNRESOLVED_BOUNDARY` | Whether the search has certified a zero-free boundary |

For root covers, inspect `complete`, `accurate`, and `isolated` independently.
A complete cover may still be too wide, while an incomplete cover cannot
establish a global root count. `complete_only=True` turns failure to meet the
requested root properties into an exception. Budget exhaustion always raises,
including in modes that otherwise return inconclusive numerical results.

The [root tutorial](../tutorial/array-roots.md#choose-the-numerical-claim-you-need)
shows how to read and refine checked evidence. The
[API reference](../reference/array-api.md#shared-types) defines `Budget`,
`PolyError`, and `PolyErrorKind`.
