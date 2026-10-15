# Handle errors and work limits

Polynomial algorithms can encounter invalid input, grow beyond representable
sizes, or need more work than an application can afford. Computational entry
points raise `PolyError` and share a caller-owned `Budget`, so you can keep
these cases separate from mathematical results.

<!-- example: docs/examples/tutorial_errors.mojo -->

The first error is deliberate work exhaustion. The second is a name that does
not belong to the ring. Neither produces a polynomial that could be mistaken
for a successful result.

## Set a budget for a workflow

Pass the same `Budget` through related calls to limit their combined work.
`Budget()` uses the maximum `UInt64` allowance; choose a finite allowance when
your application needs one. A rejected charge leaves the remaining allowance
unchanged. Algorithms charge before the corresponding work and preserve
borrowed inputs on failure.

Budget units are implementation work estimates, not milliseconds, bytes, or a
stable operation count across versions. They do not provide a hard memory
limit. Cache state can affect charges without changing the mathematical result.
Measure representative inputs when choosing an application limit.

## Respond to the error kind

| Kind | Typical cause | Useful response |
|---|---|---|
| `DOMAIN` | Zero divisor, incompatible rings, invalid nodes, or unknown variable | Correct the input or choose an appropriate domain |
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

The [numerical root tutorial](../tutorial/numerical-roots.md) shows how to
refine checked evidence. The [API reference](../reference/p1.md) defines
`Budget`, `PolyError`, and `PolyErrorKind`.
