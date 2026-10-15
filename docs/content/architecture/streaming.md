# Streaming and memory lifetimes

Large polynomial computations can retain much more than the final answer:
partial products, scalar encodings, evaluation traces, and root or CAD proofs.
The implementation reduces this storage in two ways: it bounds intermediate
polynomial work where the algorithm permits, and it lets consumers check or
write completed payloads as they arrive.

## Avoid unnecessary intermediate polynomials

Ordered sparse merge and product helpers feed terms to their consumers without
first collecting every contribution. Bounded products retain only coefficients
or terms relevant to a truncation. Occupied exponent grouping processes the
support that exists rather than allocating a full degree-sized container.

These techniques still require worklists, accumulators, or retained output.
They also preserve order: numerical coefficient accumulation must follow its
specified rounding schedule, even if a different grouping would save storage.

## Separate materialized and incremental interfaces

| Payload | Materialized route | Incremental route |
|---|---|---|
| Polynomial text | One complete string | Chunks written to a sink or read from a source |
| Enclosure evidence | A retained evaluation trace | Witness steps consumed with live-slot tracking |
| Root evidence | A result holding the cover and proofs | Ordered events followed by a receipt |
| CAD evidence | A complete decomposition and certificate | Depth-first events with active stack state |

Both routes preserve their documented mathematical claims. The incremental
route does not require callers to collect events into a second materialized
copy before checking them.

## Text transport

Streaming serializers emit the same version-one bytes as their string-returning
counterparts. The default chunk target is 4 KiB, but a header, term, or scalar
encoding may exceed that target. Scalar codecs still form one scalar string.

Readers accept arbitrary chunk boundaries and impose an individual-field
limit, 16 MiB by default. Untrusted counts do not trigger an immediate allocation
of the claimed size. A failed write can leave a partial destination; callers
must discard it. These are bounded transport operations, not a promise that
all allocations fit within the chunk target.

## Incremental proof checking

An evaluation checker retains the graph, slot metadata, and currently live
values needed to verify later steps. Completed trace payloads can be released.
Root checking retains the ordering, count, and boundary information needed to
establish global coverage. CAD checking retains shared projection data and
active sample stacks while releasing completed lifting and sign evidence.

A stream receipt summarizes production; it is not a standalone certificate.
The checker must bind the received events to the original inputs and validate
completion with `finish`. Missing, duplicate, or reordered proof events cannot
be accepted as a complete result. After a rejected event or a failed checking
computation, start with a new checker.

The [interchange guide](../guides/interchange.md) demonstrates the full lifecycle.
The numerical and toolkit references describe the event and receipt fields.

## Evaluate the tradeoff

Streaming saves retained payloads, but graphs, live values, search state,
inputs, and requested outputs still occupy memory. Some checks also perform
work earlier or repeat small computations to avoid keeping large data.

Measure producer, consumer, and retained output separately when investigating
memory. Streaming bounds what a consumer retains; it is not a universal
speedup for every call.
