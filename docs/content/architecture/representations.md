# Representations and shared kernels

The public containers keep the polynomial's mathematical data separate from
how an algorithm chooses to process it. Canonical storage supports equality,
reproducible certificates, predictable serialization, and conversion between
dense and sparse computations.

## Exact storage invariants

`DensePoly[C]` stores coefficients in ascending degree order, with a variable
ID and coefficient context. Trailing zeros are removed. The zero polynomial
has an empty coefficient list and degree -1.

`SparsePoly[C]` stores nonzero terms in increasing monomial order. A monomial
contains increasing variable IDs paired with positive exponents. Construction
merges duplicate monomials, removes zero coefficients, and retains the ring
context and chosen order. Leading terms therefore appear at the end of the
stored term list.

The order is part of the representation. `LEX`, `GRLEX`, and `GREVLEX` support
algorithms that require an admissible order. `KERNEL` is a canonical storage
order; algorithms such as Gröbner bases reject it. Exact sparse division can
use an internal admissible order and restore the caller's storage order where
its contract permits that conversion.

Variable names belong to `PolynomialRing`, which maps them to stable numeric
IDs. Algebra operates on IDs. This avoids coupling arithmetic or certificate
formats to a display spelling.

## Numerical storage has stricter semantic limits

Numerical dense and sparse containers use the same broad storage shapes, with
a numerical domain adapter rather than an exact ring trait. Constructors
accept finite point values and bounded balls. They remove only known zeros:
a ball that contains zero is still a possibly nonzero coefficient.

Consequently, `degree_bound` and `certified_degree` answer different questions.
The highest stored coefficient supplies an upper bound; a common degree for
a family needs a leading coefficient that cannot vanish.

Ordering also fixes the arithmetic schedule. Rounded addition is not
associative, so merging a duplicate group or accumulating a product in a
different order can change the output. Shared sparse streams preserve the
specified stable order, including intermediate zeros that must remain until
a duplicate group is complete.

## Shared private operations

Small scalar adapters let vector algorithms share a recurrence without making
numerical scalars pretend to form an exact field. They carry exact arithmetic
and work charging on one route, and explicit numerical contexts on the other.
Basis conversion, interpolation, calculus, and series operations use this
separation where their recurrences agree.

Other private helpers handle ordered term streams, occupied exponent groups,
assignment lookup, reusable powers, and bounded products. Their purpose is to
avoid repeated traversal and large temporary representations. A bounded
product can discard terms outside a known truncation before retaining them;
a sparse grouping routine visits occurring exponents instead of allocating
an entire degree-sized box.

Reuse must preserve the public semantics. In particular, an exact fast
multiplication or reassociation is not automatically valid for a numerical
convolution with a prescribed rounding schedule. The [numerical architecture](numerical.md)
explains that boundary.

## Prepared data and local state

`MultipointPlan` stores a product tree for repeated evaluation.
`prepare_interpolation` additionally caches inverse derivative weights for
repeated interpolation; the evaluation-only `prepare` avoids that cost.
Prepared plans are reusable values, while an operation's mutable worklists,
reduction state, and caches belong to that operation.

This distinction matters for both cost and correctness. Reusing a prepared
polynomial structure can save work. Reusing a proof node for different inputs
would change the claim and must be rejected. Cache lifetime and keys are
therefore part of the internal design even when they are not public API.

Public fields make inspection and certificate transport possible, but do not
remove the invariants. Prefer validated factories to direct field mutation.
See [representations and arithmetic](../reference/p1.md) and the
[compatibility guarantees](../reference/contract-p.md#semantic-guarantees).
