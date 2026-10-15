# Find and refine numerical root covers

For coefficient-array inputs, start with the [array workflow](array-roots.md).

A root cover tells you where roots can occur and how many lie in each region.
It also records whether the search accounted for every root and whether the
regions meet your accuracy and isolation requirements. Read these claims
separately: a complete cover can contain a cluster that still needs refinement.

## Refine a real cover

The following program finds the real roots of $x^2-2$, checks the cover, and
refines it from a coarse width target to a tighter one.

<!-- example: docs/examples/tutorial_refinement.mojo -->

`refine_real_roots_exact` takes the original polynomial, the previous result,
and both sets of options. It checks and reuses the earlier evidence. Keep the
old options alongside a saved result; they form part of the claim being
checked. Matching entry points support numerical input and complex covers.

## Read the result before using the regions

| Field | Meaning |
|---|---|
| `complete` | The cover accounts for all roots under its domain's counting convention |
| `accurate` | The counted regions satisfy the requested width bound |
| `isolated` | Each counted region contains a single root under that convention |
| `outcome` | Success or the reason a requested property remains unresolved |

Real covers include cells with count zero: these certify excluded regions.
Do not treat every cell as a root interval. A count of `-1` marks an unresolved
cell. Real root counts are **distinct** counts, so a repeated real root counts
once. Complex counts include **multiplicity**: a double root contributes two.

The root checkers verify complete coverage, counts, and the reported accuracy
and isolation claims. They reject an incomplete cover. They can accept a
complete cover whose outcome is inconclusive because accuracy or isolation
remains unresolved. Require the fields your application actually needs after
checking the certificate.

## Complex roots and clusters

Complex roots are enclosed in rectangles with certified zero-free boundaries.
Rational winding evidence establishes the count in each rectangle. Distinct
rectangles have disjoint interiors; a shared edge is allowed when the boundary
checks establish that it contains no root.

Set `require_isolation=False` when a counted cluster is acceptable. A repeated
complex root cannot be split into count-one regions, so demanding isolation
can prevent success even when its location is tightly enclosed. The
[combined numerical example](../spec/numerical.md#example) demonstrates both
real and complex covers.

## Handle uncertain input and search limits

For a numerical point polynomial, stored binary coefficients are treated as
exact point data for certification. For ball coefficients, a result must hold
uniformly for every polynomial in the family. A leading coefficient that may
vanish can leave the degree uncertain; a family can also contain polynomials
with different root counts.

`NumericOutcome` distinguishes uncertain degree, input uncertainty, precision
ceilings, subdivision limits, and unresolved boundaries. A larger budget or
higher precision will not repair a family whose members genuinely have
different roots. See [errors and work limits](../guides/errors.md) for recovery
choices.

With `complete_only=True`, the producer raises unless all requested properties
succeed, including accuracy and isolation. Budget exhaustion always raises;
it does not return a partial result. The [numerical reference](../spec/numerical.md)
defines each outcome and the width formula.
