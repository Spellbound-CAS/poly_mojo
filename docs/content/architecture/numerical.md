# Numerical arithmetic and certification

Numerical computation has two routes: point operations return rounded values,
and validated operations establish bounds for a point polynomial or a family
of polynomials. Explicit scalar contexts control arithmetic; certificates
make the bounds independently checkable.

## Domains and operation schedules

`FloatDomain` and `ComplexDomain` represent stored binary point values.
`BallDomain` and `ComplexBallDomain` represent intervals and rectangles. A
coefficient family selects one value from each enclosure and holds those
choices fixed throughout a polynomial calculation.

The numerical adapters do not implement `CoefficientRing` or
`CoefficientField`. Rounded equality and interval containment cannot support
the algebraic laws assumed by exact gcd, division, or Gröbner algorithms.
Numerical modules can reuse an algorithm's traversal only when its arithmetic
schedule and domain conditions remain valid.

For example, dense Horner evaluation descends through coefficients, convolution
visits prescribed coefficient indices, and sparse duplicate accumulation is
stable. Changing the grouping of additions can change rounded results and
ball widths. The [numerical reference](../spec/numerical.md) makes these
schedules part of the operation's semantics.

## Evaluation graphs and rational witnesses

An enclosure producer builds a deterministic computation graph from the
polynomial and assignments. It evaluates that graph with ball arithmetic and
records rational rectangles for the intermediate results.

The checker reconstructs the graph from the original inputs. For each step,
it verifies the rational relation between operand and result bounds, then
checks the reported enclosure. It does not call the floating-point producer,
repeat an adaptive search, or accept a supplied graph without binding it to
the polynomial.

Sparse evaluation reuses variable powers and assignment lookup within a call.
The resulting order of power definitions and witness steps is deterministic.
The incremental route can release completed payloads while retaining the graph
and live values needed by future steps; see [streaming](streaming.md).

## Accuracy is a separate property

A valid interval may be too wide for an application. Adaptive evaluation
doubles working precision up to a ceiling and checks the requested absolute
and relative width bound. Exact coefficients can be reenclosed at the new
precision, but uncertainty supplied in arguments or coefficient balls remains.

The accuracy checker verifies the width claim. It must be paired with the
evaluation checker, which establishes that the interval actually encloses the
answer. A precision-ceiling outcome can therefore accompany useful, valid
bounds without meeting the requested target.

## Root covers

For exact real point data, Sturm information supplies distinct-root counts.
Real ball families use uniform exclusion, sign, derivative, and interval-Newton
evidence. A global bound and ordered cells establish coverage, including
regions certified to contain no root.

Complex covers use rectangles with zero-free boundaries and rational winding
evidence. Their counts include multiplicity. Disjoint interiors and a total
count equal to the certified degree establish completeness. This explains why
a repeated complex root may have an accurate counted region without being
isolated into count-one regions.

Search and proof are separate here too. Precision or subdivision limits can
leave a requested property unresolved. Refinement checks previous evidence and
reuses satisfactory regions; it does not turn an unchecked cached result into
a premise. `complete`, `accurate`, and `isolated` retain distinct meanings.

## Numerical division and diagnostics

Point-domain division returns rounded quotient and remainder coefficients plus
bounds on the residual. Ball-domain division additionally encloses the exact
Euclidean quotient and remainder for every compatible input pair. The choices
of output coefficients can be correlated: arbitrary independent selections
from the quotient and remainder boxes need not satisfy the identity. Separate
residual bounds cover those independent selections.

Root diagnostics and Newton steps help generate and assess candidates. A small
Newton correction is not a root-distance proof. Successful polishing relies
on a separate one-root contour check against the original polynomial family.
These distinctions keep useful numerical estimates from being mistaken for
certified root counts.
