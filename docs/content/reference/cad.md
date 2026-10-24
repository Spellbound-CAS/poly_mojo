# Cylindrical algebraic decomposition

For the array workflow, see [inequalities and CAD](../tutorial/array-cad.md).

## Final-data construction

`cad_final` builds a complete CAD and independently verifies its final data
before returning. Its default combines Brown–McCallum projection reductions
with Lazard lifting,
one algebraic field per base sample, relative sections and bounded
interval-first root isolation. The verifier checks projection, specialization,
root completeness, geometry and original-input signs; it does not replay
intermediate construction proofs. Checked exact SIMD batches accelerate
eligible rational sign evaluations, with scalar and arbitrary-precision
fallbacks.

`CADFinalCertificate` is the evidence returned by `poly.cad.decompose` and
`poly.cad.cover`. It exposes metadata methods and explicit `cells(budget)` and
`sample(budget)` materialization. `check_cad_final` binds inputs, order, fixed
prefix and conjunction independently. Explicit McCallum requests and fixed
prefixes retain the compatible legacy backend and undergo its complete check.

<!-- api: src/polynomials/cad_final_types.mojo CADFinalCell -->

<!-- api: src/polynomials/cad_final.mojo CADFinalCertificate -->

<!-- api: src/polynomials/cad_final.mojo cad_final -->

<!-- api: src/polynomials/cad_final.mojo check_cad_final -->

## Legacy certificates

`cad` constructs a complete sign-invariant decomposition of real coordinate
space. Variables are supplied in **base-to-top** order; every input variable
must occur in that list. Extra variables are allowed. Inputs may share factors,
contain repeated factors, be constant, or be zero. The producer builds primitive
squarefree pairwise coprime projection bases before lifting.

McCallum projection includes each member's coefficients from the leading one
down to the first nonzero constant (every coefficient when none is constant),
discriminants and pairwise resultants. `check_cad_projection` also accepts the
full coefficient set in older certificates. Lazard projection includes leading and trailing coefficients,
discriminants and pairwise resultants. Both use `Res(f,f')` in place of the
discriminant: the extra leading-coefficient factor is already projected.
Projection certificates prove the **Sylvester determinants**, not merely
membership in the input ideal. Their checkers verify fraction-free
multiplication identities without repeating exact elimination division.

Lazard lifting removes vanishing factors before each successive substitution.
The implementation uses the first nonzero Taylor derivative at each prefix;
the result differs from division by linear factors only by a positive factorial
scalar. All earlier zero derivatives and the first nonzero derivative are
certified. This handles algebraic coordinates and nullification exactly.

McCallum mode fails with `UNSUPPORTED` on nullification except at the final
lifting level over a point base. In particular, it does not claim order
invariance for a nullified intermediate projection polynomial. Use
`lazard=True` for such inputs.

## Fix leading algebraic coordinates

`CADFiber(fixed, budget)` binds an exact prefix in the supplied variable order.
Its `decompose`, `cover`, `solve`, and `decompose_to` methods follow only those
sections, then cover or search the remaining coordinates. Include each fixed
coordinate's defining polynomial in the input family: a prefix coordinate must
be a projection section, not a point inside a sector. Invalid witnesses, excess
coordinates, and a nonsection prefix raise `DOMAIN` when encountered. Valid
witnesses for the same value may use different defining polynomials or
isolating intervals; their equality is checked independently.

Use the same `CADFiber` with `check` or `checker` to verify its results. These
methods retain the ordinary certificate and stream layouts, but their claim is
coverage over the bound prefix. Ordinary `check_cad` still requires full-space
coverage. An empty prefix recovers the ordinary traversal. All preparation,
copying, computation, and checking share the supplied budget.

The [coefficient-array interface](array-cad.md) prepares and checks this prefix
automatically for AA coefficients and user-supplied `fixed=` coordinates.
Most callers can use `poly.cad` directly; the
[workflow tutorial](../tutorial/array-cad-workflows.md#fix-coordinates-before-decomposing)
shows how to select a slice without constructing core witnesses.

<!-- api: src/polynomials/cad_fiber.mojo CADFiber -->

## Reading cells and certificates

The decomposition's root cell has index 0, level -1 and no coordinate.
Other cells have a parent, coordinate level and alternating stack position:
even positions are open sectors and odd positions are sections. Each stack
lists all distinct roots at the base sample. Its `sources` pair identifies
the basis polynomial and its root index: boundaries denote those **root
functions over the whole base cell**, not fixed numerical sample values.
In Lazard mode the root functions are those of the certified Lazard evaluation.

Following parent indices recovers a sample tuple. The geometric dimension is
the number of sector positions on that path. Terminal cells contain certified
evaluations and signs for the original input family. Projection theorems,
complete root lists, order/equality proofs and stack coverage establish that
these signs hold throughout each cell. Sample-point signs alone are insufficient.

`cad_cell` lifts only the path containing an exact query point. It uses the same
global projection closure and complete local stacks; it does not build unrelated
branches. Its checker also verifies point containment at every coordinate.

`cad_covering` decides a conjunction of polynomial sign conditions. Masks use
bit 1 for negative, 2 for zero and 4 for positive; unions express non-strict
inequalities or disequality. A nonnegative `witness_cell` identifies a satisfying
terminal cell. Otherwise every terminal cell names a violated constraint, and
the certified decomposition covers all real space. This implementation uses a
complete CAD for coverings.

`cylindrical_covering` provides the selective alternative. It retains a global
projection closure, samples uncovered intervals (zero first, then exterior
sectors and gaps), and stops at a satisfying assignment. A violated constraint
excludes the interval between its own factor roots, crossing boundaries of
unrelated polynomials. Such intervals may overlap. Constraints decided at an
earlier coordinate exclude entire cylinders, avoiding later lifting.

`CylindricalCoveringCertificate` stores only visited cells and local stacks.
Conflict intervals refer to inclusive ranges of alternating stack positions;
their parity encodes open or closed finite endpoints. Input-factor identities
and root mappings justify crossing irrelevant boundaries. For infeasibility,
the checker combines intervals and completed children from leaves to root and
requires every parent stack position to be covered. For feasibility, it checks
the complete assignment and all constraint signs. Constant contradictions need
no stack. No producer or search is invoked by this checker.

`coverings_solve` decides a conjunction without any global projection, by the
conflict-driven search over cylindrical algebraic coverings of Ábrahám,
Davenport, England and Kremer (2021). It assigns one coordinate at a time and
stops at the first satisfying point. A constraint is delineated only once it
fails at a candidate sample; a level whose samples all fail returns intervals,
each with the polynomials characterizing it. Their Lazard projection (leading
and trailing coefficients, discriminants, and resultants with interval bounds
on the sides where a polynomial has roots, and of adjacent bounds) gives the
interval excluded one level down. With `lazard=False` the characterization
uses every coefficient and McCallum nullification raises `UNSUPPORTED`.

`CADCoveringsCertificate` holds either a witness, with exact evaluations and
signs of every input, or a covering tree. Each `CoveringNode` covers one
coordinate's line above the prefix formed by its ancestors' samples, by
`CoveringInterval`s in increasing order. A `CoveringBound` names a root of a
listed polynomial at that prefix by its index. An interval is a single section
or an open interval; it is excluded either by a `constraint` whose sign at the
interval sample fails, or by a `child` covering above the sample. Supplied
factorizations of the inputs and of the projection polynomials are checked by
product. `check_coverings` recomputes every bound, requires the intervals to
cover the line, finds no characterization root inside an open interval, checks
conflict signs exactly, and for every child recomputes its Lazard projection
(requiring nonzero discriminants) and requires each factor in the parent
interval's characterization. It reuses the producer's exact root and sign
routines but none of its choices. `CADFiber.search` and its `check` overload
start the search and the check at the fiber's fixed prefix. `poly.cad.solve`
returns this evidence.

All traversals use explicit worklists and budgets. Invalid certificates return
false or a documented typed error; exhausted checking raises `LIMIT`.

<!-- example: docs/examples/cad.mojo -->

<!-- api: src/polynomials/cad_types.mojo CADCell -->

<!-- api: src/polynomials/cad_types.mojo CADCellCertificate -->

<!-- api: src/polynomials/cad_types.mojo CADCertificate -->

<!-- api: src/polynomials/cad_types.mojo CADCoveringCertificate -->

<!-- api: src/polynomials/cad_types.mojo CADPolynomialLift -->

<!-- api: src/polynomials/cad_projection.mojo CADProjectionCertificate -->

<!-- api: src/polynomials/cad_types.mojo CADSignWitness -->

<!-- api: src/polynomials/cad_types.mojo CADStackCertificate -->

<!-- api: src/polynomials/determinant.mojo DeterminantStep -->

<!-- api: src/polynomials/lazard.mojo LazardEvaluationCertificate -->

<!-- api: src/polynomials/determinant.mojo ProjectionResultantCertificate -->

<!-- api: src/polynomials/cad.mojo cad -->

<!-- api: src/polynomials/cad.mojo cad_cell -->

<!-- api: src/polynomials/cad_covering.mojo cad_covering -->

<!-- api: src/polynomials/cad_projection.mojo cad_projection -->

<!-- api: src/polynomials/cad.mojo check_cad -->

<!-- api: src/polynomials/cad.mojo check_cad_cell -->

<!-- api: src/polynomials/cad_covering.mojo check_cad_covering -->

<!-- api: src/polynomials/cad_projection.mojo check_cad_projection -->

<!-- api: src/polynomials/lazard.mojo check_lazard_evaluation -->

<!-- api: src/polynomials/determinant.mojo check_projection_resultant -->

<!-- api: src/polynomials/lazard.mojo lazard_evaluate -->

<!-- api: src/polynomials/determinant.mojo projection_resultant -->

<!-- api: src/polynomials/cylindrical_covering.mojo CADConflict -->

<!-- api: src/polynomials/cylindrical_covering.mojo CADConflictInterval -->

<!-- api: src/polynomials/cylindrical_covering.mojo CylindricalCoveringCertificate -->

<!-- api: src/polynomials/cylindrical_covering.mojo cylindrical_covering -->

<!-- api: src/polynomials/cylindrical_covering.mojo check_cylindrical_covering -->

<!-- api: src/polynomials/cad_coverings_search.mojo CoveringBound -->

<!-- api: src/polynomials/cad_coverings_search.mojo CoveringInterval -->

<!-- api: src/polynomials/cad_coverings_search.mojo CoveringNode -->

<!-- api: src/polynomials/cad_coverings_search.mojo CADCoveringsCertificate -->

<!-- api: src/polynomials/cad_coverings_search.mojo coverings_solve -->

<!-- api: src/polynomials/cad_coverings_search.mojo check_coverings -->
