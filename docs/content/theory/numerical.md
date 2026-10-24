# Numerical polynomial guarantees

Numerical polynomial operations use the public scalar operations of the pinned
apn_mojo dependency. Rounded operations follow their supplied contexts. Interval operations
enclose the set of results. The exact polynomial traits make stronger algebraic
promises and remain separate.

## Representation and enclosure evaluation

**Definition N-1.** A coefficient rectangle is the product of two closed rational
intervals. A real interval has imaginary part exactly zero. A polynomial family
chooses one coefficient in each rectangle independently, holding the choices
fixed for the whole evaluation. A singleton family denotes one polynomial.

**Proposition N-2.** Interval Horner evaluation encloses every value of every
member of the input family throughout the argument rectangle. Multiplication
of real intervals takes extrema of their four endpoint products; complex
multiplication uses real and imaginary component formulas. Induction over the
Horner steps proves containment. Sparse evaluation follows the same induction
over powers, terms and sums. Repeated arguments can introduce excess width but
cannot invalidate containment. A checker verifies each rational bound without
trusting the floating-point computation that proposed it.

### Incremental witnesses and exact endpoint extraction

**Proposition N-7.** The instruction graph fixes the witness order. Before
checking instruction $i$, retain every accepted predecessor bound with a future
use. Checking containment against exact rational input, addition or product
bounds proves the same induction as N-2. A predecessor slot can be released
after its last use. The final bound remains live until the completed scalar
receipt is checked. Missing, extra or invalid steps cannot establish a complete
proof; an invalid step makes the checker reject subsequent completion.
The producer and checker construct their own graphs and storage schedules.
Neither a receipt alone nor a successfully checked prefix is a certificate.

Streaming removes the retained list of all witness boxes. It still retains
the instruction graph, input rectangles and simultaneously live values or
bounds; its total storage is not constant in polynomial size.

**Proposition N-8.** Write each finite binary midpoint and radius as an integer
times a power of two. Let $L$ be the smaller least-significant exponent and
$H$ the larger most-significant exponent bound. Both operands are multiples
of $2^L$ with magnitude less than $2^H$. Their sum and difference fit in
$H-L+1$ bits, allowing the carry. Computing them in that binary format is
exact, including cancellation. The producer traps inexactness, underflow and
overflow and converts the results to exact rationals. Zero midpoint and zero
radius have direct exact conversions. The independent checker instead obtains
midpoint and radius rationals and computes their sum and difference using
rational arithmetic; it does not rely on the producer's binary operations.

**Proposition N-9.** For a fixed assignment and context, binary exponentiation
uses the same deterministic ladder of repeated squares. Reusing an identical
square preserves its rounded scalar and exact enclosure, provided each power's
remaining multiplication sequence is unchanged. The sparse producer records
aliases only to earlier canonical square instructions, extends their storage
lifetimes and emits every original witness position. The checker verifies the
original instructions independently of these aliases. The direct evaluator's
square cache is scoped to a single call and keyed by assignment and ladder bit.

## Real root predicates

**Proposition N-3.** An evaluation enclosure excluding zero proves root absence.
Opposite strict endpoint signs prove existence by continuity. A derivative
interval excluding zero proves strict monotonicity and hence uniqueness. These
claims hold for each polynomial in the family when all the bounds do.

**Proposition N-4.** Let I be a closed interval with midpoint m and let D enclose
every derivative value on I, with zero outside D. Every root in I lies in
N = m - F(m)/D, by the mean value theorem. Thus a disjoint Newton interval
excludes roots. Strict inclusion N inside I establishes existence and uniqueness:
for each fixed polynomial the endpoint values have opposite signs, using the
mean value theorem between m and each endpoint and the same strict inequalities.
Root coverage additionally requires that no subdivision boundary root is lost.
Exact inputs instead use Descartes certificates on the squarefree part (see the
[real roots reference](../reference/p5.md#descartes-certificates)).

## Complex root counts

**Proposition N-5.** Traverse a rational rectangle counterclockwise. Subdivide
its boundary into segments. Suppose each segment image for every permitted
polynomial lies in a convex rectangle excluding zero. Choose one fixed rational
representative polynomial from the coefficient family and join its values at
consecutive boundary vertices with straight lines. Both the polynomial image
and this polygon lie in the same zero-free convex image rectangle on each
segment. Their straight-line homotopy avoids zero and agrees at the shared
vertices. Consequently their winding numbers about zero agree.

Before enclosing a segment, polynomial Horner translation expands about its
midpoint using rational interval operations. The translated coefficient family
contains every translated original polynomial; subsequent interval evaluation
therefore preserves the homotopy premise despite any lost correlations.

The argument principle identifies this winding number with the number of roots
inside the domain rectangle, counting multiplicity. The polygon's winding number is
computed using exact rational oriented crossings and a half-open convention
at vertices. The checker must verify the entire closed boundary, image bounds,
zero exclusion, representative membership and winding count.

**Proposition N-6.** A positive lower bound L on the leading coefficient's
magnitude and upper bounds U on the others give a uniform Cauchy radius
1 + max(U/L). Disjoint certified regions whose root counts sum to the fixed
degree cover all roots, counting multiplicity. A multiple count certifies a
cluster, not coincidence of its roots. Families whose degree or relevant root
counts vary may have no requested uniform isolation; inconclusive is required.

## Resource and accuracy limits

Precision refinement reconstructs values from the original exact inputs or
preserves their original uncertainty. Increasing precision alone cannot remove
uncertainty or correlations lost by interval evaluation. Every search has an
explicit precision/subdivision ceiling and a shared charged work budget.
No claimed count or accuracy is inferred from midpoint values, small residuals,
agreement of two numerical runs, or an exhausted search.
