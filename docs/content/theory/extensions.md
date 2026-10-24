# Algebraic and CAD extensions

This page states the mathematical obligations of algebraic and CAD
certificates. The [references](references.md) provide further background.

## Exact algebraic values

**Definition A-1 (value witness).** A nonconstant monic polynomial over QQ,
an open rational interval with nonroot endpoints of opposite sign, and
Descartes bound one on that interval identify exactly one simple real root.
A linear polynomial with equal rational endpoints identifies an exact rational
root. Irreducibility is a producer normalization, not a checker assumption.

**Theorem A-2 (comparison certificate).** Ordered disjoint isolating intervals
prove strict comparison. If a nonzero common polynomial divides both defining
polynomials and has a root in their interval intersection, uniqueness forces
both values to be that root. When one operand is rational, evaluating the other
definition at that rational point and checking interval containment suffices.

Comparison production checks interval separation before polynomial gcd work.
For overlapping nonrational values, the gcd has either zero or one root in the
intersection, since each input interval contains exactly one root of its own
definition. Equality evidence retains exact common-factor identities and an
interval inside both refined witnesses where that factor changes sign and has
Descartes bound one. The checker verifies those identities and the root bound
without computing a gcd or searching for roots. Nonminimal definitions are
accepted when the selected root is simple.

Grouping specialization terms by their unassigned monomial is distributivity;
restoring first-occurrence order preserves the coefficient evaluation and
certificate schedule.

**Theorem A-3 (arithmetic certificate).** Let alpha and beta denote the operand
roots. For addition, eliminate t between `p(z-t)` and `q(t)`; for subtraction,
use `p(z+t)`; for multiplication, use the homogenization `t^deg(p) p(z/t)`.
An identity `h(z)=u(z,t) P(z,t)+v(z,t) Q(z,t)` proves that the nonzero h
vanishes at the intended result. A unique root of h in an interval containing
both the exact interval image and the proposed result proves the answer.
Negation changes alternating coefficient signs; inversion reverses the
coefficient list and requires a nonzero operand interval. These identities
give direct unary certificates.

The arithmetic checker needs a nonzero annihilator, not the exact resultant
determinant. The producer computes the resultant to find one. Cofactor
identities bind any normalized operand or result definition to its original
polynomial; input irreducibility is never trusted.

**Theorem A-4 (complete roots at algebraic coefficients).** Give each nonrational
coefficient its own variable and eliminate these variables against their
defining polynomials. Every root of the original polynomial is a root of the
resulting nonzero rational annihilator. Complete Descartes isolation of its
squarefree part inside a certified root bound, together with the binding
factor identity, accounts for all distinct real roots of the annihilator.
Certified evaluation of the original polynomial at every candidate selects
precisely the desired roots, even when
the coefficients are algebraically dependent. No candidate can be omitted.

The factory uses minimal coefficient definitions so that a nonzero leading
coefficient has no zero conjugate; consequently its norm polynomial is nonzero.
The checker only needs the final nonzero annihilator and the elimination
identities. Constant input is valid; identically zero input is rejected.

The AFP `Algebraic_Numbers` entry for Isabelle2025-2 supplies the specification
reference: `Sturm_Rat.root_cond` and `root_info_cond`, and
`Real_Algebraic_Numbers.add_3`, `mult_3`, `inverse_3` and `compare_3`.
Its root interval convention is closed; this package's nonroot open endpoints
represent the same isolated root. Isabelle totalizes inverse at zero; this
package deliberately raises `DOMAIN` at that public arithmetic boundary.

## Projection and lifting

**Theorem A-5 (determinant transcript).** In a fraction-free elimination step,
after an optional row swap, let a be the current first pivot and d the previous
pivot (initially one). Every next trailing entry B satisfies
`d*B[i,j] = a*A[i+1,j+1] - A[i+1,0]*A[0,j+1]`.
Sylvester's determinant identity proves inductively that the last entry is
the determinant, adjusted for row-swap parity. A zero first column proves
singularity. The checker verifies these multiplication identities, nonzero
pivots, shapes and termination, without division or pivot search.

The implementation maps row indices through the pivot transposition instead
of swapping or copying matrices. The previous transcript matrix supplies the
next step's input by an immutable borrow. This changes only storage ownership:
the producer computes exact quotients, while the checker independently verifies
the identities above against every supplied trailing entry. Neither keeps a
second deep copy of the trailing matrices.

**Definition A-6 (projection basis).** At each main variable, the basis consists
of normalized, primitive, squarefree, pairwise coprime nonconstant polynomials.
Product identities reduce all original and projected polynomials to these
bases. Coprimality of f with its main-variable derivative certifies both
primitivity and squarefreeness; pairwise content witnesses certify coprimality.

**Theorem A-7 (projection/lifting invariant).** McCallum's coefficients,
discriminants and pairwise resultants give delineability and sign/order
invariance under the well-orientedness hypotheses. Lazard's leading and
trailing coefficients, discriminants and pairwise resultants give the
corresponding valuation-invariant lifting theorem using Lazard evaluation.
Replacing a discriminant with `Res(f,f')` changes projection only by leading
coefficient factors already included.

The checker verifies the entire projection closure and exact determinant
values, then inducts through the cell tree. Each stack supplies complete root
lists, all merged-root bindings and strict order, sections at roots, and sectors
between consecutive roots and at both infinities. Thus it partitions each
cylinder completely. Input product identities transfer invariant basis signs
to invariant input signs; certified sample evaluations determine their values.

McCallum nullification is rejected at all intermediate levels and over
positive-dimensional bases. A nullified final-level polynomial above a point
is identically zero on that one cylinder. This deliberately conservative
condition avoids an unsupported claim of order invariance at intermediate
nullification; Lazard is available without that restriction.

**Lemma A-8 (Taylor implementation of Lazard evaluation).** If k is the first
nonzero specialization of successive derivatives in x at a, the original
polynomial is divisible by `(x-a)^k` in the specialized coefficient field.
The k-th derivative at a equals k! times the quotient evaluated at a.
Repeating in cylindrical variable order gives Lazard evaluation up to positive
factorials. Every earlier derivative is checked to vanish identically as a
polynomial in the unassigned variables, rather than merely at one full point.

**Corollary A-9 (cells and coverings).** Choosing one section or sector at each
level gives a connected sign-invariant cylindrical cell. A conjunction is
satisfied on a terminal cell precisely when every certified sign belongs to
its allowed mask. If every terminal cell violates some conjunct, complete
stack coverage proves infeasibility over all real coordinate space.

**Corollary A-9a (selective overlapping coverings).** Over a certified base
cell, an input polynomial is a product of its certified projection factors.
Between consecutive roots of those factors its sign cannot change; roots of
unrelated factors may be crossed. Identically vanishing specializations remain
zero. Thus a conflict at a stack sample excludes every consecutive stack
position up to the nearest relevant factor roots. At a relevant root the
section alone is used. A zero allowed-sign mask excludes the whole cylinder.
Constraints using only assigned variables exclude every extension of that
base cell. Verified conflict intervals and recursively excluded children may
overlap; a bottom-up union covering every stack position proves infeasibility.
A fully checked satisfying assignment suffices for feasibility, so unrelated
branches need not be constructed. Global projection remains the delineability
certificate; only lifting and conflict sampling are selective.

## Further factor and gcd algorithms

**Lemma A-10 (monic transformation).** For a primitive polynomial f in x of
degree n with leading coefficient a in QQ[Y],
`F(X)=a^(n-1) f(X/a)` belongs to QQ[Y][X] and is monic. Factoring F, substituting
X=a*x, and removing parameter content recovers the irreducible factors of f
by Gauss's lemma and the invertible change over QQ(Y).

The implementation uses deterministic integer evaluation points with
squarefree specializations, factors the specialized univariate polynomial,
and lifts through increasing total parameter degree. Complement inverses
modulo the specialized coprime factors solve each homogeneous correction.
The producer compares three distinct good evaluation points, preferring fewer
specialized factors and then smaller maximum coefficient magnitude. A single
irreducible monic specialization immediately proves irreducibility globally.
Otherwise each true monic factor corresponds uniquely to a subset of the
coprime specialized factors. Total parameter degree is additive under products,
so it bounds every true factor's degree. Truncated homogeneous lifting to this
bound reconstructs each true subset exactly. Subsets are tested in increasing
cardinality, and exact divisibility accepts the first genuine factors; all
proper subgroups have already been excluded. Exhausting subsets up to half the
remaining cardinality proves that the remaining factor is irreducible. Hence
split specializations are recombined, without depending on a later Hilbert
specialization. The final product is verified exactly.

**Lemma A-11 (sparse reconstruction).** After separating main-variable content,
a good modular specialization gives an upper bound on the main degree of the
characteristic-zero gcd. Sparse interpolation uses a previously discovered
monomial support and equations equating each specialized coefficient to the
corresponding monic univariate-gcd coefficient times the unknown leading
coefficient. A full-rank solution and modular trial division validate the
interpolant. Fixed-prime CRT lifts coefficients scaled by the gcd of the input
leading coefficients. An exact common divisor with the modular degree bound
is the full primitive gcd. Restoring content and checking the coprimality cofactor
certificate supplies independent maximality evidence.

Unlucky degrees, support changes and singular interpolation are rejected or
rediscovered; they never become negative gcd answers. Support discovery follows
Z3's recursive modular construction. Treat the last lexicographic variable t
as the parameter, remove content in Fp[t], and compute the univariate gcd gamma
of the primitive inputs' leading coefficients in the other variables. Only
samples preserving both leading coefficients are accepted. Recursion removes
t; multiplying each monic image gcd by gamma(t) fixes its scale for sparse
polynomial Newton interpolation. A smaller leading monomial resets prior
samples. Remove the interpolant's Fp[t] content and require exact division of
both primitive inputs and the minimal observed leading monomial. The gcd's
leading monomial survives each good specialization, so these conditions force
maximality. Field exhaustion rejects the prime. Integer content uses the same
sparse gcd recursively on coefficient polynomials with one fewer variable;
only occurring exponents are grouped. Coprimality certificate generation remains an
independent, existing proof-producing path.

At each recursive variable set, the scaled interpolant's monomial support is
retained as a hint. Future samples solve its coefficient equations over the
finite field; the same leading-coefficient, content, degree and exact-division
checks apply before acceptance. Failed or singular support interpolation
continues the independently accumulated Newton interpolant. Cache state cannot
justify an answer without those checks, including when support changes with
a prime or specialization.

**Theorem A-12 (trace lattice recombination).** For each p-adic factor, Newton
sums supply scaled root traces. True rational-factor indicators have bounded
integer traces. Removing low p-adic digits above that bound puts each true
indicator into a short vector of the trace-knapsack lattice. Exact LLL preserves
the integer lattice. Dropping only a suffix whose Gram–Schmidt squared norms
exceed the short-vector bound preserves every true indicator in the remaining
span. If the projected row space has a partition basis and all reconstructed
groups divide f, that partition equals the irreducible rational factorization.
Otherwise trace information or precision is increased under the budget.

This is van Hoeij's 2002 trace-knapsack argument, sections 2.1–2.2. The
implementation uses exact rational Gram–Schmidt calculations, delta=3/4,
incremental exact reduction/swap identities, and integer Hensel lifting.
No floating-point decisions or external lattice library are used.

**Lemma A-13 (scoped reuse).** Hensel correction inverses depend only on the
fixed input and its prime factors, so increasing precision can resume from a
previous modulus and lifted factors. After exact singleton extraction the
remaining lifted residues are retained and inverses are recomputed for the
new input. Algebraic normalization caches factorizations and Sturm sequences
by exact defining polynomial; each requested interval still receives its own count and
factor selection. Thus conjugates are never identified by a cache key. Caches
are local to one calculation, budget their storage and lookup work, and do not
alter certificate relations or checker behavior.

That scope extends to all polynomial sign evaluations at one CAD point.
Projection also reuses factorizations under exact polynomial, order and ring
keys. Each scheduled certificate occurrence receives owned factors, including
repeated and zero coefficients in the McCallum schedule. Lazard projection
extracts its required trailing and leading coefficients directly. Cache hits
never omit a checker obligation or permit one returned witness to mutate another.
