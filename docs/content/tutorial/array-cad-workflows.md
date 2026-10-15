# Select the CAD work your calculation needs

Use `cad.project` to compute projection polynomials, `cad.cell_at` to construct
the cell containing an exact point, and `cad.lazard_evaluate` to prepare one
lifting polynomial. Each returns retained evidence for `cad.check` and accepts
ordinary sparse coefficient arrays.

<!-- example: docs/examples/array_parity_cad.mojo -->

## Keep axes and coordinates distinct

`axes=[1, 0]` changes the base-to-top CAD traversal order. Input and output
coordinates still use the polynomial's original axis order: `(2, -3)` assigns
$x=2$ and $y=-3$. Lists accept homogeneous exact values; coordinate tuples can
mix integers, rationals, and `algebraic.Real` values.

`cad.signs(cell)` reports the constant sign of each input polynomial on the
selected cell. `cad.sample(cell)` returns the exact query point. The retained
certificate establishes both sign invariance and query-point containment.
Use `cad.projection_polynomials(projection, axis=0)` to copy the projection
basis for an original axis, or omit `axis` to collect every level.

## Fix coordinates before decomposing

Pass `fixed=[(axis, value), ...]` to study a slice of your problem. For example,
fixing $y=\sqrt{2}$ in $x-y$ leaves three cells along the $x$ direction:
$x<\sqrt{2}$, $x=\sqrt{2}$, and $x>\sqrt{2}$.

<!-- example: docs/examples/array_fixed_cad.mojo -->

The same `fixed=` argument works with `decompose`, `solve`, `cover`,
`decompose_to`, and `checker`. Values can be integers, rationals, or exact
real algebraic numbers. Use a list for one value type, or a tuple such as
`fixed=((0, 1), (2, a))` to mix types. Use each selected axis once; axis
numbers must fit the polynomial's power shape. An empty iterable or `fixed=()` leaves every
original coordinate free.

By default, CAD traverses the fixed axes in increasing order, then the free
axes in increasing order. If you supply `axes=`, include every original axis
and put the fixed ones first. For example, fixing axis 1 of a two-variable
polynomial requires `axes=[1, 0]`. The order of pairs in `fixed` does not change
this traversal. A streaming checker takes the same `fixed=` and `axes=` as
its producer.

Samples retain **all original coordinates**, including the fixed values, in
original axis order. Signs still correspond only to your input polynomials.
Cell dimensions count the remaining free directions. Fixing every axis gives
one point; `solve` and `cover` decide whether that point satisfies the
constraints. A negative answer proves infeasibility on the requested slice.
Use `cad.fixed_coordinates(result)` to copy the retained bindings.

Fixed coordinates also work with algebraic coefficients. Lifting visits only
the requested coefficient values and coordinate values, while projection and
local root isolation still process the expanded rational problem. Shared
budgets include preparation, traversal, and checking. Stream completion requires
an accepted receipt; a partial stream after an error does not prove coverage.

Checked archives preserve the exact bindings and their root witnesses through
`io.dumps`/`loads` and incremental `io.dump`/`load`. Restoring a result checks
both the geometry and its coordinate scope. Equivalent exact values can use
different defining witnesses in a stream producer and checker.

## Lift through a vanishing specialization

The polynomial $x^2(y-1)$ becomes zero when you substitute $x=0$. Lazard
evaluation removes the vanishing factor before lifting in $y$. It records the
least nonvanishing derivative order, two in this example, and produces
$2(y-1)$. The positive factorial factor preserves the roots and signs used by
lifting.

Assignments are ordered `(axis, exact_value)` pairs. Leave the lifting axis
unassigned and assign every other occurring axis. `cad.vanishing_orders` keeps
that assignment order. `cad.lifting_coefficients` returns ascending exact-real
coefficients, which `dense.real_roots` can consume directly.

## Choose a feasibility workflow

`cad.solve` uses selective cylindrical covering to search a conjunction of sign
constraints. `cad.cover` builds a full CAD decomposition and retains either a
satisfying cell or a covering of conflicts. Use the latter when you need that
complete decomposition evidence; `cad.cells(result)` returns all terminal cells
with their dimensions, samples, and signs. This can require substantially more
work and memory. Both use the same constraints, `satisfiable()`, `cad.sample`, and
`cad.check` conventions. `cad.decompose` remains the direct operation for all
sign-invariant cells.

These operations accept integer, rational, and real algebraic coefficients.
Query points and assignments can also contain exact algebraic values. Complex
and modular coefficients are rejected. Lazard projection is the default; `lazard=False` selects McCallum projection. All preparation and
checking share a supplied `Budget`, or use `max_steps=` for a single call.
See the [CAD reference](../reference/array-cad.md) for signatures.


## Put algebraic values in coefficients

Use a `Real` value wherever a coefficient belongs. For example, the same
constraint helpers can express $x=\sqrt{2}$ and $x>0$:

<!-- example: docs/examples/array_algebraic_cad.mojo -->

Integer, rational, and algebraic constraints can share one list. A homogeneous
sparse polynomial uses one coefficient type; `sparse.asarray(((power, value), ...))` prepares a tuple with mixed
coefficient types. All existing workflows accept `AA`, including projection,
cell queries, Lazard evaluation, selective solving, full covering, and streaming.
`projection_polynomials` returns AA arrays when its input coefficients are AA.

Samples and signs keep their original meaning. Samples contain only your
polynomial axes, in their original order; signs contain one entry per input
polynomial. `cells` returns only the cells for the specified coefficient values.
An algebraic coefficient's other conjugates do not introduce extra user cells.
`vanishing_orders` likewise lists only the assignments you supplied.

Decomposition follows only the coefficient values you supplied. For example,
$x-\sqrt{2}$ gives three cells on the real line: below, at, and above
$\sqrt{2}$. Lifting does not explore the coefficient's negative conjugate or
intervals between coefficient roots. This applies to `decompose`, `cover`,
`solve`, and `decompose_to` automatically.

The implementation still builds a rational projection problem with auxiliary
variables for distinct irrational coefficients. Rational values need no extra
axis; equal or opposite irrational values share one. Projection and exact root
arithmetic can remain expensive, especially with many distinct coefficients.
Use a work budget, and choose `solve`, `project`, or `cell_at` when a full
decomposition would give you more than you need.

## Check or consume algebraic CAD evidence

Use `cad.check(result)` to check the geometric certificate and its exact
connection to your coefficients and fixed coordinates. `.evidence()` borrows
the underlying rational certificate, including auxiliary variables and
isolating constraints. Its tree covers the selected coefficient and coordinate
values, then every remaining free polynomial axis. The ordinary core
`check_cad` requires full-space coverage; use the facade checker for this
restricted claim.

`decompose_to` streams the restricted traversal without retaining its full cell
history. Prepare its checker with the same polynomials, `fixed=`, and `axes=`,
then finish with the returned receipt. For AA inputs, the checker also enforces
the auxiliary-first axis order and selected coefficient sections.

Custom sinks receive original input signs first, followed by triples for each
coefficient's defining polynomial and two strict isolating bounds. Those
triples have signs `(0, 1, 1)`. Any user-fixed axes add one defining-polynomial
sign each, in traversal order; these signs are zero. Auxiliary variable IDs
begin at the original axis count, and the projection header places them before
the original axes.

A stack event contains every local root and its proof, but only the selected
section has descendants on a fixed axis. Receipt cell counts equal the number
of terminal cells returned by `cad.cells`; stack counts include both coefficient
and user-fixed prefixes. For $x-\sqrt{2}$ with its original axis left free,
the receipt records two stacks and three cells.
