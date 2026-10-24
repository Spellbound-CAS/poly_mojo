# Find roots from coefficient arrays

Pass coefficients to `dense.real_roots` or `dense.complex_roots` to get exact
algebraic numbers you can calculate with. Use `dense.real_root_cover` or
`dense.complex_root_cover` when you want enclosures with an explicit accuracy
target. The matching `sparse` functions accept `(power, coefficient)` terms.

## Calculate with an exact root

The coefficients `[-2, 0, 1]` represent $x^2-2$. Its two real roots are
irrational, but their squares are exactly 2.

<!-- example: docs/examples/array_solver_roots.mojo -->

Exact roots are sorted in increasing order, and each distinct root appears
once. Rational and irrational roots use the same `algebraic.Real` value type.
Index the result or iterate over it; neither operation copies the full proof.
`roots.check` verifies that the retained certificate accounts for every real
root of the original polynomial.

Use `dense.root_between(p, low, high)` when you need one selected root. An open
interval must contain exactly one distinct root and have nonroot endpoints.
Equal endpoints can select an exact root, including an algebraic endpoint. `dense.count_roots` uses
ordered nonroot endpoints and returns a result whose `count()` gives the number
of distinct roots and whose evidence can be checked with `roots.check`.

Real algebraic coefficients and endpoints work with the same calls.
`sign_on_interval(p, low, high)` returns a sign only when the closed interval
is root-free; either endpoint being a root returns `None`. Integer coefficients
promote to rationals. Nonzero
constants have no roots; the zero polynomial raises `DOMAIN` because its root
set is not finite. Sparse inputs must have a one-variable power shape and are
converted to dense storage under the same work budget.

## Reuse a Sturm chain

If you need root counts in many intervals for one rational polynomial, prepare
its Sturm chain once. `dense.sturm_sequence` accepts integer or rational
coefficients and returns a shared immutable plan. Each subsequent query
evaluates the chain at its endpoints without rebuilding it.

<!-- example: docs/examples/array_access_sturm.mojo -->

For $x^2-2$, the chain is $x^2-2$, $2x$, $2$. At zero its signs are negative,
zero, positive; `roots.sign_changes(plan, 0)` ignores the zero and returns one.
`roots.polynomials(plan)` returns the chain as ordinary rational coefficient
arrays when you want to inspect it.

`roots.count_roots(plan, low, high)` counts distinct roots in the open interval;
endpoints must be ordered integers or rationals and must not be roots. Repeated
roots count once. The result shares the plan, retains both endpoints, and can be
verified with `roots.check`. Checking verifies the recurrence and count
independently; it does not repeat the polynomial divisions used to build the
chain. `roots.check(plan)` checks the chain itself.

Use `sparse.sturm_sequence(terms, axis=...)` for sparse input with one active
axis. Other axes may remain in the power shape but must have zero exponents.
Sparse preparation charges the full dense extent before allocating it.
Nonzero constants have no roots; zero raises `DOMAIN`. For real algebraic
coefficients or endpoints, use the direct `dense.count_roots` or
`sparse.count_roots` workflows instead of a rational Sturm plan.

Plans and count results support string and incremental [serialization](array-serialization.md).
The archive retains the quotient witnesses and validates the complete chain on
load. All preparation, queries, inspection, checking, and IO accept a shared
`Budget` or `max_steps=`.

## Find all exact complex roots

`dense.complex_roots(p)` returns every distinct complex root as an
`algebraic.Complex`. Integer, rational, Gaussian-rational, `AA`, and `QQbar`
coefficients are accepted. You don't choose a working precision: each result
has exact real and imaginary algebraic components.

<!-- example: docs/examples/array_complex_roots.mojo -->

The result orders roots by real part, then by imaginary part. This is a
consistent presentation order, not an ordering of the complex field. A repeated
root appears once; `result.multiplicity(i)` gives its multiplicity. The sum of
these multiplicities equals the input degree. Nonzero constants return an empty
result; zero raises `DOMAIN` because it has infinitely many roots.

Index or iterate over the result to reuse its roots as ordinary coefficients,
call `dense.from_roots(result)` for their monic squarefree product, or evaluate
another polynomial at them. `dense.from_roots` consumes each distinct value
once; expand according to the multiplicities when reconstructing a repeated
factor. For approximations, use `algebraic.bounds`, `refine`, or `approximate`
on `root.real()` and `root.imag()` separately.

`sparse.complex_roots(terms, axis=0)` accepts a selected axis in any power shape,
provided no other variable remains. For a parameterized polynomial, specialize
its parameters first. Conversion to dense storage is charged to the same budget.

`roots.check(result)` verifies the component witnesses, strict presentation
order, positive multiplicities, and the full polynomial reconstruction.
`roots.check_complex(p, result.evidence())` makes the same check against an
independently supplied dense coefficient sequence. It rejects missing roots,
duplicates, incorrect multiplicities, and values that do not reconstruct the
supplied polynomial.
Both routines accept `Budget` or `max_steps=`. Borrowed evidence is immutable;
copy it before editing or handing it to code that needs mutable storage.

The solver uses squarefree decomposition, real-coordinate elimination, and
exact real-root isolation. These calculations can grow substantially for high
degrees or unrelated algebraic coefficients. Budget exhaustion raises `LIMIT`;
the call never presents an incomplete list as a complete result. Numerical root
covers remain useful when you need certified regions rather than algebraic
values. Use `complex_roots_to` when you want to consume exact roots as they arrive;
`complex_roots` retains and sorts the complete result.

## Keep arithmetic exact until you choose an approximation

`algebraic.add`, `subtract`, `multiply`, `divide`, `negative`, and `reciprocal`
accept exact real values, integers, and rationals. They return new values and
retain evidence for their arithmetic steps. `algebraic.check` verifies the
retained step; division retains both its inverse and multiplication proofs.
`algebraic.compare` returns -1, 0, or 1 using exact ordering.

`algebraic.bounds(value)` returns rational endpoints. They agree for rational
values. `algebraic.refine(value, width)` returns the same exact value with a
narrower root enclosure, leaving the input unchanged. The refined result keeps
a root witness; the earlier arithmetic derivation stays with the original
value. `algebraic.approximate(value, width)` returns the enclosure's rational
midpoint with absolute error below `width / 2`. This approximation does not
silently convert the result to a machine float.

A `Real` is a selected exact root and a coefficient in the inferred `AA` field.
Use it in ordinary dense or sparse arithmetic, evaluation, specialization, and
the [algebraic-coefficient root workflows](#continue-from-algebraic-coefficients-or-parameters)
below. See [algebraic coefficients](algebraic-coefficients.md) for exact complex
values and mixed-domain promotion.

## Choose the numerical claim you need

Root covers retain the input polynomial and the options used to produce them.
`roots.refine(cover, new_options)` keeps the original input with the new request;
you do not have to keep a separate copy of the original options.

| Result | Interpretation |
|---|---|
| `complete()` | Every root is accounted for |
| `accurate()` | Counted regions meet the requested width target |
| `isolated()` | Every counted region has count one |
| `outcome()` | Success or the reason a requested property remains unresolved |

`roots.intervals` returns `(low, high, count)` tuples. For a real cover of
integer or rational coefficients there is one `(low, high, 1)` per distinct
root, increasing: an open interval with the root strictly inside, or
`low == high` for an exact rational root; unresolved pieces keep count -1.
Real covers of `AA` coefficients and numerical arrays list cells, omitting
zero-count cells. Real counts are distinct-root counts.
`roots.rectangles` returns `(real_low, real_high, imag_low, imag_high, count)`;
complex counts include multiplicity. A repeated complex root therefore needs
`NumericOptions(require_isolation=False)` if a certified cluster is acceptable.

The default `complete_only=False` lets you inspect unresolved outcomes. Set it
to `True` to raise unless every requested property succeeds. A checker rejects
an incomplete cover, but can accept a complete cover whose requested accuracy
or isolation remains unresolved. Inspect the flags your application needs.

Real covers accept integers, rationals, and `AA` coefficients. Complex covers
also accept Gaussian rationals and `QQbar`. Covers accept prepared
[numerical arrays](array-numeric.md), including ball coefficient families. Complex
rational components enter certification exactly; algebraic coefficients retain
checked conversion witnesses. The
[numerical guarantees](../theory/numerical.md) explain the underlying enclosure
methods.

## Process witnesses as they arrive

Numerical root-cover `_to` functions feed a sink directly. They do not build a complete result and then iterate over it. The following program uses an
independent checker as the consumer:

<!-- example: docs/examples/array_solver_streams.mojo -->

Keep the receipt and call `finish` after successful production. A receipt alone
is not a proof, and events received before an exception do not establish a
complete result. Budget exhaustion always raises. Streaming avoids retaining
the full witness history; input storage, algorithm state, and checker state
still consume memory. Complex stream checkers also retain rectangle metadata
for disjointness checks.

See the [root result reference](../reference/array-roots.md) and
[algebraic value reference](../reference/array-algebraic.md) for the functional
interfaces and borrowed certificate access.


## Evaluate and specialize at exact roots

`dense.evaluate(p, x)` returns an exact `algebraic.Real` when `x` is algebraic.
Sparse evaluation also accepts a coordinate tuple mixing exact integers,
rationals, and algebraic values. `algebraic.check(value)` verifies the retained
evaluation chain as well as the returned value's root witness.

<!-- example: docs/examples/array_algebraic_points.mojo -->

`specialize` takes `(axis, value)` assignments and preserves the original axis
numbering. Its iterable result contains `(powers, algebraic_coefficient)`
groups in the core's traversal order, including groups that evaluate to zero.
`algebraic.check(result)` verifies every group. When one free axis remains,
`algebraic.coefficients(result, axis=...)` supplies an ascending sequence for
further univariate evaluation. Missing powers become exact zeros.

## Recover multiplicities and inspect a candidate

<!-- example: docs/examples/array_root_diagnostics.mojo -->

A real cover counts distinct roots. `roots.multiplicities(cover)` adds exact
multiplicities for an isolated complete cover of an integer or rational
polynomial. Its `counts()` align with `roots.intervals(cover)`. `roots.check`
verifies the squarefree product and that each root's block vanishes there.
Complex cover counts already include multiplicity, but a region containing
several roots does not identify each root's individual multiplicity.

`roots.diagnostics(p, point)` returns residual, derivative, coefficient
backward-error, and local sensitivity bounds. These are evaluated at the
specified point. Sensitivity is a root condition bound only when the point
is a root; an unavailable derivative-based bound is `None`.

`roots.polish` takes a complex arithmetic context, a positive exact `radius`,
and `NumericOptions`. Newton iteration proposes a candidate, then an independent
contour calculation attempts to prove that its rectangle contains exactly one
root. A good-looking candidate alone is not a certificate. Inspect `outcome()`
and use `roots.check(result)` to validate the optional one-root evidence;
the checker returns `False` when that evidence is absent. A checked enclosure
can still have `SUBDIVISION_LIMIT` status if its width misses the requested
accuracy. Diagnostics and polishing accept exact coefficient iterables and
prepared numerical dense arrays. Convert a univariate sparse input with
`poly.to_dense` first.


## Continue from algebraic coefficients or parameters

`dense.real_roots` also accepts a sequence of `algebraic.Real` coefficients.
Use `algebraic.asreal` for rational entries in that sequence. This lets an exact
value from an earlier calculation become a coefficient without approximation.

<!-- example: docs/examples/array_extension_algebraic_roots.mojo -->

Both equations have the exact root $\sqrt{2}$. `sparse.real_roots` accepts a
prepared specialization and selects the remaining root `axis=`. `sparse.roots_at`
combines specialization and root finding when you do not need the intermediate
coefficients. Assign every other occurring variable and leave the root axis
unassigned. A zero specialization raises `DOMAIN`; a nonzero constant has no roots.

The result is iterable and sorted, with repeated roots listed once. `roots.check`
verifies the specialization, binds its coefficients to the root problem, and
checks completeness of the returned root list. `specialization()` returns the
retained intermediate value and its coefficient proofs.

## Stream numerical root witnesses

Prepared numerical arrays also support `real_root_cover_to`, `complex_root_cover_to`, and their matching `*_checker` factories. A checker can serve directly as the producer's sink; validate the final receipt with `checker.finish(receipt, budget)`. Root options, family uncertainty, and completeness requirements have the same meaning as for retained covers.

Univariate sparse numerical inputs require `context=` for conversion to dense coefficients. That conversion shares the work budget and allocates storage up to the highest power. For very high sparse degrees, check whether the resulting dense representation fits your application before choosing this workflow.

## Diagnose an algebraic polynomial numerically

Keep an `AA` or `QQbar` polynomial exact and pass an explicit ball context to
`roots.diagnostics(p, point, context=BallContext(precision=128))`. The result
retains enclosures of every algebraic coefficient and their refinement
witnesses. `roots.check` verifies those conversions before checking the
residual and derivative bounds.

`roots.polish` accepts the same coefficient domains and uses the precision in
its `ComplexContext` to enclose them. It retains the conversion evidence with
the polishing result. A precision limit can still leave root certification
unresolved; inspect the result and call `roots.check` before treating it as
certified. More working precision may be needed for tightly clustered roots.

`sparse.roots_at` also accepts `AA` coefficients. Assign every parameter axis,
leave the requested root axis free, and use the returned exact roots as usual.
Both the rational and `AA` paths retain the specialization transcript and
complete root evidence. `roots.check(result)` binds the original coefficients
and assignments to the specialized polynomial and then to its roots.
`result.specialization()` gives you the checked intermediate terms.

<!-- example: docs/examples/array_workflow_specialization.mojo -->

Here the assignment changes $\sqrt{2}x+\sqrt{2}y$ into $\sqrt{2}x+2$.
The result keeps the reason its root is $-\sqrt{2}$, including after an archive
round trip. A standalone `sparse.specialize` result is iterable; pass it to
ordinary sparse operations or check it with `algebraic.check`.

## Cover algebraic coefficients directly

Pass an `AA` polynomial directly to `real_root_cover`, or an `AA`/`QQbar`
polynomial to `complex_root_cover`. The calls keep the exact coefficients and
the witnesses connecting them to enclosing balls. `roots.check` verifies that
connection before checking the numerical regions. The matching `_to` and
`_checker` calls support the same domains.

<!-- example: docs/examples/array_workflow_covers.mojo -->

The cover encloses the root of $x-\sqrt{2}$. Coefficients are enclosed at
`options.max_precision`, giving the producer and incremental checker the same
fixed coefficient family. Root predicates start at `initial_precision` and
respect the usual precision and subdivision limits. This preparation costs
work even when a coarse root enclosure would suffice; choose a reasonable
ceiling and share a `Budget` across production and checking.

`roots.refine` validates the retained conversion, tightens the original
algebraic coefficients under the new options, and builds a fresh cover. Higher
precision can reduce conversion uncertainty, but repeated roots or search
boundaries can still leave a numerical cover unresolved. Exact root solving
remains available when you need exact algebraic values instead of regions.

## Consume exact complex roots incrementally

`dense.complex_roots_to` and `sparse.complex_roots_to` send one
`roots.ExactComplexRootEvent` per distinct root. Each event contains `value`
and `multiplicity`. The returned receipt records the event count and total
multiplicity. Feed the stream to `complex_roots_checker`, then verify the
receipt with `finish`.

<!-- example: docs/examples/array_workflow_root_stream.mojo -->

The two roots are $i$ and $-i$, each with multiplicity two. The checker validates
their witnesses and reconstructs the original polynomial, rejecting repeated
events, incorrect multiplicities, and incomplete streams. A custom consumer
implements `roots.ExactComplexRootSink.write(event, budget)`.

Events follow deterministic factor and coordinate traversal order, which need
not be sorted. The stream discards completed output witnesses; the solver still
retains elimination and coordinate-solving state, and the checker retains the
input polynomial and a growing product. A consumer error or exhausted budget
can leave partial output, which is not a complete root result. Use the retained
`complex_roots` result when you need indexing, sorted output, or a complete
checked archive.
