# CAD and sign constraints

Use `poly.cad` to express real polynomial sign conditions, find witnesses, and decompose space from sparse coefficient arrays over integers, rationals, or AA.

Every computational operation has a shared `Budget` overload and a convenience
overload with `max_steps=`. Borrowed evidence is immutable. Decomposition and
full-covering evidence now uses `CADFinalCertificate`; copies share immutable
storage. Call its metadata methods instead of reading the old `.projection`,
`.cells` or `.stacks` fields. `cad.cells` explicitly materializes absolute
geometry, and `cad.sample` materializes only the requested satisfying path.
Both costs are charged to their caller. For AA, decomposition and solving
follow only the supplied coefficient values. The borrowed tree retains rational
projection and local stack evidence, with one selected section per coefficient
axis. Use `cad.check` or `cad.checker` to verify that restricted claim and its
binding to the original coefficients; ordinary full-space core checkers have a
different coverage requirement. The
[workflow tutorial](../tutorial/array-cad-workflows.md) explains the raw stream
format and auxiliary axes.

Decomposition and full covering default to Brown–McCallum projection reductions
with Lazard lifting and the verified compact backend. `lazard=True` requests
plain Lazard on that backend;
`lazard=False` and restricted coordinate prefixes use checked legacy evidence.
`solve` and the streaming API retain their existing implementations.

Pass `fixed=[(axis, value), ...]` to restrict decomposition or solving to exact
coordinate values. `decompose`, `solve`, `cover`, `decompose_to`, and `checker`
accept homogeneous iterables or mixed tuples of axis/value pairs. Fixed axes
lead the traversal; samples retain every original axis and signs omit synthetic
binding polynomials. `cad.check` verifies coverage on that slice, including
its exact root choices. Use the same `fixed=` and `axes=` for a stream producer
and checker. Checked archives retain these bindings.

<!-- api: src/poly/cad.mojo Constraint -->

<!-- api: src/poly/cad.mojo negative -->

<!-- api: src/poly/cad.mojo zero -->

<!-- api: src/poly/cad.mojo positive -->

<!-- api: src/poly/cad.mojo nonpositive -->

<!-- api: src/poly/cad.mojo nonnegative -->

<!-- api: src/poly/cad.mojo nonzero -->

<!-- api: src/poly/cad.mojo Decomposition -->

<!-- api: src/poly/cad.mojo Solution -->

<!-- api: src/poly/cad.mojo decompose -->

<!-- api: src/poly/cad.mojo solve -->

<!-- api: src/poly/cad.mojo check -->

<!-- api: src/poly/cad.mojo sample -->

<!-- api: src/poly/cad.mojo Cell -->

<!-- api: src/poly/cad.mojo cells -->

<!-- api: src/poly/cad.mojo fixed_coordinates -->

<!-- api: src/poly/cad.mojo decompose_to -->

<!-- api: src/poly/cad.mojo checker -->

## Projection, selected cells, and lifting

The [specialized workflows tutorial](../tutorial/array-cad-workflows.md) explains axis order and how to consume lifting coefficients.

<!-- api: src/poly/cad.mojo Projection -->

<!-- api: src/poly/cad.mojo CellAt -->

<!-- api: src/poly/cad.mojo Lazard -->

<!-- api: src/poly/cad.mojo Covering -->

<!-- api: src/poly/cad.mojo project -->

<!-- api: src/poly/cad.mojo cell_at -->

<!-- api: src/poly/cad.mojo lazard_evaluate -->

<!-- api: src/poly/cad.mojo cover -->

<!-- api: src/poly/cad.mojo projection_polynomials -->

<!-- api: src/poly/cad.mojo signs -->

<!-- api: src/poly/cad.mojo vanishing_orders -->

<!-- api: src/poly/cad.mojo lifting_coefficients -->
