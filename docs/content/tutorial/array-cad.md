# Solve inequalities with CAD

Sparse coefficient arrays can also describe real sign constraints. Use CAD to
find real points satisfying inequalities or to describe the cells on which
your polynomials have constant signs. For equations and ideal membership, see
[Gröbner bases](array-systems.md).

## State real sign constraints directly

For $x^2+y^2-1$, negative values describe the inside of the unit circle.
The following example finds an interior point, checks a contradictory request,
and constructs the complete circle decomposition.

<!-- example: docs/examples/array_solver_cad.mojo -->

`cad.negative`, `zero`, `positive`, `nonpositive`, `nonnegative`, and `nonzero`
each compare a polynomial to zero. To compare two polynomials, form their
difference with `sparse.subtract` first. Constraints accept integer and rational
coefficients and promote to rational storage, so differently typed constraints
can share the same collection.

`cad.solve` treats the collection as a conjunction. Its result reports
`satisfiable()` and retains either witness evidence or a covering proof of
failure. `cad.sample` returns exact real coordinates in polynomial axis order,
or `None` when no point satisfies the conjunction. Work exhaustion raises;
it is not an unsatisfiable result.

`cad.decompose` constructs the full sign-invariant decomposition. `cad.cells`
returns terminal cells with their dimensions, exact sample points, and input
signs. A sample identifies a cell; the certificate establishes the signs
throughout that cell. `cad.check` verifies either kind of result independently.
A selective solution is evidence for its conjunction, not a full decomposition.

## Control variable order and memory

CAD's `axes` parameter orders variables **from base to top**. This differs from
the highest-priority-first interpretation used for Gröbner bases. Both require
a permutation of every axis, including unused axes. Changing the order can
change the amount of work and the shape of a decomposition. Returned points
always use the original polynomial axis order.

By default `cad.decompose` and `cad.cover` use Brown–McCallum projection
reductions with Lazard lifting, interval-first root isolation and independent
final verification. `lazard=True` uses
plain Lazard projection with the same compact backend. `lazard=False` requests
the checked legacy McCallum backend, which can raise `UNSUPPORTED` on
nullification outside its supported cases. Fixed-coordinate and algebraic
coefficient prefixes also use checked legacy evidence. `cad.solve`,
`cad.decompose_to` and the other CAD operations default to Lazard projection:
a stream cannot fall back once McCallum refuses, because the sink has already
accepted its events.

For large outputs, use `cad.decompose_to` and `cad.checker`, as shown in the
[streaming example](array-roots.md#process-witnesses-as-they-arrive). They consume
projection and local cell witnesses incrementally. They still retain the
projection basis and active traversal state; streaming does not make the
underlying CAD computation constant-memory.

See [sparse arrays](../reference/array-sparse.md) and the
[CAD array reference](../reference/array-cad.md) for the complete signatures.
