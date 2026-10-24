# Example gallery

Every program listed here is checked against its expected output. The pages
include its source, a download link, and a run command. Begin with the first
tutorial, then use the later examples to see how operations fit together.

## Tutorials and guides

| Program | Walkthrough |
|---|---|
| `tutorial_first.mojo`, `tutorial_steps.mojo` | [Your first polynomial](../tutorial/first-polynomial.md) |
| `array_api.mojo`, `array_operations.mojo`, `array_surface_scalars.mojo`, `array_surface_format.mojo` | [Polynomial arrays](../tutorial/arrays.md) |
| `array_access_inspection.mojo`, `array_exact_algebra.mojo`, `array_factor_fields.mojo`, `array_exact_samples.mojo` | [Exact algebra](../tutorial/array-exact.md) |
| `array_parity_checked.mojo`, `array_exact_choices.mojo` | [Check exact algebra results](../tutorial/array-checked.md) |
| `array_solver_roots.mojo`, `array_access_sturm.mojo`, `array_complex_roots.mojo`, `array_solver_streams.mojo` | [Roots](../tutorial/array-roots.md) |
| `array_algebraic_points.mojo`, `array_root_diagnostics.mojo`, `array_extension_algebraic_roots.mojo` | [Evaluate at exact roots and inspect candidates](../tutorial/array-roots.md#evaluate-and-specialize-at-exact-roots) |
| `array_workflow_specialization.mojo`, `array_workflow_covers.mojo`, `array_workflow_root_stream.mojo` | [Algebraic coefficients in root searches](../tutorial/array-roots.md#diagnose-an-algebraic-polynomial-numerically) |
| `array_algebraic_coefficients.mojo`, `array_algebraic_factor.mojo`, `array_algebraic_multivariate_factor.mojo` | [Algebraic coefficients](../tutorial/algebraic-coefficients.md) |
| `array_solver_ideals.mojo` | [Gröbner bases and membership](../tutorial/array-systems.md) |
| `array_solver_cad.mojo` | [Inequalities and CAD](../tutorial/array-cad.md) |
| `array_parity_cad.mojo`, `array_fixed_cad.mojo`, `array_algebraic_cad.mojo` | [Specialized CAD workflows](../tutorial/array-cad-workflows.md) |
| `array_toolkit_series.mojo` | [Formal series](../tutorial/array-series.md) |
| `array_toolkit_basis.mojo` | [Alternative bases](../tutorial/array-bases.md) |
| `array_toolkit_rational.mojo`, `array_algebraic_rational.mojo`, `array_fraction_enclosures.mojo`, `array_workflow_fraction.mojo` | [Rational functions](../tutorial/array-rational.md) |
| `array_numeric.mojo`, `array_access_dtype.mojo`, `array_algebraic_numeric.mojo`, `array_extension_numerical.mojo`, `array_extension_representations.mojo` | [Numerical coefficients](../tutorial/array-numeric.md) |
| `array_extension_certified.mojo`, `array_parity_stream.mojo`, `array_algebraic_stream.mojo` | [Certified evaluation](../tutorial/array-certified.md) |
| `array_batching.mojo`, `array_batch_reductions.mojo`, `array_workflow_accumulate.mojo`, `array_batch_interop.mojo` | [Batch calculations](../tutorial/array-batching.md) |
| `array_value_serialization.mojo`, `array_checked_serialization.mojo`, `array_serialization.mojo`, `array_algebraic_serialization.mojo`, `array_surface_kernel.mojo` | [Save and restore values](../tutorial/array-serialization.md) |
| `tutorial_errors.mojo` | [Typed errors and budget exhaustion](../guides/errors.md) |

## Core API examples

These programs use the core `polynomials` package directly. They are for
developers extending the library; applications use the programs above.

| Program | Where to find it |
|---|---|
| `array_exact_core.mojo`, `array_access_interop.mojo`, `array_surface_interop.mojo` | [Exchange arrays with core values](../reference/core-api.md) |
| `arithmetic.mojo` | [Arithmetic and the kernel codec](../reference/p1.md#example) |
| `gcd.mojo`, `tutorial_division.mojo` | [Division, gcds, and cancellation certificates](../reference/p2.md#example) |
| `algebra.mojo`, `tutorial_factor.mojo` | [Factorization, bases, and roots](../reference/p3.md) |
| `tutorial_groebner.mojo`, `complex_groebner.mojo` | [Gröbner bases and membership](../reference/p4.md#example) |
| `tutorial_roots.mojo` | [Rational and irrational real roots](../reference/p5.md#example) |
| `toolkit.mojo`, `tutorial_calculus.mojo`, `tutorial_interpolation.mojo`, `tutorial_rational.mojo`, `tutorial_bases.mojo`, `tutorial_series.mojo` | [Polynomial toolkit](../reference/toolkit.md) |
| `streaming.mojo` | [Incremental text and proof checking](../reference/toolkit.md#incremental-text-and-cad-transport) |
| `numerical.mojo`, `tutorial_numerical.mojo`, `tutorial_refinement.mojo` | [Numerical certificates and refinement](../spec/numerical.md#example) |
| `algebraic_numbers.mojo` | [Exact real algebraic arithmetic](../reference/algebraic.md) |
| `cad.mojo` | [CAD and selective coverings](../reference/cad.md) |

To run a downloaded example, [install the package](../start/installation.md)
and run it with `pixi run mojo run`, or run it with `mojo -I build` inside a
checkout's Pixi environment.
