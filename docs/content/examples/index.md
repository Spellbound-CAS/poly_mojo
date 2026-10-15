# Example gallery

Every program listed here is checked against its expected output. The pages
include its source, a download link, and a run command. Begin with the focused
tutorials, then use the broader examples to see how operations fit together.

## Focused tutorials

| Program | Walkthrough |
|---|---|
| `array_api.mojo`, `array_operations.mojo` | [Polynomial arrays](../tutorial/arrays.md) |
| `array_exact_algebra.mojo`, `array_exact_samples.mojo`, `array_exact_core.mojo` | [Exact algebra and core interoperability](../tutorial/array-exact.md) |
| `array_toolkit_series.mojo` | [Formal series from arrays](../tutorial/array-series.md) |
| `array_toolkit_basis.mojo` | [Basis coefficients and conversion](../tutorial/array-bases.md) |
| `array_toolkit_rational.mojo` | [Rational functions from arrays](../tutorial/array-rational.md) |
| `tutorial_named.mojo` | [Named expression input](../tutorial/named-polynomials.md) |
| `tutorial_first.mojo` | [Construct, evaluate, and differentiate](../tutorial/first-polynomial.md) |
| `tutorial_division.mojo` | [Division with a checked remainder](../tutorial/division.md) |
| `tutorial_calculus.mojo` | [Integration and substitution](../tutorial/calculus.md) |
| `tutorial_interpolation.mojo` | [Recover a polynomial from samples](../tutorial/interpolation.md) |
| `tutorial_bases.mojo` | [Bernstein and Chebyshev representations](../tutorial/bases.md) |
| `tutorial_series.mojo` | [A reciprocal as a formal series](../tutorial/series.md) |
| `tutorial_rational.mojo` | [Partial fractions and a residue](../tutorial/rational-functions.md) |
| `tutorial_factor.mojo` | [Factorization with multiplicities](../tutorial/factorization.md) |
| `tutorial_roots.mojo` | [Rational and irrational real roots](../tutorial/real-roots.md) |
| `tutorial_groebner.mojo` | [A basis and ideal membership](../tutorial/groebner.md) |
| `tutorial_numerical.mojo` | [Rounded values and rigorous bounds](../tutorial/numerical.md) |
| `tutorial_refinement.mojo` | [Refine a checked numerical root cover](../tutorial/numerical-roots.md) |
| `tutorial_errors.mojo` | [Typed errors and budget exhaustion](../guides/errors.md) |

## Further examples

| Program | Where to find it |
|---|---|
| `arithmetic.mojo` | [Arithmetic and the kernel codec](../reference/p1.md#example) |
| `gcd.mojo` | [Gcd and cancellation certificates](../reference/p2.md#example) |
| `algebra.mojo` | [Factorization, bases, and roots](../reference/p3.md) |
| `complex_groebner.mojo` | [Exact complex Gröbner bases](../reference/p4.md#exact-complex-coefficients) |
| `algebraic_numbers.mojo` | [Exact real algebraic arithmetic](../tutorial/algebraic-numbers.md) |
| `cad.mojo` | [CAD and selective coverings](../tutorial/cad.md) |
| `numerical.mojo` | [Real and complex numerical certificates](../spec/numerical.md#example) |
| `toolkit.mojo` | [Combined toolkit workflows](../reference/toolkit.md) |
| `streaming.mojo` | [Incremental text and proof checking](../guides/interchange.md) |

To run a downloaded example, [install the package](../start/installation.md)
and run it with `pixi run mojo run`, or run it with `mojo -I build` inside a
checkout's Pixi environment. Contributors can run the complete example
matrix with the commands in [maintaining the docs](../contributing/documentation.md).

Find [roots from coefficient arrays](../tutorial/array-roots.md), reuse a
[Gröbner basis or solve sign constraints](../tutorial/array-systems.md), and
check root and CAD witnesses incrementally through the array interface.

See [numerical arrays](../tutorial/array-numeric.md), [batch calculations](../tutorial/array-batching.md), [serialization](../tutorial/array-serialization.md), and [algebraic points and root diagnostics](../tutorial/array-roots.md#evaluate-and-specialize-at-exact-roots) for the corresponding runnable programs.
