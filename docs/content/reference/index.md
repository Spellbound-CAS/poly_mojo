# API reference

Import `poly` for polynomial calculations. This reference groups operations by
task; declarations and docstrings are rendered directly from source.
Start with a [tutorial](../tutorial/first-polynomial.md) if you are learning
the workflow, or choose an operation below when you need its precise contract.

| Topic | Reference |
|---|---|
| Construction, shared conventions, budgets, and errors | [Array API](array-api.md) |
| Arithmetic, calculus, division, gcds, factorization, interpolation, and elimination | [Dense arrays](array-dense.md), [sparse arrays](array-sparse.md) |
| Retained evidence for division, gcds, resultants, factorization, and interpolation | [Checked exact results](array-exact-results.md) |
| Exact real and complex roots, root counts, and numerical root covers | [Root results and refinement](array-roots.md) |
| Real and complex algebraic numbers | [Exact algebraic values](array-algebraic.md) |
| Gröbner bases, reduction, and ideal membership | [Sparse arrays](array-sparse.md#roots-and-polynomial-systems) |
| Sign constraints, samples, and cylindrical decompositions | [CAD and sign constraints](array-cad.md) |
| Formal series, basis representations, and rational functions | [Series](array-series.md), [bases](array-basis.md), [rational functions](array-rational.md) |
| Floating-point and ball coefficients | [Numerical coefficients](array-numeric.md) |
| Certified evaluation and accuracy targets | [Certified evaluation](array-evaluation.md) |
| Maps, broadcasting, and reductions over batches | [Batch maps and broadcasting](array-batch.md) |
| Save complete values and shaped batches | [Value serialization](array-io.md) |

## Read a signature in context

Each computational function has an overload that takes a mutable `Budget` and
one that takes `max_steps=`; both raise `PolyError`. See
[errors and work limits](../guides/errors.md). Numerical operations take an
explicit `context=`. The [domain guide](../guides/domains.md) lists the
coefficient domains each kind of operation accepts.

Certificate checkers establish the relations described by their own contracts.
Read [certificate claims](../guides/certificates.md) before treating an accepted
product, identity, or enclosure as a stronger statement.

Developers extending the library can use the [core API](core-api.md), which
documents the containers and algorithms underneath `poly`.
