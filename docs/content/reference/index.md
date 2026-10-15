# API reference

Import `poly` for coefficient arrays and `polynomials` for the core algebra
and certificate APIs. This reference groups operations by task; declarations and docstrings are rendered directly from source.
Start with a [tutorial](../tutorial/first-polynomial.md) if you are learning
the workflow, or choose an operation below when you need its precise contract.

| Topic | Reference |
|---|---|
| Coefficient arrays, exact algebra, roots, Gröbner bases, CAD, and core conversion | [Array API](array-api.md) |
| Save complete values and shaped batches | [Value serialization](array-io.md) |
| Formal series, basis representations, and rational functions from coefficient arrays | [Series](array-series.md), [bases](array-basis.md), [rational functions](array-rational.md) |
| Domains, contexts, containers, orders, arithmetic, errors, and binary codec | [Representations and arithmetic](p1.md) |
| Division, gcds, resultants, fractions, and identity certificates | [Division and elimination](p2.md) |
| Squarefree, finite-field, integer, rational, and Gaussian factorization | [Univariate factorization](p3.md) |
| Bases, reduction, and ideal membership | [Gröbner bases](p4.md) |
| Sturm sequences, root counts, isolation, and refinement | [Exact real roots](p5.md) |
| Calculus, interpolation, bases, series, rational utilities, and text transport | [Polynomial toolkit](toolkit.md) |
| Point and ball arithmetic, enclosures, root covers, and incremental witnesses | [Numerical polynomials](../spec/numerical.md) |
| Exact algebraic arithmetic and specialization | [Real algebraic numbers](algebraic.md) |
| Sparse modular gcds and multivariate factors | [Multivariate factorization](advanced-factorization.md) |
| Projection, cells, complete decomposition, and selective coverings | [CAD](cad.md) |

## Read a signature in context

Exact operations are generic over a coefficient ring or field where supported.
The [domain guide](../guides/domains.md) explains the distinction. Numerical
operations take explicit arithmetic contexts. Most computational calls accept
a mutable `Budget` and raise `PolyError`; see [errors and work limits](../guides/errors.md).

Certificate checkers establish the relations described by their own contracts.
Read [certificate claims](../guides/certificates.md) before treating an accepted
product, identity, or enclosure as a stronger statement.

The exact public surface is recorded in [Contract P](contract-p.md). The
numerical declarations have a separate snapshot. Historical P1–P5 labels in
source paths and design records name implementation groups; they are not
steps a user needs to complete.

The array API also includes [numerical coefficient arrays](array-numeric.md) with explicit precision and [batch maps and broadcasting](array-batch.md).

The array interface also includes [certified evaluation results](array-evaluation.md).
