# Architecture overview

`polynomials` separates the representation of a polynomial, the algorithm that
computes a result, and the evidence used to check that result. This lets the
library use fast search and reconstruction methods while keeping certificate
checking smaller and focused on explicit mathematical claims.

The package is standalone pure Mojo. Its runtime imports are `std` and
apn_mojo, which supplies the arbitrary-precision scalar types. It does not
require a C algebra library or a Python interpreter to run a compiled
consumer.

`poly` adds a coefficient-first interface with inferred exact domains, explicit
numerical contexts, immutable dense and sparse arrays, and shaped batch maps. It depends on public `polynomials` exports
and APN scalar types. See [representations](representations.md) for its storage
and [budgets and errors](../guides/errors.md) for work accounting. The module table below describes the core package.

## From coefficients to certified results

```text
apn_mojo scalars
    |
    +-- exact coefficient rings and fields
    |       +-- DensePoly / SparsePoly
    |       +-- exact algebra, real algebraic values, CAD
    |
    +-- numerical domain adapters + explicit contexts
            +-- NumericDensePoly / NumericSparsePoly
            +-- rounded operations, enclosures, root covers

algorithms --> results and certificates --> independent checkers
                         |
                         +--> materialized values or incremental sinks
```

Exact and numerical containers share monomials, ordering, work budgets, and
applicable traversal helpers. They have different scalar contracts: rounded
arithmetic and interval families do not satisfy the exact coefficient-field
laws. Numerical modules may use exact helpers, but exact algorithms do not
depend on numerical modules.

## Where to look in the source

| Responsibility | Main modules |
|---|---|
| Public imports and errors | `__init__.mojo`, `errors.mojo` |
| Exact scalar adapters | `coefficients.mojo`, `exact_complex.mojo` |
| Representation and construction | `monomial.mojo`, `dense.mojo`, `sparse.mojo`, `builders.mojo` |
| Division and elimination | `division.mojo`, `gcd.mojo`, `resultant.mojo`, `fractions.mojo` |
| Factorization and ideals | `squarefree.mojo`, `finite_field.mojo`, `factor.mojo`, `multivariate_factor.mojo`, `gaussian_factor.mojo`, `field_multivariate_factor.mojo`, `sparse_squarefree.mojo`, `sparse_gcd.mojo`, `groebner.mojo` |
| Exact roots and real geometry | `roots.mojo`, `algebraic*.mojo`, `cad*.mojo`, `lazard.mojo`, `cylindrical_covering.mojo` |
| General polynomial operations | `calculus.mojo`, `interpolation.mojo`, `bases.mojo`, `series.mojo`, `rational_tools.mojo` |
| Numerical operations and evidence | `numeric_*.mojo` |
| Interchange and streams | `interchange.mojo`, `numeric_interchange.mojo`, `text_stream.mojo`, `numeric_stream.mojo`, `root_stream*.mojo`, `cad_stream.mojo` |
| Exact identity evidence | `certificates.mojo` and the specialized certificate modules |
| Kernel integration | `codec.mojo` |

Modules prefixed with `_` hold private algorithms and shared implementation
helpers. Consumer code imports the package root; private module paths are not
part of the supported interface.

## Producers and checkers

A producer finds factors, a basis, a root region, or a decomposition. Its
checker verifies the relevant evidence against the original inputs. The
checker does not repeat the search. For example, membership checking verifies
$f=\sum_i h_i g_i$ without computing a Gröbner basis, and numerical enclosure
checking verifies rational bounds without rerunning floating-point evaluation.

The claim determines the proof format. A factor product is sufficient to
prove recomposition but not irreducibility. A CAD projection needs determinant
evidence stronger than a resultant ideal-membership identity. The
[certificate guide](../guides/certificates.md) summarizes these boundaries.

## Determinism, budgets, and ownership

Canonical forms and prescribed ordering make results reproducible. Modular
algorithms use fixed prime sequences; randomized finite-field splitting uses
explicit deterministic seeds. Generator order remains significant where it
indexes certificate columns. Cache state may affect work charges, but cannot
change the mathematical result.

A caller-owned `Budget` travels through each computation, including nested
helpers. Work is charged before it is performed; representational growth is
checked before allocation. Borrowed input values remain unchanged on failure.
Internal code must also respect Mojo ownership: a view into a list cannot
survive a call that may grow that list.

## Integration and further reading

The binary codec bridges exact rational polynomials to a proof kernel's
`FORMAL_POLYNOMIAL` byte format. Library certificate acceptance and numerical
enclosure checking do not imply acceptance by that kernel.

Read [representations](representations.md) for storage and shared kernels,
[exact algorithms](algorithms.md) for computation and proof choices,
[numerical certification](numerical.md) for rounding and family semantics, and
[streaming](streaming.md) for memory lifetimes.
