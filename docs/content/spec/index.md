# Implementation contracts

These pages define representation invariants, algorithm requirements, and
verification obligations for the implemented library. They are intended for
maintainers investigating behavior or reviewing a change. Use the
[API reference](../reference/index.md) for current source declarations and
[tutorials](../tutorial/first-polynomial.md) to learn a workflow.

The [theory](../theory/index.md) supplies definitions and proofs. The
[architecture](../architecture/overview.md) explains the current module
boundaries and data flow. Public changes must preserve the compatibility
requirements in [Contract P](../reference/contract-p.md) and the separate
numerical snapshot, or include an agreed migration.

## How to read it

**MUST**, **MUST NOT**, **SHOULD**, and **MAY** express requirements with their
RFC 2119 meanings. Identifiers remain stable so source comments and test
names can refer to a precise requirement.

| Identifier | Meaning | Reference |
|---|---|---|
| G and PG | Shared and polynomial-specific engineering requirements | [Ground rules](ground-rules.md) |
| P1–P5 | Exact-core contract groups | Topic pages below |
| R, V, G, E, F, Z, B, T, A, N | Mathematical definitions and results | [Theory index](../theory/index.md) |

Text-fenced signatures on the contract pages describe the design. The
generated API reference and declaration snapshots record the actual Mojo
spelling, including generic constraints and typed errors.

## Repository map

| Contract | Subject |
|---|---|
| [Representations and arithmetic](p1-arithmetic.md) | Domains, containers, orders, budgets, and codec |
| [Division and elimination](p2-gcd.md) | Division, gcd, resultants, and fractions |
| [Factorization](p3-factorization.md) | Squarefree, finite-field, integer, rational, and Gaussian factors |
| [Gröbner bases](p4-groebner.md) | Reduction and membership |
| [Exact real roots](p5-roots.md) | Sturm counts, isolation, and refinement |
| [Numerical operations](numerical.md) | Rounded arithmetic, enclosures, and root covers |
| [Validation](acceptance.md) | Validation obligations and cross-checks |

Later exact extensions are specified by the [algebraic](../reference/algebraic.md),
[multivariate factorization](../reference/advanced-factorization.md),
[CAD](../reference/cad.md), and [toolkit](../reference/toolkit.md) references.
The current source map is in the [architecture overview](../architecture/overview.md).

## Conventions

- **Certificates are values.** Every function that returns a certified
  result returns a struct holding the result and its certificate; the
  matching `check_*` function recomputes the certificate's relation.
- **Budgets.** Every public function whose work depends on its input
  takes `mut budget: Budget` and charges it before the work it measures.
  This includes queries, negation, coefficient arithmetic, and the codec.
  Constant-time accessors and fixed-signature language protocols are the
  only exceptions; protocols used within algorithms are precharged by their
  caller.
- **Errors.** Public functions raise `PolyError`, never plain `Error`:
  apn_mojo's errors are caught at the call site and converted.
- **Rings.** Every polynomial carries the context of its coefficient ring
  next to its monomial order; constants are built from the context, and
  results inherit both from their inputs.
- **Pure Mojo.** No foreign code and no copyleft sources (ground rule
  PG7).
