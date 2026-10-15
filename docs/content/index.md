# Polynomial algebra in Mojo

<div class="hero" markdown="1">

Poly Mojo brings exact algebra and arbitrary-precision numerical computation
together in one pure Mojo library. It includes polynomial arithmetic, calculus,
interpolation, factorization, Gröbner bases, real and complex roots, exact
algebraic numbers, and cylindrical algebraic decomposition. Build polynomials
from coefficient sequences or sparse power tuples, and use certificates
and interval enclosures to check the claims your calculation depends on.

The library is pure Mojo and builds on apn_mojo, with no foreign algebra
library. Start with an exact calculation, then add numerical precision or
certificate checking where your application needs it.

<div class="hero-actions" markdown="1">
[Get started](start/installation.md){ .button .primary }
[First polynomial](tutorial/arrays.md){ .button }
[API reference](reference/index.md){ .button }
</div>

</div>

## Start with a calculation

The [array tutorial](tutorial/arrays.md) constructs polynomials from coefficients,
evaluates them, and differentiates sparse terms. It introduces inferred domains
and optional work limits through a complete, checked program.

| What you want to do | Where to go |
|---|---|
| Choose exact, rounded, or interval coefficients | [Coefficient domains](guides/domains.md) |
| Choose storage and understand coefficient order | [Dense and sparse polynomials](tutorial/representations.md) |
| Divide and find common factors | [Division and gcds](tutorial/division.md) |
| Simplify fractions or find residues | [Rational functions](tutorial/array-rational.md) |
| Factor a polynomial | [Factorization](tutorial/factorization.md) |
| Differentiate, integrate, and substitute | [Calculus and transformations](tutorial/calculus.md) |
| Reconstruct a polynomial from samples | [Interpolation](tutorial/interpolation.md) |
| Use Newton, Chebyshev, or Bernstein coefficients | [Alternative bases](tutorial/array-bases.md) |
| Compute a truncated expansion | [Formal power series](tutorial/array-series.md) |
| Isolate real roots exactly | [Real roots](tutorial/real-roots.md) |
| Evaluate with arbitrary precision and interval bounds | [Numerical coefficient arrays](tutorial/array-numeric.md) |
| Apply a function across polynomial batches | [Maps and broadcasting](tutorial/array-batching.md) |
| Save and restore coefficients | [Canonical serialization](tutorial/array-serialization.md) |
| Find and refine validated real or complex roots | [Numerical roots](tutorial/numerical-roots.md) |
| Reduce polynomial systems | [Gröbner bases and membership](tutorial/groebner.md) |
| Compute with algebraic values or decompose real space | [Algebraic numbers](tutorial/algebraic-numbers.md), [CAD](tutorial/cad.md) |

## Understand the result

An exact polynomial represents a single mathematical polynomial. A numerical
point polynomial uses rounded coefficients; a ball polynomial represents a
family of possible coefficients. These choices affect which algorithms apply
and what an answer means. The [domain guide](guides/domains.md) helps you choose.

A certificate contains evidence for a particular claim. For example, a
factorization product checker verifies that the proposed factors multiply
back to the input; it does not prove their irreducibility. The
[certificate guide](guides/certificates.md) explains these boundaries and the
separate guarantees of numerical enclosures and root covers.

Computational operations take a shared `Budget`. Numerical searches also
accept accuracy and subdivision settings. See [budgets and errors](guides/errors.md)
for handling work exhaustion and inconclusive numerical outcomes.

## Learn how it works

The [architecture overview](architecture/overview.md) explains the package's
layers and data flow. The [mathematical chapters](theory/index.md) give
algorithms, proofs, and certificate relations. Maintainers can start with the
[maintenance guide](contributing/handoff.md), then consult the contracts and
recorded measurements for the component they are changing.

See the [example gallery](examples/index.md) for complete programs and
[support and limitations](guides/status.md) for the current scope.

Continue with [numerical arrays](tutorial/array-numeric.md), [batch maps and broadcasting](tutorial/array-batching.md), or [saving coefficient data](tutorial/array-serialization.md).
