# Rational functions from arrays

Build normalized rational functions over QQ, QQi, GF(p), AA, or QQbar from dense coefficients or sparse terms, then calculate, decompose, and verify their normalization and partial fractions.

`check(f)` verifies cancellation, coprimality, and the monic denominator.
`check(parts)` also verifies irreducible denominator factors, proper nonzero
terms without duplicate factor/exponent pairs, and recomposition. Laurent
checks verify source normalization and a formal truncated-product identity.

Every computational call accepts a shared `Budget` or `max_steps=`. See the
[tutorial](../tutorial/array-rational.md) for examples, domain restrictions,
and the meaning of retained metadata.

For numerical fractions, pass `context=` to the constructor and each arithmetic
operation. `NumericFunction` stores a numerator and denominator without
approximate cancellation. `enclose` produces independently checked quotient
bounds; `check` verifies the polynomial evaluations, their common point, and
the final division. A denominator enclosure containing zero raises `DOMAIN`.
Numerical partial fractions and Laurent expansions are not exact operations
and are not provided by these overloads.

Save fractions, partial fractions, Laurent expansions, and checked enclosures
with [poly.io](array-io.md).

<!-- api: src/poly/rational.mojo Function -->

<!-- api: src/poly/rational.mojo asarray -->

<!-- api: src/poly/rational.mojo from_terms -->

<!-- api: src/poly/rational.mojo numerator -->

<!-- api: src/poly/rational.mojo denominator -->

<!-- api: src/poly/rational.mojo add -->

<!-- api: src/poly/rational.mojo subtract -->

<!-- api: src/poly/rational.mojo multiply -->

<!-- api: src/poly/rational.mojo divide -->

<!-- api: src/poly/rational.mojo negative -->

<!-- api: src/poly/rational.mojo reciprocal -->

<!-- api: src/poly/rational.mojo evaluate -->

<!-- api: src/poly/rational.mojo derivative -->

<!-- api: src/poly/rational.mojo compose -->

<!-- api: src/poly/rational.mojo substitute -->

<!-- api: src/poly/rational.mojo check -->

<!-- api: src/poly/rational.mojo PartialFractions -->

<!-- api: src/poly/rational.mojo partial_fractions -->

<!-- api: src/poly/rational.mojo poles -->

<!-- api: src/poly/rational.mojo Expansion -->

<!-- api: src/poly/rational.mojo laurent -->

<!-- api: src/poly/rational.mojo residue -->

<!-- api: src/poly/rational.mojo to_core -->

<!-- api: src/poly/rational.mojo from_core -->

<!-- api: src/poly/_numeric_fraction.mojo NumericFunction -->

<!-- api: src/poly/_numeric_fraction.mojo Enclosure -->

<!-- api: src/poly/rational.mojo power -->

<!-- api: src/poly/rational.mojo enclose -->

`adaptive_evaluate` accepts `NumericOptions` and increases precision for the
whole quotient. Inspect `has_value()`, `accurate()`, and `outcome()` separately:
a valid certificate can describe a precision ceiling without a finite quotient.

<!-- api: src/poly/rational.mojo adaptive_evaluate -->

<!-- api: src/poly/_numeric_fraction.mojo AdaptiveEnclosure -->
