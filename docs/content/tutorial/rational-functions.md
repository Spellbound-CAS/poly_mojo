# Rational functions

For coefficient iterables and functional calls, start with the [array interface](array-rational.md).

Use `RationalFunction` for a ratio of polynomials over `QQ`. It normalizes
common factors and denominator scaling, and works with evaluation,
differentiation, substitution, partial fractions, and local Laurent expansions.

## Construct a normalized fraction

`RationalFunction.normalize(numerator, denominator, work)` returns the
normalized function and a cancellation certificate. The denominator must be
nonzero. Construct the two polynomials separately; the named polynomial
parser accepts scalar division, not a general rational expression.

Normalization preserves the fraction wherever its original denominator is
nonzero. For example, $(x^2-1)/(x-1)$ simplifies to $x+1$, but the original
expression still excludes $x=1$. Keep that original domain condition when
using the result as an identity of functions.

The [gcd example](division.md#compute-and-check-a-gcd) shows normalization and
checks the cancellation evidence.

## Decompose a univariate fraction

For $f(x)=x^2/(x-1)$, polynomial division gives

$$
f(x)=x+1+\frac{1}{x-1}.
$$

Its value at $2$ is $4$, and its residue at the pole $1$ is $1$.

<!-- example: docs/examples/tutorial_rational.mojo -->

`partial_fractions` returns a polynomial part and proper fraction terms.
Each term stores a numerator, a base denominator factor, and a positive
exponent. `check_partial_fractions` checks properness and recomposition;
it does not independently certify the factors' irreducibility. Use
`poly.rational.check` on an [array partial-fraction result](array-rational.md#read-the-decomposition-and-its-evidence)
when you also need that guarantee.

## Differentiate, substitute, and expand

`rational_derivative` differentiates with respect to a variable.
`rational_substitute` performs simultaneous polynomial substitution, while
`rational_compose` composes two univariate rational functions. Evaluation
at a pole raises a domain error.

`rational_poles` represents complex poles by exact rational denominator
factors and multiplicities. It does not return a list of numerical complex
root approximations.

`laurent_expansion` computes coefficients around a rational center. The
result records the center and the starting exponent, so negative exponents
retain the pole's principal part. `rational_residue` extracts the coefficient
of $(x-a)^{-1}$. A finite Laurent expansion carries its algebraic certificate;
it does not itself supply an analytic tail-error estimate.

See the [rational-function reference](../reference/toolkit.md#rational-functions)
and [normalization theory](../theory/rational-functions.md) for the precise
identities and supported boundaries.
