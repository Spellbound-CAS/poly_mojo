# Simplify and expand rational functions

Use `poly.rational` for a ratio of coefficient polynomials. Construction
cancels common factors and makes the denominator monic. Arithmetic,
differentiation, and substitution return normalized results, while partial
fractions and Laurent expansions retain their original function for checking.

## Start with a numerator and denominator

`rational.asarray([-1,0,1], [-1,1])` constructs $(x^2-1)/(x-1)$ and reduces it
to $x+1$. The following example also decomposes $1/(x(1-x))$ and expands it
around its pole at zero.

<!-- example: docs/examples/array_toolkit_rational.mojo -->

`rational.asarray(coefficients)` embeds a polynomial with denominator one.
Use `rational.from_terms(numerator, denominator)` for sparse input, including
multivariate exponent tuples. Both inputs must have the same power shape.
Integer coefficients promote to `QQ`. Gaussian rationals select `QQi`, real
algebraic coefficients select `AA`, and complex algebraic coefficients select
`QQbar`. Mixing `AA` with `QQi` also selects `QQbar`. Modular coefficients
select `GF(p)`; both inputs must use the same prime. Use
`poly.coefficients((...))` when a dense input mixes scalar types.

Evaluation uses the reduced function, so the first result evaluates to 2 at
$x=1$. The original expression still excluded that point. Keep any original
domain exclusions separately when simplifying a partially defined expression.
A remaining pole raises `DOMAIN`.

## Combine functions and substitute variables

Use `add`, `subtract`, `multiply`, `divide`, `negative`, and `reciprocal` on
prepared functions. Division by a zero function raises `DOMAIN`. Extract the
reduced polynomials with `rational.numerator` and `rational.denominator`; both
return sparse arrays with the original power shape.

`evaluate` takes a scalar for one variable or a coordinate tuple for several.
`derivative(f, axis=...)` differentiates along a zero-based array axis.
`substitute(f, replacements, axes=...)` accepts sparse replacement polynomials
and applies them simultaneously. Without `axes=`, supply one replacement per
axis. Repeated axes and a substituted denominator that becomes identically
zero raise `DOMAIN`.

`compose(outer, inner)` composes univariate rational functions. Partial
fractions, poles, Laurent expansions, and residues also require a single
occurring variable. Their `axis=` option selects that variable in a larger
power shape; unused axes remain part of the shape.

## Use algebraic poles and expansion centers

The field follows your data throughout the calculation. You can evaluate a QQ
function at an algebraic point, combine AA and QQbar functions, or expand a QQ
function around $\sqrt{2}$ without specifying a ring.

<!-- example: docs/examples/array_algebraic_rational.mojo -->

The example expands $1/(x^2-2)$ at $\sqrt{2}$. Its valuation is -1, and its
residue is $1/(2\sqrt{2})$. The checker verifies the exact coefficients by
multiplication. Arithmetic, composition, and polynomial substitution use the
same promotion rules as construction.

## Read the decomposition and its evidence

`partial_fractions(f)` returns `.polynomial` and `.terms`. Each term is
`(numerator, denominator_factor, exponent)`, using dense arrays. `rational.check`
checks complete cancellation in the source function, irreducibility of every
monic denominator factor, and exact recomposition. Terms must have nonzero
numerators of lower degree than their factor. Each `(factor, exponent)` pair
appears at most once; omitted zero terms are allowed.

The example also checks $1/(x^4+1)$ over QQ. Its denominator remains one
irreducible quartic. Changing the coefficient field can change the decomposition:
over AA it splits into two quadratics, and over QQbar into four linear factors.

`poles(f)` returns a denominator scalar and `(factor, multiplicity)` entries.
Each factor represents its complex zeros. Over `AA`, factors are linear or
irreducible quadratic; over `QQbar`, every factor is linear. Multiplicities
refer to the reduced denominator, so cancelled poles are absent.

`laurent(f, center, terms=n)` returns `.center`, `.valuation`, and `.coefficients`.
The first coefficient multiplies $(x-center)^{valuation}$; each following
coefficient increases the exponent by one. The example produces valuation -1
and four ones: $x^{-1}+1+x+x^2+O(x^3)$. An algebraic center promotes the computation to `AA` or `QQbar`.
`residue(f, center)` extracts the coefficient of exponent -1 directly.
`rational.check(expansion)` verifies the valuation and supplied coefficients
by a truncated multiplication identity. No analytic remainder bound is implied.

On a function itself, `rational.check(f)` checks the retained cancellation
identity, proves that no common factor remains, and verifies the monic
denominator. This also certifies the normalized zero function as $0/1$.
`.evidence()` borrows the cancellation identity; `.coprimality_evidence()`
borrows the core proof for the returned pair. Keep the complete function for
`rational.check`. These checks also work for multivariate fractions such as
$x/y$, whose numerator and denominator are coprime.

Irreducibility checking over QQ first tries prime reduction. When that cannot
settle the question, a finite integer-interpolation search checks all possible
proper divisors. This fallback can be expensive. QQi uses a shifted conjugate
norm and the same QQ checks; AA uses linear factors or quadratics with negative
discriminant, and QQbar uses linear factors. All verification work shares your
budget. `LIMIT` means checking did not finish, never that a factor is irreducible.

All computational calls accept `max_steps=` or a shared `Budget`. The first
factorization may also populate the core prime cache; later calls reuse it
and can consume less work. Results share
owned storage, and core conversion copies it so later core mutation cannot
change the prepared function. See the [rational API reference](../reference/array-rational.md)
for signatures and [normalization theory](../theory/rational-functions.md) for
the normalization guarantee.

## Exchange values with the core

`to_core(f)` returns a core `RationalFunction` for QQ. For QQi, GF(p), AA, and QQbar,
it returns `(numerator, denominator)` as public `SparsePoly` values, since the
core rational-function container is QQ-specific. Both forms own independent
storage. Import either form with `from_core[ndim=n](value)`; `variables=` maps
core variable IDs into array-axis order. Zero-axis constants are supported.

## Work over a finite field

Use `modint(value, prime)` coefficients to infer `GF(p)`. The fraction keeps
that modulus through arithmetic, differentiation, substitution, composition,
partial fractions, Laurent expansion, and serialization. Both operands and
any evaluation or expansion point must belong to the same field. Integers
inside a finite-field coefficient sequence should first be converted with
`modint`; incompatible moduli raise `DOMAIN`.

Partial fractions retain an irreducibility certificate for each denominator
factor. The checker uses Frobenius powers and polynomial GCDs in the stored
field. `poles` returns denominator factors with multiplicities; over `GF(p)`,
those factors describe roots in an algebraic extension, not complex numbers.

## Choose numerical precision explicitly

Pass `context=` when building a numerical fraction, either from coefficients
or by converting an exact fraction. Arithmetic, derivatives, integer powers,
composition, and substitution return new values under the context you supply.
No approximate GCD is used to cancel factors. `to_core` and `from_core` copy
numerator and denominator as a pair of numerical sparse polynomials.

For a certified value, call `rational.enclose(f, point, context=...)` and
`rational.check(result)`. The checker verifies both polynomial evaluations at
the same point and bounds their quotient with exact rational interval
arithmetic. It rejects a denominator enclosure containing zero. That can
mean a pole, uncertainty in ball coefficients, or insufficient precision;
raising precision helps only when the underlying inputs separate from zero.

<!-- example: docs/examples/array_fraction_enclosures.mojo -->

Use `result.value()` to copy the certified ball. Save the whole result with
`poly.io` when you also need the evidence. Rounded `evaluate` remains useful
for ordinary approximate calculations; `enclose` is the operation that retains
a checkable enclosure claim.

## Request accuracy for the quotient

Use `rational.adaptive_evaluate(f, point, options=...)` when you know the
required accuracy but not the working precision. Exact inputs use `QQ`, `QQi`,
`AA`, or `QQbar`; finite fields have no ball-enclosure interpretation. It
evaluates the numerator and denominator at the same point, increases precision as needed, and tests
the width of the final quotient. Exact algebraic coefficients and points are
refined on each attempt; supplied numerical uncertainty is preserved.

<!-- example: docs/examples/array_workflow_fraction.mojo -->

Check `has_value()` before extracting a finite ball, and `accurate()` before
relying on the requested width. A `PRECISION_CEILING` result may still contain
a useful, wider enclosure. If its denominator could not be separated from
zero, `has_value()` is false and `value()` raises `DOMAIN`. That outcome alone
does not distinguish a pole from uncertainty or insufficient precision.

`rational.check` verifies both polynomial traces, their common point, the
quotient bounds when present, and the accuracy outcome. An honestly unresolved
result can pass that check; it is not a successful accuracy claim. `poly.io`
preserves this evidence, including outcomes without a finite quotient. All
attempts share the caller's work budget. Scalar points and coordinate tuples
work for both exact and numerical fractions.

## Map variables at the core boundary

`rational.to_core(f, variables=ids)` assigns one core variable ID to each
array axis, including unused axes. Restore with
`rational.from_core[ndim](core, variables=ids)`. IDs may be signed and need
not increase; their positions define array-axis order. Duplicate IDs or
unmapped variables raise `DOMAIN`.

The same mapping works for exact and numerical fractions. Numerical copies
preserve stored formats and coefficient bounds without a new context. Imports
still validate both polynomials and require a denominator family that excludes
the identically zero polynomial. Exported storage is independent of `f`.
