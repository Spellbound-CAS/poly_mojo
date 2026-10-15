# Factor a polynomial

Factorization separates a polynomial into a scalar and powers of normalized
factors. The coefficient domain determines which factors are available:
$x^2+1$ has no linear rational factor, although it has complex roots.

For coefficient and term inputs through `poly`, see
[exact array algebra](array-exact.md). This walkthrough uses the core types
and their explicit contracts.

## Factor over the rationals

The polynomial $(x-1)^2(x+2)$ has two distinct rational factors, one with
multiplicity two.

<!-- example: docs/examples/tutorial_factor.mojo -->

`factor` accepts `DensePoly[QQ]`. The result contains a scalar `unit` and
`factors`, a list of `(polynomial, multiplicity)` pairs. The product of the
unit and these powers reconstructs the input.

Rational factors are monic. Integer factorization uses `factor_integer`;
its factors are primitive with positive leading coefficients, with the
remaining scalar and sign stored in the unit. The result order is
deterministic, which makes serialization and comparison reproducible.

## Check the claim you need

`check_factorization` independently checks the product equation. It does not
prove that each supplied factor is irreducible. A correct reconstruction
can still contain a factor that splits further.

`squarefree` groups factors by multiplicity. Its product checker likewise
checks reconstruction, without separately proving squarefreeness or pairwise
coprimality. The [certificate guide](../guides/certificates.md) explains how
to distinguish an algorithm's result contract from the narrower relation
verified by an individual checker.

Zero and constant polynomials have empty factor lists; the unit retains the
input constant. This convention lets reconstruction use the same form for
all inputs.

## Other domains and several variables

`factor_finite_field` accepts arbitrary prime-field polynomials, extracts
the leading scalar, and returns irreducibles with their multiplicities.
`squarefree` handles `ModInt` directly, including polynomials with a zero
derivative. `factor_mod_p` remains the lower-level splitter for monic
squarefree inputs.

`factor_gaussian` works over `QQi`: for example, $x^2+1$ becomes $(x-i)(x+i)$.
The coefficients remain exact Gaussian rationals throughout the computation.
Rounded and ball coefficients do not enter these exact factorization routines.

For sparse multivariate inputs, use `factor_multivariate` for `QQ`,
`factor_multivariate_integer` for `ZZ`, `factor_multivariate_gaussian` for
`QQi`, or `factor_multivariate_finite_field` for `ModInt`.
`squarefree_multivariate` groups multiplicities without splitting every
irreducible. In the coefficient-array API, `sparse.factor` and
`sparse.squarefree` infer the corresponding field and choose these routines.

The Gaussian and prime-field sparse factorization paths use a dense
substitution image and exact subset recombination. Give large inputs a work
budget; sparse input alone does not imply low-memory factorization.

The functional API also factors univariate and sparse multivariate polynomials
over real and complex algebraic numbers. In one variable, `AA` produces linear
factors and real quadratics with negative discriminant; `QQbar` produces only
linear factors. Multivariate irreducibles can have higher degree in either field. See
[algebraic coefficients](algebraic-coefficients.md#factor-over-real-or-complex-algebraic-numbers)
for choosing the field and checking irreducibility.

See the [factorization reference](../reference/p3.md),
[multivariate reference](../reference/advanced-factorization.md), and
[factorization algorithms](../theory/factorization.md).
