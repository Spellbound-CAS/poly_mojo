# Calculus and transformations

Differentiate, integrate, or substitute directly on polynomial coefficients.
Over rational coefficients these operations are exact, including the fractions
introduced by integration.

For coefficient and term inputs through `poly`, see
[exact array algebra](array-exact.md). This walkthrough uses the core types
and their explicit contracts.

## Integrate and differentiate

For $p(x)=3x^2+2x$, an antiderivative with constant term zero is $x^3+x^2$.
Its integral from $0$ to $1$ is exactly $2$.

<!-- example: docs/examples/tutorial_calculus.mojo -->

`integrate_dense` returns the antiderivative with zero integration constant.
The assertion differentiates it back to the original polynomial. To select
another constant, add a constant polynomial afterward.

An integer polynomial may need rational coefficients after integration.
Use `integrate_integer_dense` or `integrate_integer_sparse` for the explicit
promotion from `ZZ` to `QQ`. Over a finite field, a required coefficient
division must be invertible; an antiderivative is not available in every
characteristic.

For sparse polynomials, `derivative` and `integrate_sparse` take the variable
ID. Differentiating with respect to a variable absent from a polynomial gives
zero. [Variable IDs and storage](representations.md) explains the common
representation rules.

## Substitute another polynomial

`compose_dense(outer, inner, work)` computes $\text{outer}(\text{inner}(x))$.
The example substitutes $x+1$ into $3x^2+2x$ and checks the expansion
$3x^2+8x+5$.

Use `affine_transform` for $p(ax+b)$ and `taylor_shift` for $p(x+a)$.
`jet_evaluate` computes a local Taylor jet. Its entries are Taylor/Hasse
coefficients: the coefficient of $t^k$ in $p(a+t)$, rather than the raw
$k$th derivative. In characteristic zero, multiply by $k!$ to obtain that
derivative.

`substitute_many` performs simultaneous sparse substitution. Each replacement
is interpreted against the original variables, so replacing `x` with an
expression containing `y` does not apply a second replacement inside the
first one. This matters for transformations such as swapping two variables.

## Use numerical coefficients

The numerical counterparts, such as `numeric_integrate_dense`,
`numeric_compose_dense`, and `numeric_substitute_many`, take an explicit
apn_mojo context. Point domains round their arithmetic; ball domains enclose
the represented inputs. They do not inherit exact identities merely because
the corresponding symbolic operation is exact.

See the [calculus reference](../reference/toolkit.md#calculus-and-transformations)
for signatures and the [interpolation tutorial](interpolation.md) for building
a polynomial from values or derivative data.
