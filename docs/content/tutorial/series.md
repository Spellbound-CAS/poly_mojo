# Formal power series

For coefficient iterables and functional calls, start with the [array interface](array-series.md).

A truncated series records coefficients below a known order. An order of
$n$ means that calculations take place modulo $x^n$: coefficients at degree
$n$ and above are discarded, without an assumption about their values.

## Compute a reciprocal

Modulo $x^4$, the reciprocal of $1-x$ is $1+x+x^2+x^3$. Multiplying them
gives $1-x^4$, whose retained part is exactly $1$.

<!-- example: docs/examples/tutorial_series.mojo -->

`from_polynomial` keeps the coefficients below the requested order. The
returned `TruncatedSeries` carries both the polynomial and the order.
`.inverse` computes the reciprocal to that same order.

This is a formal statement about coefficients. It provides no analytic tail
bound for evaluating an infinite series at a nonzero argument. Use interval
methods and a separately justified remainder bound when that is the claim
you need.

## Keep the order explicit

Binary series operations require equal orders. Use `.truncate` to reduce a
series to a common order; truncation cannot recover unknown higher-order
coefficients. Differentiation lowers the known order by one, and integration
raises it by one with zero integration constant.

| Operation | Required constant or leading coefficient |
|---|---|
| Reciprocal or division | The divisor's constant coefficient is invertible |
| Composition | The inner series has exactly zero constant coefficient |
| Logarithm | Constant coefficient is one |
| Exponential | Constant coefficient is zero |
| Compositional inverse (`revert`) | Constant coefficient is zero and linear coefficient is nonzero |

Exact series use a coefficient field. Finite characteristic can prevent
required divisions by integers, so a valid constant coefficient alone does
not guarantee that every requested series operation exists at every order.

## Use numerical series deliberately

`NumericSeries` carries the same order information with numerical
coefficients. Its operations take a context. Point results round;
ball results enclose the operations on the represented coefficients.
The formal-order convention still supplies no analytic tail estimate.

See the [series reference](../reference/toolkit.md#formal-series) for the full
API, including differentiation, integration, logarithms, exponentials, and
reversion.
