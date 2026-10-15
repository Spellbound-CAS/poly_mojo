# Compute with formal series

Use `poly.series` to calculate a fixed number of exact coefficients without
expanding higher powers. A series retains its order, so a reciprocal,
composition, or derivative carries the precision information into the next call.

## Expand a reciprocal and an exponential

The reciprocal of $1-x$ begins $1+x+x^2+\cdots$. At order six, you know
coefficients through $x^5$ and leave the rest as $O(x^6)$.

<!-- example: docs/examples/array_toolkit_series.mojo -->

`asarray([1, -1], order=6)` treats the input as a polynomial and supplies zeros
for its missing coefficients below order six. `reciprocal` returns six ones.
The product is $1+O(x^6)$. For $\exp(x)$, the coefficient of $x^5$ is exactly
$1/120$; differentiation reduces the known order from six to five.

## Keep the unknown tail visible

A prepared series supports indexing, iteration, and `len`. Its length is its
order, including known trailing zeros. The backing polynomial stays compact;
`series.coefficients(s)` materializes all known coefficients under the work
budget. Use `series.to_polynomial(s)` when you deliberately want to discard
the unknown tail and work with the known polynomial part.

Raw coefficient iterables require `order=`. Operations on a prepared series
reuse its order. Two prepared inputs must have equal orders unless you pass
an explicit lower `order=`. `series.truncate(s, n)` lowers the order, and no
operation accepts a request to recover unknown higher coefficients.

This is a formal coefficient calculation. It supplies no analytic error bound
for evaluating an infinite series away from zero.

## Choose an operation

| Function | Required input |
|---|---|
| `add`, `subtract`, `multiply`, `negative` | Compatible exact fields and a common order |
| `divide`, `reciprocal` | An invertible divisor constant |
| `compose` | A zero-constant inner series and a known constant term |
| `derivative` | Returns order `max(0, n-1)` |
| `integral` | Returns order `n+1` with zero constant; required divisions must exist |
| `log` | A known constant equal to one |
| `exp` | A known constant equal to zero |
| `revert` | Order at least two, zero constant, and nonzero linear coefficient |

Series algorithms use fields. Integer input promotes to `QQ`; rational,
Gaussian rational, and modular inputs retain their inferred fields. Mixed
scalar sequences use `poly.coefficients`, as in the [array tutorial](arrays.md).
Modular coefficients can mix with raw integer coefficients. A prepared `QQ`
series does not implicitly convert into a prime field. Different moduli raise
`DOMAIN`, including for zero series.

Finite characteristic can prevent the divisions needed by integration,
logarithms, and exponentials. A nonzero constant alone does not guarantee
that every operation exists at the requested order.

Every calculation accepts `max_steps=` or a shared `Budget`. Input collection,
promotion, the core algorithm, and result storage consume the same allowance.
A failure leaves prepared inputs usable. `to_core` and `from_core` copy storage
for interoperability with the mutable core API.

See the [series API reference](../reference/array-series.md) for signatures and
[numerical series](series.md#use-numerical-series-deliberately) for context-based
floating-point and ball operations.

For floating-point or ball coefficients, pass `context=` when constructing the
representation and when computing with it. The [numerical series and basis
example](array-numeric.md#numerical-series-and-alternative-bases) covers both workflows.

## Map a core variable

`to_core(p, variable=42)` exports an independent core value with variable ID
42. `from_core` accepts any core variable ID and normalizes it to array axis
zero. Both exact and numerical values support this boundary. Numerical
imports take `context=` and preserve the stored coefficients' formats and
bounds. See the [representation example](array-numeric.md#exchange-representation-metadata).
