# Your first polynomial

A dense polynomial is a sequence of coefficients, starting with the constant
term. `[1, 1]` means $1+x$; `[1, 2, 3]` means $1+2x+3x^2$.
Pass those coefficients to `poly.dense` functions to calculate with exact
integers and fractions, without selecting a ring or parsing an expression.

[Install the package](../start/installation.md), then run this program. It expands
$(1+x)^3$ and evaluates the result and its derivative at $x=3$.

<!-- example: docs/examples/tutorial_first.mojo -->

`power` returns the coefficients `[1, 3, 3, 1]`. The input integers select
integer arithmetic; `evaluate` returns an apn_mojo `Integer`. All intermediate
values and both answers are exact.

## Reuse a polynomial

Every arithmetic function accepts a coefficient iterable. Call
`poly.dense.asarray(values)` when you want to normalize an input once and reuse
it. A prepared result supports `len`, indexing, iteration, `domain()`, and
`to_list()`. Copies share immutable storage; `to_list()` gives you a mutable
list. Trailing zero coefficients are removed, so the zero polynomial has an
empty coefficient sequence.

Rational coefficients select `QQ`; exact complex coefficients select `QQi`.
Use `poly.coefficients((...))` for a tuple of mixed scalar types, or pass an
ordinary homogeneous list directly. The [array tutorial](arrays.md) covers
promotion, sparse terms, modular coefficients, and evaluation at many points.

## Bound a calculation

Each function accepts an optional `max_steps` allowance. For example,
`poly.dense.power([1, 1], 100, max_steps=10000)` raises `PolyError` with kind
`LIMIT` if it cannot finish within that allowance. A default call uses
`UInt64.MAX`. The allowance measures work, not elapsed time or a hard RAM limit.

For a sequence of calculations, import `Budget` from `polynomials` and pass
the same mutable budget as the last positional argument to each call. See
[budgets and errors](../guides/errors.md) for accounting and failure
behavior.

Continue with [dense and sparse arrays](arrays.md). The
[named-polynomial interface](named-polynomials.md) is available when your input
already contains variable names or mathematical expressions.
