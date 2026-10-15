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

Each function accepts an optional `max_steps` limit. As a call computes, the
library counts steps: one for each term or coefficient it processes and for
each element it allocates, and more for arithmetic on large numbers, growing
with their size in bytes. If a call would exceed its limit, it stops and
raises `PolyError` with kind `LIMIT`. The default, `UInt64.MAX`, means no
limit.

The count is exact and repeatable: the same call on the same input takes the
same number of steps on any machine. To see how many steps a call takes, pass
it a `Budget` from `polynomials` as the last argument, and read how much of
the budget it used:

<!-- example: docs/examples/tutorial_steps.mojo -->

Computing $(1+x)^{100}$ takes 329,732 steps, so a limit of 10,000 stops it and
a limit of 1,000,000 lets it finish.

Steps do not convert to time. On one laptop core, different operations ran at
anywhere from about 100 million to over 10 billion steps per second, because
the count charges arithmetic on large numbers more than it costs to run.
Measure the steps your own inputs take, and set the limit with room to spare.
To bound a sequence of calculations together, pass the same `Budget` to each
call. [Budgets and errors](../guides/errors.md#what-a-step-costs) lists what
each kind of work costs.

Continue with [dense and sparse arrays](arrays.md).
