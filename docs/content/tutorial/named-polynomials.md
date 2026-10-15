# Named polynomials and expression input

Start with a familiar expression: $(x+1)^3$. We will construct it, evaluate it
at $x=3$, and evaluate its derivative at the same point. With rational
coefficients, every step is exact.

[Build the library](../start/installation.md) first. The following program is
complete and checks its answers against $64$ and $48$.

<!-- example: docs/examples/tutorial_named.mojo -->

## Choose a coefficient ring

`QQ` is the polynomial library's wrapper for apn_mojo `Rational` values. It
provides the coefficient operations that polynomial algorithms need.
`NoContext()` means that this coefficient ring has no extra runtime setting,
such as a finite-field modulus. It does not select numerical precision.

`PolynomialRing[QQ].create(["x"], ...)` names the variables and selects the
coefficient domain. Variable IDs are their positions in the list, so `x` is
variable `0`. A ring created with `["x", "y"]` uses IDs `0` and `1`.
The default monomial order is lexicographic.

## Parse an expression

`ring.parse` expands the expression into a canonical sparse polynomial.
The parser accepts integer literals, variable names, parentheses, unary
signs, addition, subtraction, explicit multiplication, division by scalar
constants, and nonnegative integer powers written with `^`.

Write `2*x`, not implicit multiplication such as `2x`. Over `QQ`, `x/3`
has the exact coefficient $1/3$. Division by another polynomial belongs in
the [division](division.md) or [rational-function](rational-functions.md) APIs;
the expression parser does not perform it. Functions such as `sin(x)` are
outside this polynomial grammar.

## Evaluate and differentiate

`evaluate` takes a list of assignments from variable IDs to coefficients.
The example assigns the rational number `3` to variable `0`. The returned
`QQ` stores the apn_mojo rational in its `.value` field, which is what the
program prints.

`derivative(polynomial, 0, work)` differentiates with respect to `x`. Its
result is another polynomial, which can be evaluated, multiplied, or
substituted like the original. Here $p'(x)=3(x+1)^2$, so $p'(3)=48$.

## Share the work budget

`Budget(1000000)` gives these operations a shared allowance of abstract work
units. Each call consumes part of it. If an operation cannot pay for its
work, it raises a typed error instead of claiming a result it did not finish.
The budget measures work, not seconds or a hard process-memory limit.

`Budget()` uses the largest representable allowance. Start with it when
exploring a small calculation, then choose a finite allowance when your
application needs to bound a request. The [error guide](../guides/errors.md)
explains recovery and the separate limits used by numerical searches.

Next, learn when to use [dense or sparse storage](representations.md), or
continue directly to [division and gcds](division.md).
