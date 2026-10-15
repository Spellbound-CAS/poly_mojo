# Work in another polynomial basis

Use `poly.basis` when your coefficients belong to a Newton, Chebyshev, or
Bernstein basis. Construct a representation directly, or convert a polynomial
from its ordinary ascending coefficients. Evaluation uses the retained basis;
arithmetic and calculus return a new representation.

## Represent and subdivide a quadratic

The polynomial $x^2$ has Chebyshev coefficients $(1/2,0,1/2)$ because
$T_2(x)=2x^2-1$. Its degree-two Bernstein coefficients on $[0,1]$ are $(0,0,1)$.

<!-- example: docs/examples/array_toolkit_basis.mojo -->

`basis.asarray(values, basis=...)` interprets the supplied values as coefficients
in that basis. `basis.from_polynomial(values, basis=...)` converts monomial
coefficients. Both are useful, but they answer different questions. The example
uses conversion for Chebyshev and direct construction for Bernstein.

Splitting the Bernstein polynomial at $1/2$ returns representations of
$p(x/2)$ and $p(1/2+x/2)$. Each uses its own parameter on $[0,1]$. They meet at
$1/4$, and both retain the original representation degree.

## Retain the data that defines the basis

| Basis name | Required metadata | Evaluation |
|---|---|---|
| `"newton"` | Exactly one fewer node than directly supplied coefficients; repeated nodes are allowed | Nested multiplication |
| `"chebyshev"` | First-kind polynomials $T_k$ | Clenshaw recurrence |
| `"bernstein"` | Coefficient count minus one is the retained degree | de Casteljau evaluation |

For Newton conversion, pass enough `nodes=` to cover the polynomial degree.
Conversion fills any remaining basis coefficients with zero. For Bernstein
conversion, `degree=` can retain a larger degree; omitting it chooses the
polynomial degree, with degree zero for the zero polynomial.

Trailing zeros are retained during direct construction. In particular,
`basis.asarray([0,0,0], basis="bernstein")` is a degree-two representation of
zero. An empty Bernstein coefficient list raises `DOMAIN`.

Read coefficients by indexing or iteration, or use `basis.coefficients` and
`basis.nodes` for independent scalar sequences. `basis.to_polynomial` returns
a canonical dense array. Core import and export copy mutable storage.

## Calculate without changing the input

`add`, `subtract`, and `multiply` accept prepared representations of the same
basis kind and promote compatible coefficient fields. `negative`, `derivative`,
and `integral` act on one representation; integration uses a zero constant.
These operations convert through exact monomial arithmetic and return to the
chosen basis. Direct evaluation and subdivision use the core basis recurrences.

Bernstein addition retains the larger degree, multiplication adds the retained
degrees, and calculus adjusts the degree by one, with degree at least zero.
Newton results retain the left operand's nodes. If the result needs a higher
degree, the call raises `DOMAIN`; prepare enough nodes first, or perform the
operation with `poly.dense` and convert back with the desired `nodes=`.
Different basis kinds require an explicit conversion.

Integer input promotes to `QQ`. Rational, Gaussian rational, and modular
coefficients use their inferred fields; Newton nodes participate in that
inference. Finite-field conversions can fail when a required factor is not
invertible. Use an explicit `max_steps=` or a shared `Budget` to bound work,
including conversion and output allocation.

The [basis API reference](../reference/array-basis.md) lists every operation.

For floating-point or ball coefficients, pass `context=` when constructing the
representation and when computing with it. The [numerical series and basis
example](array-numeric.md#numerical-series-and-alternative-bases) covers both workflows.
