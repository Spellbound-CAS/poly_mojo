# Choose a coefficient domain

Start with the mathematical meaning of your coefficients. Exact domains are
suited to identities, divisibility, and symbolic transformations. Numerical
point domains give rounded values at a chosen precision. Ball domains retain
bounds on coefficient uncertainty and rounding error.

## Exact coefficients

| Domain | Values | Construct it from | Choose it for |
|---|---|---|---|
| `ZZ` | Arbitrary-precision integers | Native integers or `Integer` | Integral arithmetic and integer factorization |
| `QQ` | Exact rational numbers | `Rational` | Fractions, exact division, roots, and general symbolic work |
| `GF(p)` | Residues in a prime field | `poly.modint(value, p)` or `modulus=p` | Polynomial arithmetic over a finite field |
| `QQi` | Gaussian rationals, with rational real and imaginary parts | `ExactComplex` | Exact complex arithmetic, factorization and Gröbner bases |
| `AA` | Real algebraic numbers | `poly.algebraic.Real` | Exact irrational coefficients and real roots |
| `QQbar` | Complex algebraic numbers | `poly.algebraic.Complex` | Exact complex coefficients with algebraic components |

`poly` infers the domain from the coefficient types; the
[array tutorial](../tutorial/arrays.md#let-coefficient-types-choose-the-domain)
lists the promotion rules. Indexing and evaluation return apn_mojo `Integer`,
`Rational`, and `ExactComplex` values for `ZZ`, `QQ`, and `QQi`. A `ModInt`
carries its modulus; coefficients from different moduli cannot be mixed in
one polynomial. The supported prime modulus is below $2^{62}$.

`ZZ` is a ring; the other exact domains are fields, which support division by
nonzero coefficients. This difference explains why an operation over `ZZ` may
promote to `QQ`: integrating $x$ requires a coefficient of $1/2$, and ordinary
polynomial division may also introduce fractions. The array functions perform
that promotion. Use `divide_exact` or `pseudo_divmod` when the result must stay
in `ZZ`.

### Match the algorithm to the domain

| Operation | Coefficient scope |
|---|---|
| Arithmetic, calculus, and evaluation | `ZZ`, `QQ`, `GF(p)`, `QQi`, `AA`, `QQbar` |
| Dense and sparse factorization | `ZZ`, `QQ`, `QQi`, `GF(p)`, `AA`, `QQbar`, all with repeated factors |
| Squarefree decomposition | `QQ`, `QQi`, `GF(p)`, `AA`, `QQbar`; integer inputs promote to `QQ` |
| Gröbner bases and membership | `QQ`, `QQi`, `GF(p)`, `AA`, `QQbar`; integer inputs promote to `QQ` |
| Exact real roots | `ZZ`, `QQ`, `AA` |
| Exact complex roots | `ZZ`, `QQ`, `QQi`, `AA`, `QQbar`, with multiplicities |
| Rational functions | `QQ`, `QQi`, `GF(p)`, `AA`, `QQbar`; integer inputs promote to `QQ` |
| CAD | Integer, rational, and `AA` polynomial coefficients |
| Algebraic serialization | `AA`, `QQbar`, with checked exact root witnesses |
| Algebraic numerical conversion and certified evaluation | `AA`, `QQbar`, with explicit precision contexts |
| Structured bases and formal series | Exact fields, subject to invertible conversion factors |

`poly` infers `AA` and `QQbar` from immutable algebraic values and supports
factorization, rational functions, serialization, and explicit numerical
conversion in both fields. CAD accepts real coefficients, including AA; complex
coefficients have no ordered sign interpretation. See [algebraic coefficients](../tutorial/algebraic-coefficients.md)
for examples and costs.

## Numerical coefficients

| Domain | apn_mojo value | Meaning |
|---|---|---|
| `FloatDomain` | `Float` | A stored real point value |
| `ComplexDomain` | `Complex` | A stored complex point value |
| `BallDomain` | `Ball` | A real interval |
| `ComplexBallDomain` | `ComplexBall` | A rectangle in the complex plane |

Numerical arrays store one of these domains. Build them with `dense.asarray`
or `sparse.asarray` and an explicit apn_mojo context: `ArithmeticContext`,
`ComplexContext`, or `BallContext` (for either ball domain). Every numerical
operation takes `context=`. The context controls the arithmetic performed by
the operation. It cannot recover information already lost when the input
coefficients were rounded.

Point operations round at their documented arithmetic steps. They do not
claim a correctly rounded exact result for an entire multi-step polynomial
calculation. Ball operations enclose the represented coefficient family;
a wider result can reflect both input uncertainty and accumulated enclosure
error.

Increasing working precision can reduce arithmetic error. It cannot remove
uncertainty supplied in the coefficients. An uncertain leading coefficient
can even leave the degree undecided. See [numerical root covers](../tutorial/array-roots.md#choose-the-numerical-claim-you-need)
for interpreting such outcomes.

## Keep exact and numerical workflows explicit

Pass an exact array to `asarray` with a `context=`, or use `astype`, to obtain
numerical coefficients; `astype` also changes a numerical array's domain or
precision. There is no implicit route from rounded coefficients into exact
gcds or Gröbner bases.

You can also retain the original exact array and call `dense.enclose` or
`sparse.enclose` to evaluate it with rigorous bounds. This is useful when
the exact coefficients are known and only the evaluation point or the desired
output is numerical.

For AA and QQbar arrays, point conversion resolves the requested rounding,
while ball conversion produces an outward enclosure. `dense.enclose` and
`sparse.enclose` retain coefficient-witness evidence alongside their numerical
trace. See [numerical arrays](../tutorial/array-numeric.md) and
[certified evaluation](../tutorial/array-certified.md).
