# Choose a coefficient domain

Start with the mathematical meaning of your coefficients. Exact domains are
suited to identities, divisibility, and symbolic transformations. Numerical
point domains give rounded values at a chosen precision. Ball domains retain
bounds on coefficient uncertainty and rounding error.

## Exact coefficients

| Domain | Values | Context | Choose it for |
|---|---|---|---|
| `ZZ` | Arbitrary-precision integers | `NoContext` | Integral arithmetic and integer factorization |
| `QQ` | Exact rational numbers | `NoContext` | Fractions, exact division, roots, and general symbolic work |
| `ModInt` | Residues in a prime field | `Modulus` | Polynomial arithmetic over a finite field |
| `QQi` | Gaussian rationals, with rational real and imaginary parts | `NoContext` | Exact complex arithmetic, factorization and Gröbner bases |
| `AA` (`poly.algebraic.Real`) | Real algebraic numbers | `NoContext` | Exact irrational coefficients and real roots |
| `QQbar` (`poly.algebraic.Complex`) | Complex algebraic numbers | `NoContext` | Exact complex coefficients with algebraic components |

`ZZ`, `QQ`, and `QQi` wrap the corresponding apn_mojo types. `ModInt` carries
its modulus; coefficients from different moduli cannot be mixed in one ring.
The supported prime modulus is below $2^{62}$.

`CoefficientRing` supplies exact arithmetic. `CoefficientField` additionally
supports division by nonzero coefficients. This difference explains why an
operation over `ZZ` may need promotion to `QQ`: integrating $x$ requires a
coefficient of $1/2$, and ordinary polynomial division may also introduce
fractions. Use the integer-specific promotion helpers or pseudo-division
when appropriate.

### Match the algorithm to the domain

| Operation | Coefficient scope |
|---|---|
| Exact container arithmetic | `ZZ`, `QQ`, `ModInt`, `QQi`, `AA`, `QQbar` |
| Dense and sparse factorization | `ZZ`, `QQ`, `QQi`, prime `ModInt`; `AA` and `QQbar` through `poly`, all with repeated factors |
| Squarefree decomposition | `QQ`, `QQi`, prime `ModInt`, `AA`, `QQbar`; the array API promotes `ZZ` to `QQ` |
| Gröbner bases and membership | `QQ`, `QQi`, prime `ModInt`, `AA`, `QQbar` |
| Exact complex roots through `poly` | `ZZ`, `QQ`, `QQi`, `AA`, `QQbar`, with multiplicities |
| Core Sturm isolation | `QQ` |
| Rational functions through `poly` | `QQ`, `QQi`, `AA`, `QQbar`; integer inputs promote to `QQ` |
| CAD through `poly` | Integer, rational, and `AA` polynomial coefficients |
| Algebraic serialization | `AA`, `QQbar`, with checked exact root witnesses |
| Algebraic numerical conversion and certified evaluation | `AA`, `QQbar`, with explicit precision contexts |
| Structured bases and formal series | Exact fields, subject to invertible conversion factors |

The array API infers `AA` and `QQbar` from immutable algebraic values. Their types
also work with field-generic core containers and algorithms. The facade extends
algebraic factorization, rational functions, serialization, and explicit numerical
conversion to both fields. CAD accepts real coefficients, including AA; complex
coefficients have no ordered sign interpretation. See [algebraic coefficients](../tutorial/algebraic-coefficients.md)
for examples and costs.

## Numerical coefficients

| Domain | apn_mojo value | Meaning |
|---|---|---|
| `FloatDomain` | `Float` | A stored real point value |
| `ComplexDomain` | `Complex` | A stored complex point value |
| `BallDomain` | `Ball` | A real interval |
| `ComplexBallDomain` | `ComplexBall` | A rectangle in the complex plane |

Use these domains with `NumericDensePoly` and `NumericSparsePoly`. Numerical
operations take an explicit apn_mojo context: `ArithmeticContext`,
`ComplexContext`, or `BallContext` (for either ball domain).
The context controls the arithmetic performed by the operation. It cannot
recover information already lost when the input coefficients were rounded.

Point operations round at their documented arithmetic steps. They do not
claim a correctly rounded exact result for an entire multi-step polynomial
calculation. Ball operations enclose the represented coefficient family;
a wider result can reflect both input uncertainty and accumulated enclosure
error.

Increasing working precision can reduce arithmetic error. It cannot remove
uncertainty supplied in the coefficients. An uncertain leading coefficient
can even leave the degree undecided. See [numerical roots](../tutorial/numerical-roots.md)
for interpreting such outcomes.

## Keep exact and numerical workflows explicit

Use `numeric_from_dense` or `numeric_from_sparse` to convert an exact
polynomial for numerical work. Use `convert_numeric_dense` or
`convert_numeric_sparse` to change a numerical domain or context under that
conversion's contract. There is no implicit route from rounded coefficients
into exact gcds or Gröbner bases.

You can also retain the original exact polynomial and call `enclose_dense`
or `enclose_sparse` to evaluate it with rigorous bounds. This is useful when
the exact coefficients are known and only the evaluation point or the desired
output is numerical.


For AA and QQbar arrays, use `dense.asarray` or `sparse.asarray` with an explicit
numerical `context=`. Point conversion resolves the requested rounding, while
ball conversion produces an outward enclosure. `dense.enclose` and
`sparse.enclose` retain coefficient-witness evidence alongside their numerical
trace. See [numerical arrays](../tutorial/array-numeric.md) and
[certified evaluation](../tutorial/array-certified.md).
