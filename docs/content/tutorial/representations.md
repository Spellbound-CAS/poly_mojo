# Dense and sparse polynomials

Choose a representation according to the shape of the polynomial and the
operation you need. Dense storage suits a univariate polynomial with many
occupied degrees. Sparse storage handles several variables and avoids
storing coefficients for absent monomials.

| Representation | Stores | Typical use |
|---|---|---|
| `DensePoly[C]` | A variable ID and coefficients in increasing degree | Univariate division, factorization, interpolation, and roots |
| `SparsePoly[C]` | Nonzero coefficient/monomial pairs and an explicit order | Multivariate arithmetic, substitution, Gröbner bases, and CAD |
| `NumericDensePoly[D]` | Dense coefficients in a numerical domain | Rounded or interval univariate calculations |
| `NumericSparsePoly[D]` | Sparse terms in a numerical domain | Multivariate numerical evaluation and arithmetic |

The [coefficient guide](../guides/domains.md) explains `C` and `D`.

## Read dense coefficients from the constant term upward

The coefficient list `[1, 2, 3]` represents $1+2x+3x^2$. It is not ordered
from the leading coefficient downward. Construction removes trailing exact
zeros; the zero polynomial has an empty coefficient list and degree `-1`.
A zero or constant polynomial still carries its variable and coefficient
context.

Use `DensePoly.from_coefficients` to construct a dense value. When a named
expression is easier to read, parse it with `PolynomialRing` and convert
with `DensePoly.from_sparse`. Supply `variable` explicitly for constants and
zero, whose terms do not identify a variable. Conversion rejects a polynomial
that depends on another variable.

The [calculus example](calculus.md) uses this route. It builds an expression
with named variables, converts it to dense storage, then integrates it.

## Keep sparse terms canonical

Sparse construction combines repeated monomials, removes zero coefficients,
and orders the remaining terms. Use `SparsePoly.from_terms`, the named ring
builder, or other public factories to establish those invariants.

A monomial records nonzero exponents as variable/exponent pairs. Variable IDs
identify variables independently of their printed names. Combining values
requires compatible coefficient contexts and monomial orders; two independent
name lists do not make their same-numbered IDs different variables.

The `terms` field is public for inspection. If you modify it directly, you
are responsible for preserving its canonical form. Prefer construction and
arithmetic methods when building a new value.

## Choose a monomial order deliberately

| Order | Compares first | Useful for |
|---|---|---|
| `LEX` | Exponents in variable-priority order | Elimination and lexicographic Gröbner bases |
| `GRLEX` | Total degree, then lexicographic order | Degree-oriented polynomial calculations |
| `GREVLEX` | Total degree, then a reverse comparison | Gröbner computations where graded reverse lexicographic order is suitable |
| `KERNEL` | The kernel's canonical storage convention | Exact kernel interchange |

Sparse terms are stored in increasing order, so the leading term is last.
`KERNEL` is a storage order, not an admissible Gröbner order. Use LEX, GRLEX,
or GREVLEX for Gröbner bases and ordered list division. `reorder` makes an
explicit conversion while preserving the polynomial's coefficients.

## Convert without changing the domain

`DensePoly.from_sparse` and `DensePoly.to_sparse` change storage, not the
coefficient domain. Promoting integer coefficients to rational coefficients
is a separate operation: use `rational_dense` or `rational_sparse` when a
field operation requires fractions.

Numerical conversion is also explicit. `numeric_from_dense` and
`numeric_from_sparse` convert exact inputs into a chosen numerical domain
with a context. A conversion to a point type rounds coefficients; a conversion
to a ball type encloses them. See [numerical calculations](numerical.md).

For the canonical invariants and shared internal operations, see
[representations in the architecture](../architecture/representations.md).
