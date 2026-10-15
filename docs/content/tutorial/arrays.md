# Calculate with polynomial arrays

Use `poly.dense` for ascending coefficients and `poly.sparse` for nonzero terms.
Both infer exact coefficient domains and return immutable, iterable results.
The same function names cover arithmetic, powers, calculus, division, and
evaluation in either representation.

## Start with coefficients or terms

<!-- example: docs/examples/array_api.mojo -->

`[1, 2, 3]` represents $1+2x+3x^2$. Sparse terms `[(1, 2), (2, 3)]`
represent $2x+3x^2$. For several variables, use power tuples:
`[((1, 0), 2), ((0, 2), 3)]` represents $2x+3y^2$.
Tuple positions identify variables; every term must use the same power shape.
Sparse construction combines repeated powers and removes zero coefficients.

Lists, fixed arrays, spans, and finite iterables work as inputs. Dense
construction also accepts one-dimensional APN integer and rational batches.
An owning iterator can be consumed with `asarray(iterator^)`. The resulting
polynomial stays in memory; iteration does not imply constant-memory storage.

## Let coefficient types choose the domain

| Coefficients | Domain | Values returned by indexing and evaluation |
|---|---|---|
| Native integers or `Integer` | `ZZ` | `Integer` |
| Any `Rational`, with optional integers | `QQ` | `Rational` |
| Any `ExactComplex`, with optional integers or rationals | `QQi` | `ExactComplex` |
| `ModInt`, with optional integers | A common prime field | `ModInt` |
| `algebraic.Real`, with optional integers or rationals | `AA` | `algebraic.Real` |
| `algebraic.Complex`, or `Real` mixed with `ExactComplex` | `QQbar` | `algebraic.Complex` |

Promotion follows types, including zero coefficients. A rational zero still
selects `QQ`; an exact complex value with zero imaginary part still selects
`QQi`. Algebraic arrays likewise retain `AA` or `QQbar`, even when a result
is rational or zero. See [algebraic coefficients](algebraic-coefficients.md).
Floats and `Bool` require an explicit conversion before they can enter
this exact interface.

Mojo must choose a list's element type before calling a function. When scalar
types differ, use `poly.coefficients((...))`, as above. For sparse inputs,
`poly.sparse.asarray(((0, 3), (2, Rational(1, 2))))` directly accepts a tuple of
mixed terms. Prepare that result before passing it to an arithmetic function.

`poly.modint(value, prime)` constructs a finite-field coefficient.
Alternatively, `asarray(values, modulus=7)` reduces integer coefficients into
$\mathbb F_7$. All modular inputs must agree on the prime; rational and complex
values have no implicit modular conversion. An empty modular input needs
`modulus=` unless it is already a prepared result retaining its field.

## Divide, integrate, and evaluate at many points

<!-- example: docs/examples/array_operations.mojo -->

`divmod` returns `(quotient, remainder)`. Integer coefficients promote to exact
rationals, so dividing $x^2+1$ by $2x+1$ returns $x/2-1/4$ and $5/4$.
`integral` makes the same promotion and chooses zero as the integration
constant. In a prime field, integration raises `DOMAIN` when a nonzero term
would require division by zero.

`evaluate_many` evaluates one polynomial at each point, in input order. Its
`Values` result retains all zeros, including trailing zeros. For a multivariate
sparse polynomial, each point is a tuple of coordinates. The polynomial is
prepared and promoted once for the call. Use [batch maps](array-batching.md)
for broadcasting collections of polynomials or constructing Cartesian grids.

`from_roots` constructs the monic product of $(x-r)$ factors, keeping repeated
and zero roots. It accepts known exact roots; use `real_roots` to find them from coefficients
or a [root cover](array-roots.md) for certified numerical bounds.

## Work with scalar constants

Pass a scalar directly to `add`, `subtract`, or `multiply`, in either operand
position. `dense.subtract(4, p)` computes $4-p$. `divide(p, scalar)` scales by
the scalar's reciprocal; integer coefficients promote to `QQ`, so no precision
is lost. Use `divmod` for polynomial quotient and remainder, or
`poly.rational` to retain a polynomial denominator.

<!-- example: docs/examples/array_surface_scalars.mojo -->

Exact coefficients and scalars follow the usual domain promotion rules.
Prime-field operands must share a modulus and can mix with integers.
Division by zero raises `DOMAIN`, including for a zero polynomial. Prepared
inputs remain unchanged.

For numerical arrays, pass `context=` to each operation. The scalar converts
to the array's stored domain; use `astype` first if you need a complex domain.
A ball divisor must exclude zero. Sparse scaling visits stored terms without
allocating coefficients for missing powers.

## Display variable names

Keep coefficients as data and choose names when displaying a result:

<!-- example: docs/examples/array_surface_format.mojo -->

`format` supports `ZZ`, `QQ`, `QQi`, and prime fields. It prints every nonzero
term with its exact coefficient in decreasing lexicographic order. Dense
polynomials and one-axis sparse polynomials default to `x`; multiple axes
default to `x0`, `x1`, and so on. Supply one distinct ASCII identifier per
sparse axis, including unused axes. The name `i` is reserved over `QQi`.

Display text omits domain, modulus, and shape metadata. Use
[serialization](array-serialization.md) to save a value, including algebraic
or numerical coefficients.

## Choose the operation

| Task | Dense | Sparse |
|---|---|---|
| Construct or reuse | `asarray` | `asarray` |
| Add, subtract, multiply (polynomial or scalar) | `add`, `subtract`, `multiply` | Same |
| Divide by a scalar | `divide` | Same; preserves sparse storage |
| Display named variables | `format(..., variable="t")` | `format(..., variables=["u", "v"])` |
| Negate or raise to a nonnegative power | `negative`, `power` | Same |
| Highest power | `degree` | `degree`; use `axis=` for one variable |
| Differentiate or integrate | `derivative`, `integral` | Same, with `axis=` |
| Quotient and remainder | `divmod` | `divmod` under lexicographic order |
| Evaluate | `evaluate`, `evaluate_many` | Same; coordinate tuples for several variables |
| Substitute polynomials | `compose` | `compose` for one axis; `substitute` for several |
| Exact algebra | `gcd`, `extended_gcd`, `content`, `primitive_part`, `resultant`, `subresultants` | Same; some operations take `axis=` or require one axis |
| Exact division or several divisors | `divide_exact`, `pseudo_divmod`, `divide_many` | Same; pseudo-division requires one axis |
| Irreducible factors or multiplicity blocks | `factor`, `squarefree` | `factor`; `squarefree` requires one axis |
| Recover from samples | `interpolate`, `hermite_interpolate`, `prepare_points` | Convert the dense result with `poly.to_sparse` |
| Use core algorithms | `to_core`, `from_core` | Same, with explicit variable maps |
| Find exact real roots | `real_roots`, `count_roots`, `root_between` | Same, for one variable |
| Certified root covers | `real_root_cover`, `complex_root_cover` | Same, for one variable |
| Prepare and query ideals | Use sparse polynomials | `groebner`, `reduce`, `contains`, `membership` |
| Construct from known roots | `from_roots` | Convert the dense result with `poly.to_sparse` |

`poly.to_sparse` drops zero coefficients. `poly.to_dense` fills missing powers
of a univariate polynomial with zeros. Both preserve the domain and modulus.
Dense storage needs one slot for every power through the degree; set
`max_steps` when converting a sparse input with potentially large powers.

The zero polynomial has degree `-1`. Raising any polynomial, including zero,
to power zero returns one. Sparse degree defaults to total degree; sparse
output terms use increasing lexicographic order with tuple position zero
having the highest comparison priority.

For precise signatures and shared conventions, see the
[array API reference](../reference/array-api.md); [budgets and errors](../guides/errors.md)
explains work limits. Continue with [exact array algebra](array-exact.md) for factorization,
interpolation, elimination, and core interoperability. The `poly` namespaces
also expose root finding, Gröbner bases, CAD, and their certificate checks.

Continue with [root finding and algebraic values](array-roots.md), or use
[sparse equations and inequalities](array-systems.md) for Gröbner bases and CAD.

Use [formal series](array-series.md) for truncated calculations,
[basis arrays](array-bases.md) for Newton, Chebyshev, or Bernstein coefficients,
and [rational functions](array-rational.md) for normalized polynomial ratios.
Each representation retains the metadata needed by its functional operations.

For floating-point or ball coefficients, choose an [explicit numerical context](array-numeric.md). You can also [save canonical coefficient records](array-serialization.md).
