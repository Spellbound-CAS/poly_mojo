# Work with floating-point and ball coefficients

Use the same ascending coefficients and sparse power tuples for numerical
polynomials. Supply a precision context when you construct them, then pass the
context to each arithmetic operation. Choose floating-point coefficients for
rounded calculations or balls when you need an enclosure of every possible
result.

## Choose the arithmetic

<!-- example: docs/examples/array_numeric.mojo -->

The first calculation evaluates $1+2x+3x^2$ at $x=2$, giving 17. The second
starts from the exact rational $1/3$ and encloses its conversion and evaluation
in a ball. A ball's width includes rounding error and any uncertainty you put
in the coefficients.

| Context | Coefficient values |
|---|---|
| `ArithmeticContext(format=FloatFormat(precision=80))` | APN `Float` |
| `ComplexContext(real_context)` | APN `Complex` |
| `BallContext(precision=80)` | APN `Ball`, or `ComplexBall` for complex input |

Construct a context explicitly; a `BallContext` without a precision raises
`DOMAIN`. Without `context=`, `dense.asarray` and `sparse.asarray` retain their
exact integer, rational, Gaussian rational, algebraic, and prime-field behavior.

Numerical constructors accept exact integers, rationals, real and complex
algebraic values, native floating-point values, and APN numerical scalars. A native float contributes its stored binary
value: converting `Float64(0.1)` does not recover the decimal rational $1/10$.
Use `Rational(1, 10)` when that fraction is your intended input. Modular values
and booleans do not convert implicitly.

For dense mixed scalar types, pass a tuple, for example
`dense.asarray((1, Rational(1, 3)), context=context)`. An explicitly complex
context can also embed real coefficients. Ball construction encloses conversion
error; converting a ball to a point type requires all enclosed values to round
to the same point.

## Choose the destination type

Use `dense.astype[Float](values, context=...)` or its `sparse` counterpart when
you want an explicit coefficient type. The supported destinations are APN
`Float`, `Complex`, `Ball`, and `ComplexBall`, paired with the contexts in the
table above. Both ball types use `BallContext`.

<!-- example: docs/examples/array_access_dtype.mojo -->

The first conversion rounds the stored 100-bit coefficients to 30 bits and
leaves the original array unchanged. The second requests complex balls even
though all input coefficients are real. Raw homogeneous sequences, mixed-type
tuples, and prepared arrays use the same calls. Sparse mixed input is a tuple
of `(power, coefficient)` pairs with a common power shape.

`astype` converts every coefficient, including within the same domain. A real
destination requires exactly zero imaginary parts; an uncertain imaginary
interval is rejected even when it contains zero. A ball-to-point conversion
succeeds only when the entire enclosure rounds to one value. Ball destinations
round bounds outward, preserving uncertainty. Choose an exact rational input
when you want conversion from that exact value instead of from an already
rounded approximation.

## Convert algebraic coefficients

Pass an `AA` or `QQbar` array to `asarray` with a context to obtain numerical
coefficients. The exact input stays unchanged. Real and complex floating-point
conversion refines exact root bounds until every enclosed value rounds to the
same result with the same rounding status. It honors the context's rounding
mode and traps, and raises `DOMAIN` if `max_precision` cannot establish a
finite rounded result. Work exhaustion separately raises `LIMIT`.

<!-- example: docs/examples/array_algebraic_numeric.mojo -->

With `BallContext`, conversion refines the real components at a scale relative
to their magnitude, then rounds their bounds outward. `QQbar` selects
`ComplexBall`; an explicitly real context rejects a nonzero imaginary part.
Complex contexts can choose separate real and imaginary precisions.

A numerical array contains rounded coefficients or coefficient enclosures.
Use `enclose` or `adaptive_evaluate` on the exact array when you also need an
independently checkable connection to the original algebraic coefficients.
Those results retain root witnesses alongside their arithmetic traces; see
[certified evaluation](array-certified.md).

## Calculate with prepared arrays

Numerical arrays support arithmetic, calculus, evaluation, interpolation, and
transformations. Dense and univariate sparse arrays also support composition
and division.
Every call takes `context=`. Sparse calculus uses `axis=`, and a multivariate
evaluation takes a coordinate tuple in axis order.

Construction from raw values applies the context. Calling `asarray` on an
already prepared numerical array of the selected domain shares its stored
coefficients; it does not round them again. Arithmetic applies the context you
pass to that operation. Use `astype` to convert every stored coefficient into
another type or precision context.

Arrays support indexing, iteration, and `to_list()`. `poly.to_dense` and
`poly.to_sparse` convert numerical storage with an explicit context. Dense
conversion allocates gap zeros, so use a work limit when sparse exponents can
be large.

## Exchange numerical polynomials with the core

`to_core` and `from_core` make independent copies while preserving every stored
numerical format and ball bound. They need no arithmetic context. Dense export
accepts `variable=`; import normalizes that ID to the array's axis zero.

<!-- example: docs/examples/array_access_interop.mojo -->

For sparse polynomials, `variables=` maps array axes to core variable IDs in
both directions. Use `from_core[ndim=N]` to preserve the full power shape,
including unused axes. Import accepts any valid core monomial order and returns
increasing lexicographic array storage. Signed IDs and a different axis mapping
are allowed; duplicate map entries or an unmapped used variable raise `DOMAIN`.
Changing the exported core polynomial leaves the original array unchanged.

## Distinguish a degree bound from a known degree

Only certain zero coefficients are removed. If the leading ball contains zero
and nonzero values, `degree_bound(p)` returns its power, while
`certified_degree(p)` returns `None`. A zero polynomial has degree -1.
This distinction prevents uncertain leading coefficients from disappearing
during construction or arithmetic.

Numerical root covers use `dense.real_root_cover` or
`dense.complex_root_cover`. Inspect `complete()`, `accurate()`, `isolated()`,
and `outcome()` before using a cover. `roots.check` verifies its claims, and
`roots.refine` returns a new cover under new `NumericOptions`. Ball coefficients
describe a family of polynomials; a method can leave cells unresolved when one
uniform root claim does not hold across that family.

All computations accept `max_steps=` or a shared `Budget`, including input
conversion and output allocation. See [root workflows](array-roots.md),
[saving arrays](array-serialization.md), and the
[numerical array reference](../reference/array-numeric.md).


## Interpolate, transform, and divide

Numerical `interpolate` and `hermite_interpolate` accept ordinary sequences of
nodes and samples with `context=`. `evaluate_many` returns every sample, including
zeros at the end. `jet_evaluate(p, x, order, context=...)` returns `order` Taylor
coefficients at `x`; these are derivatives divided by their factorials.

<!-- example: docs/examples/array_extension_numerical.mojo -->

The samples recover $(x+1)^2$. Dividing by $x+1$ produces the expected quotient
and zero remainder. For general floating-point data, a rounded quotient and
remainder can leave a residual. `divide_checked` retains bounds on that residual,
while `check_division` verifies them independently. For ball inputs, the evidence
also establishes a division enclosure for every represented input pair. The
leading divisor coefficient must exclude zero.

Use `divmod` when you only need the quotient/remainder pair. `affine_transform`
returns $p(ax+b)$, `taylor_shift` returns $p(x+b)$, and `definite_integral` computes
an endpoint difference. Sparse transforms and definite integration select an
`axis=` and retain the original power shape. `sparse.substitute` takes
`(axis, replacement)` pairs and replaces them simultaneously. Univariate sparse
`divmod`, `compose`, and `jet_evaluate` use a dense conversion and charge its full
extent before allocation.

For interpolation with balls, node enclosures must be separated so their
differences exclude zero. `hermite_interpolate` treats each row as a value followed
by consecutive derivatives; pass `derivatives=False` for Taylor coefficients.

## Numerical series and alternative bases

<!-- example: docs/examples/array_extension_representations.mojo -->

The series has five known coefficients, through $x^4$, with an unknown formal
tail $O(x^5)$. Its balls enclose coefficient values; they do not bound the analytic
remainder of $e^x$ at a chosen point. Numerical series provide the same arithmetic,
calculus, composition, logarithm, exponential, and reversion operations as exact
series. Construct with `series.asarray(..., order=..., context=...)`, then pass a
context to each computation. Different formal orders require an explicit common
lower `order=`; increasing a prepared series' known order raises `DOMAIN`.

`basis.asarray` interprets supplied coefficients in a named basis;
`basis.from_polynomial` converts monomial coefficients. Both accept `context=`
for Newton, Chebyshev, and Bernstein representations. Coefficients retain their
representation extent, including zeros. Bernstein subdivision preserves that
degree. Newton arithmetic keeps the left operand's nodes and raises if they
cannot represent the resulting degree.

See [formal series](array-series.md), [alternative bases](array-bases.md), and
[certified evaluation](array-certified.md) for the mathematical contracts.

## Exchange representation metadata

Series and basis arrays can exchange core values with arbitrary variable IDs.
Use `to_core(..., variable=...)` to choose an ID; `from_core` normalizes it to
array axis zero. Fraction conversions take `variables=`, one core ID per array
axis, for both numerator and denominator.

<!-- example: docs/examples/array_surface_interop.mojo -->

These conversions copy stored coefficients and Newton nodes without rounding
them. Numerical series and basis imports still take `context=`: for a series,
it supplies the format of implicit zero coefficients. It does not change
stored nonzero coefficients. Formal order, basis kind, coefficient extent,
and unused fraction axes survive the conversion. Mutating an exported core
value cannot change the prepared input.

Use `add(p, scalar, context=...)`, `subtract`, `multiply`, or `divide` for
[scalar arithmetic](arrays.md#work-with-scalar-constants). Sparse scalar division
keeps sparse storage, even for very large powers; division by a polynomial
continues to use the quotient/remainder operations described above.
