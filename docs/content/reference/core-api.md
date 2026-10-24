# Core API

Applications use `poly`: its functions take coefficient sequences, infer
domains, and return immutable arrays and checkable results. They are built on
the core package, `polynomials`, which provides the mutable containers,
coefficient rings, algorithms, and certificate checkers. Use the core API when
you extend the library, add an algorithm, or integrate with a system that
uses the [core compatibility contract](contract-p.md). Both packages share `Budget`, `PolyError`,
and the option types such as `NumericOptions`.

| Topic | Reference |
|---|---|
| Domains, contexts, containers, orders, arithmetic, errors, and binary codec | [Representations and arithmetic](p1.md) |
| Division, gcds, resultants, fractions, and identity certificates | [Division and elimination](p2.md) |
| Squarefree, finite-field, integer, rational, and Gaussian factorization | [Univariate factorization](p3.md) |
| Bases, reduction, and ideal membership | [Gröbner bases](p4.md) |
| Sturm sequences, root counts, isolation, and refinement | [Exact real roots](p5.md) |
| Calculus, interpolation, bases, series, rational utilities, and text transport | [Polynomial toolkit](toolkit.md) |
| Point and ball arithmetic, enclosures, root covers, and incremental witnesses | [Numerical polynomials](../spec/numerical.md) |
| Exact algebraic arithmetic and specialization | [Real algebraic numbers](algebraic.md) |
| Sparse modular gcds and multivariate factors | [Multivariate factorization](advanced-factorization.md) |
| Projection, cells, complete decomposition, and selective coverings | [CAD](cad.md) |

Exact operations are generic over a coefficient ring or field where supported.
Numerical operations take explicit arithmetic contexts. Most computational
calls accept a mutable `Budget` and raise `PolyError`; see
[errors and work limits](../guides/errors.md). The exact public surface is
described by the [compatibility contract](contract-p.md). Numerical operations
also follow the [precision and certificate guarantees](../spec/numerical.md).

## Exchange exact arrays with core values

<!-- example: docs/examples/array_exact_core.mojo -->

`to_core` gives you independent mutable `DensePoly` or `SparsePoly` storage.
`from_core` validates and copies it back into an immutable array. Changing
either side afterwards cannot change the other. Dense import normalizes its
single variable ID to zero. Sparse import accepts any valid core term order
and returns the array convention of increasing lexicographic order.

Sparse axes and core variable IDs need not coincide. Supply the same
`variables=` map in each direction, and set `from_core[ndim=2]` for two axes,
including axes unused by the polynomial. The map is in array-axis order;
unknown used IDs and duplicate map entries raise `DOMAIN`.

The example uses the core factorization result with its independent checker.
That checker establishes the factor product, not irreducibility. Array users
can retain and check the same evidence through the
[checked array operations](../tutorial/array-checked.md), including
`factor_checked`. Use the core bridge when you need mutable core values or an
existing core algorithm interface. The example shares one allowance across
construction, conversion, computation, and checking.

## Exchange numerical arrays

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

## Exchange series, bases, and rational functions

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

`rational.to_core(f)` returns a core `RationalFunction` for QQ. For QQi, GF(p),
AA, and QQbar, it returns `(numerator, denominator)` as public `SparsePoly`
values, since the core rational-function container is QQ-specific. Both forms
own independent storage. Numerical fractions export the pair as numerical
sparse polynomials. Import any of these forms with `from_core[ndim=n](value)`;
zero-axis constants are supported.

`rational.to_core(f, variables=ids)` assigns one core variable ID to each
array axis, including unused axes. Restore with
`rational.from_core[ndim](core, variables=ids)`. IDs may be signed and need
not increase; their positions define array-axis order. Duplicate IDs or
unmapped variables raise `DOMAIN`. The same mapping works for exact and
numerical fractions. Numerical copies preserve stored formats and coefficient
bounds without a new context. Imports still validate both polynomials and
require a denominator family that excludes the identically zero polynomial.

## Algebraic coefficients in core containers

`algebraic.Real` and `algebraic.Complex` implement `CoefficientField`.
`dense.to_core` therefore returns a `DensePoly` parameterized by the
corresponding type, and sparse conversion works the same way. The mutable
polynomial buffer is independent; coefficient copies share immutable witnesses.

`algebraic.asreal(core_witness)` validates and copies a mutable `RealAlgebraic`
witness. It accepts valid nonminimal defining polynomials and normalizes the
value for equality; later edits to the input witness cannot change the scalar.
See the [algebraic array reference](array-algebraic.md) for scalar helpers.
