# Polynomial array API

Import `poly` for coefficient-first polynomial calculations. Start with
[dense](array-dense.md) and [sparse](array-sparse.md) for ordinary polynomials.
The sections below cover scalar helpers, series and basis representations,
rational functions, roots, and polynomial systems.
See the [array tutorial](../tutorial/arrays.md) and
[exact algebra walkthrough](../tutorial/array-exact.md) for complete examples.

Computational functions have a `mut budget: Budget` overload and a
convenience overload with `max_steps: UInt64 = UInt64.MAX`. A shared budget is
the last positional argument; `modulus` and `axis` are keyword arguments.
Inputs are borrowed unless an owning iterator is explicitly moved.

The compiler infers generic parameters from input types. Internal kind values
in rendered signatures mean `0=ZZ`, `1=QQ`, `2=QQi`, `3=GF(p)`, `4=AA`, and `5=QQbar`; callers do
not need to supply them. Exact arrays support iteration, indexing, `len`,
`domain()`, `modulus()`, and `to_list()`. Scalar results are APN values or `ModInt`.
Unsupported scalar types fail compilation. Invalid powers, moduli, axes, or
zero divisors raise `PolyError` with kind `DOMAIN`; budget exhaustion raises
`LIMIT`, and unrepresentable allocation sizes raise `CAPACITY`. Algorithms
without support for an inferred domain raise `UNSUPPORTED`; the operation
references list those boundaries.

`dense.to_core` and `sparse.to_core` return independent mutable core storage.
Their `from_core` counterparts validate and copy it back. Use
`sparse.from_core[ndim=2](p, variables=[7, 2])` to map core IDs 7 and 2 onto
array axes zero and one. This explicit dimension is needed even for zero or
constant input; coefficient domains are still inferred.

<!-- api: src/poly/_common.mojo Values -->

<!-- api: src/poly/_common.mojo coefficients -->

<!-- api: src/poly/_common.mojo modint -->

<!-- api: src/poly/_conversion.mojo to_dense -->

<!-- api: src/poly/_conversion.mojo to_sparse -->

## Root and system workflows

- [Root results](array-roots.md): exact real roots, counts, numerical covers,
  independent checks, and functional refinement.
- [Algebraic values](array-algebraic.md): exact arithmetic, comparisons, bounds,
  and explicit rational approximation.
- [Sparse bases](array-sparse.md#roots-and-polynomial-systems): reusable Gröbner
  bases, normal forms, membership, and identity witnesses.
- [CAD](array-cad.md): readable sign constraints, exact samples, full
  decompositions, and streamed checking.

## Series, bases, and rational functions

- [Formal series](array-series.md): exact coefficients with an explicit unknown-tail order.
- [Alternative bases](array-basis.md): Newton, Chebyshev, and Bernstein coefficients, conversion, arithmetic, and subdivision.
- [Rational functions](array-rational.md): normalized fractions, calculus, partial fractions, and Laurent expansions.

For rounded and enclosing coefficients, use [numerical arrays](array-numeric.md). [Batch maps](array-batch.md) apply reusable functions across whole polynomials. The [serialization guide](../tutorial/array-serialization.md) covers complete values and existing coefficient records; [poly.io](array-io.md) provides string and streaming entry points.

[Certified evaluation](array-evaluation.md) retains independent enclosure evidence and explicit adaptive accuracy outcomes.

[Checked exact results](array-exact-results.md) retain algebra evidence for functional access and independent checking. [Specialized CAD workflows](../tutorial/array-cad-workflows.md) expose projection, selected cells, Lazard evaluation, and full covering.


Exact [queries and algorithm choices](../tutorial/array-checked.md#exact-queries-and-algorithm-choices)
include closed-interval signs, determinant-certified resultants, and sparse modular
GCD selection. Sparse division also accepts monomial order and axis priority.
