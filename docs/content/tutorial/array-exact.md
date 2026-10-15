# Exact algebra with polynomial arrays

Compute GCDs, factor polynomials, recover exact coefficients from samples, and
change variables through `poly.dense` and `poly.sparse`. Start with coefficient
or term sequences; the functions infer domains and return arrays you can use
in the next calculation. If these inputs are new to you, begin with
[polynomial arrays](arrays.md).

## Inspect sparse terms

Use `sparse.leading_term(p)` to select the greatest term in lexicographic order,
with axis zero first. Select `order="grlex"` or `order="grevlex"` for graded orders,
and pass `axes=` as a complete permutation when another axis should take
priority. The query leaves the array's canonical storage order unchanged.

<!-- example: docs/examples/array_access_inspection.mojo -->

`sparse.variables(p)` returns sorted axis numbers that actually occur. The
example uses a three-axis power shape, but only axes zero and one occur in its
terms. Constants and zero return an empty list; zero has no leading term and
raises `DOMAIN` for that query. Exact raw terms are combined before inspection,
so canceled terms do not contribute to the result.

These queries also accept prepared numerical arrays without a context. For
balls they describe stored terms: a coefficient containing both zero and
nonzero values still contributes a term. Its presence does not prove that every
polynomial in the family has the same leading term.

## Find and remove common factors

<!-- example: docs/examples/array_exact_algebra.mojo -->

The integer GCD includes coefficient content, so the example returns
$2(x-1)$. Over a field, `gcd` returns a monic polynomial. `extended_gcd` also
returns multipliers: `(g, u, v)` satisfies $up+vq=g$. It promotes integer
inputs to `QQ`, where the nonzero GCD can always be made monic.

Choose division according to the result you need:

| Operation | Result | Integer inputs |
|---|---|---|
| `divmod(p, q)` | Quotient and remainder | Promote to `QQ` |
| `divide_exact(p, q)` | An optional exact quotient | Stay in `ZZ` |
| `pseudo_divmod(p, q)` | Scale, quotient, and remainder | Stay in `ZZ` |
| `divide_many(p, divisors)` | Quotients in divisor order, then a remainder | Promote to `QQ` |

For example, `divide_exact([1, 1], [2])` returns `None`: $(1+x)/2$ has no
integer coefficient representation. Use a rational coefficient when you want
the quotient over `QQ`. Invalid inputs and work exhaustion still raise errors;
they are not negative divisibility answers. Pseudo-division instead multiplies the dividend by
$\operatorname{lc}(q)^{\max(\deg p-\deg q+1,0)}$. The example gives
$4(x^2+1)=(2x-1)(2x+1)+5$.

Both namespaces provide `content` and `primitive_part`. Dense content is a
scalar; `sparse.content(p, axis=0)` treats the other axes as coefficient
variables and returns a polynomial in them. Multiplying content by primitive
part reconstructs the input.

`divide_many` accepts a finite collection of polynomials of one Mojo input
type. A list of prepared arrays works well; for nested coefficient lists,
annotate the collection as `List[List[Int]]` or the corresponding scalar type.
Sparse reduction defaults to lexicographic order, with axis zero first. Use
`order=` and `axes=` to choose another [reduction order](array-checked.md#choose-the-order-for-sparse-division). The
remainder can depend on both monomial and divisor order; use the [prepared Gröbner API](array-systems.md) when
you need reduction by a Gröbner basis.

## Factor or group repeated roots

`factor` returns `(scalar, factors)`, with each factor stored alongside its
positive multiplicity. The product in the example checks reconstruction.
Integer factors are primitive with positive leading coefficient; rational
factors are monic. Zero and constants return an empty factor list.

Use `squarefree` when you need multiplicity groups without splitting every
irreducible factor. Its monic blocks are ordered by multiplicity, and integer
inputs promote to `QQ`. A block can contain several irreducible factors.

| Domain | `factor` | `squarefree` |
|---|---|---|
| `ZZ` | Primitive integer factors | Monic blocks over `QQ` |
| `QQ`, `QQi` | Monic factors in the inferred field | Monic blocks in that field |
| Prime field | Monic factors with multiplicities | Monic blocks, including inseparable inputs |
| `AA` | Monic real algebraic factors; degree one or two for univariate inputs | Monic blocks in `AA` |
| `QQbar` | Monic complex algebraic factors; linear for univariate inputs | Monic blocks in `QQbar` |

Both operations accept dense univariate inputs and sparse multivariate inputs
over every listed domain. The
coefficient field matters: $x^2+1$ stays irreducible over `QQ`, but splits into
$(x-i)(x+i)$ over `QQi`. Choose exact complex coefficients to request that field.

<!-- example: docs/examples/array_factor_fields.mojo -->

Use `algebraic.asreal(1)` or `algebraic.ascomplex(1)` in the coefficient sequence
to choose `AA` or `QQbar`, even when the polynomial has rational coefficients.
The [algebraic factorization example](algebraic-coefficients.md#factor-over-real-or-complex-algebraic-numbers)
shows how the field changes the factors and how `exact.check` verifies their
irreducibility in these two domains.

Finite-field factorization extracts the leading scalar and handles repeated
factors automatically. Its `seed=0` default fixes the splitting sequence; other
seeds change the work, while the sorted factors remain the same. The example
also shows why a zero derivative need not mean a constant polynomial: in
characteristic three, the derivative of $x^3+1$ vanishes.

Sparse squarefree decomposition uses partial derivatives directly. General
QQi, prime-field, `AA` and `QQbar` sparse factorization uses a dense substitution image,
whose size grows with the product of the per-variable degree ranges, followed
by a potentially exponential subset search. Set `max_steps` for large inputs;
an exhausted budget raises `LIMIT` before the corresponding work or allocation.
See [the algorithm notes](../theory/factorization.md#sparse-field-factorization)
for the cost and correctness argument.

Sparse `extended_gcd` and `pseudo_divmod` still require a one-axis power shape
and use dense conversion internally.

## Recover coefficients from exact samples

<!-- example: docs/examples/array_exact_samples.mojo -->

`dense.interpolate(points, values)` performs a single interpolation.
`prepare_points` builds a product tree you can reuse for several polynomials;
`interpolation=True` also caches interpolation weights. Samples keep their
order and all zeros. Distinct nodes and matching sample counts determine a
unique polynomial of degree below the node count.

The plan retains the field inferred from its nodes. Integer nodes select
`QQ`. To work with exact complex coefficients, prepare the plan from nodes
containing an `ExactComplex`, even if its imaginary part is zero. Modular
nodes must share a prime. Copies of the plan share its immutable tree.
Repeated nodes are allowed for evaluation when `interpolation=False`, but
interpolation rejects them.

Hermite interpolation adds derivative constraints. Each inner sequence starts
with the value, followed by consecutive derivatives. In the example, the two
jets describe $(x+1)^2$ at zero and one. Set `derivatives=False` to supply
Taylor coefficients instead. `jet_evaluate` returns those coefficients
$p^{(j)}(a)/j!$, including requested trailing zeros. It computes them without
factorial division, which also makes it useful in finite fields.

## Change variables and eliminate an axis

The sparse example swaps $x$ and $y$ simultaneously. Replacements keep the
original variables, so replacing $x$ with $y$ and $y$ with $x$ turns $x^2y+y$
into $xy^2+x$. To replace only one axis, use
`substitute(p, [replacement], axes=[0])`.

`affine_transform(p, scale, shift)` computes $p(\text{scale}\,x+\text{shift})$;
`taylor_shift(p, point)` computes $p(x+\text{point})$. Sparse versions take
`axis=`. `definite_integral` accepts exact scalar endpoints and promotes
integer inputs to `QQ`; sparse integration leaves exponent zero on the
integrated axis and preserves the full power shape.

For elimination, `resultant(p, q)` returns a scalar in the dense namespace.
`sparse.resultant(p, q, axis=0)` returns a polynomial in the remaining axes,
retaining the input shape. The sign follows the Sylvester determinant
convention. Two nonzero constants in the selected variable raise `DOMAIN`.
`subresultants` returns the indexed family $S_0,\ldots,S_{m-1}$, where $m$ is
the smaller input degree in that variable; zero entries retain their indices.

## Continue with checked results

You can retain and check exact algebra evidence directly through the
[checked array operations](array-checked.md), including `factor_checked`.
Roots, Gröbner bases, and CAD also have direct [array workflows](array-systems.md).

Every operation above accepts `max_steps=` or an explicit shared `Budget`.
See the [dense](../reference/array-dense.md) and
[sparse](../reference/array-sparse.md) references for signatures and domain
requirements.
