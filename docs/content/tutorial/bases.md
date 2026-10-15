# Newton, Chebyshev, and Bernstein bases

For coefficient iterables and functional calls, start with the [array interface](array-bases.md).

A change of basis represents the same polynomial with different coefficients.
`BasisPolynomial` supports Newton, first-kind Chebyshev, and Bernstein bases
over an exact coefficient field. `NumericBasisPolynomial` provides the
corresponding numerical operations with an explicit context.

## Represent a quadratic in two ways

The polynomial $x^2$ is the degree-two Bernstein basis element $B_{2,2}$,
so its Bernstein coefficients are $0,0,1$. In the Chebyshev basis,

$$
x^2=\tfrac12 T_0(x)+\tfrac12 T_2(x),\qquad T_2(x)=2x^2-1.
$$

<!-- example: docs/examples/tutorial_bases.mojo -->

Both conversions expand back to the original dense polynomial. The example
also evaluates the Bernstein representation at $1/2$, obtaining exactly
$1/4$ without first converting back to monomial coefficients.

## Choose the representation

| Basis | Data you supply | Evaluation method |
|---|---|---|
| Newton | Nodes defining successive factors $(x-a_0),(x-a_1),\ldots$ | Nested multiplication |
| Chebyshev | A polynomial to express in the first-kind basis | Clenshaw recurrence |
| Bernstein | A retained degree, at least the polynomial degree | de Casteljau evaluation |

Newton nodes may repeat; there must be enough nodes for the polynomial's
degree. A Newton basis conversion is a representation change, distinct from
reconstructing a polynomial from sample values by interpolation.

Bernstein polynomials use the interval $[0,1]$. The representation retains
its chosen degree, even when the polynomial has a smaller degree. Increasing
that degree changes the coefficient list while preserving the polynomial.

## Subdivide a Bernstein polynomial

`subdivide(a, work)` returns representations of $p(ax)$ and
$p(a+(1-a)x)$ at the retained degree. For $0<a<1$, these parameterize the
left and right portions of the original interval. Both returned polynomials
again use a local parameter on $[0,1]$.

Exact field operations need the conversion factors to be invertible.
Chebyshev or Bernstein conversion can therefore reject some
finite-characteristic inputs. Promote integer polynomials to rationals when
you need a general exact field calculation.

See the [basis reference](../reference/toolkit.md#alternative-bases) for the
constructors and numerical variants.
