# Keep and check exact algebra evidence

A coefficient array gives you an answer. A checked result also keeps the inputs
and evidence needed to verify a stated relation. Use the `_checked` operations
when a later calculation needs to validate the result independently.

<!-- example: docs/examples/array_parity_checked.mojo -->

The first calculation establishes

$$
1+x^2=(1+2x)\left(-\frac14+\frac12x\right)+\frac54.
$$

`exact.check` verifies the identity and the dense remainder's degree bound.
The GCD check additionally establishes maximality through Bézout or coprimality
evidence. `exact.value` returns the ordinary array, so subsequent calculations
use the same functional API.

## Choose the evidence you need

Both namespaces offer checked GCD, resultant, division, ordered list division,
factorization, and squarefree decomposition. Dense interpolation and
Hermite interpolation can retain their nodes and all value or derivative
constraints. `coprimality(a, b)` returns optional evidence for coprime operands.

A check proves the relation documented for that operation. Factor checks prove
recomposition; for `AA` and `QQbar`, they also verify irreducibility, monicity,
and exact multiplicities. Other factorization domains keep the product-only check.
`resultant_checked` proves an ideal-membership identity. Choose `resultant_certified`
when you need verification of the complete Sylvester determinant over any
supported exact field.
The [checked-results reference](../reference/array-exact-results.md) lists every
claim so you can decide whether it meets your application's needs.

Share one `Budget` across preparation, production, output extraction, and
checking when the entire workflow needs a work limit. Result copies share
storage, while extracted polynomial arrays and `to_core` results own independent
mutable data.

## Exact queries and algorithm choices

Use `sign_on_interval(p, low, high)` to establish a strict sign throughout a
closed interval. It returns `1` or `-1` when the interval contains no root, and
`None` when a root occurs anywhere, including either endpoint. A repeated root
also prevents a strict sign, even if the polynomial never changes sign. Equal
endpoints are allowed; reversed endpoints raise `DOMAIN`. Inputs and endpoints
can be integers, rationals, or real algebraic values. Sparse input needs a
one-variable power shape.
This query returns an exact answer without retaining a certificate.

<!-- example: docs/examples/array_exact_choices.mojo -->

`resultant_certified` returns a coefficient array and a complete determinant
transcript. Integer inputs promote to `QQ`; `QQi`, `GF(p)`, `AA`, and `QQbar`
retain their fields. Sparse inputs accept `axis=` to select the
variable to eliminate, and the output retains the original power shape. As with
other exact results, use `exact.value`, `exact.check`, and an immutable
`evidence()` borrow. Keeping the elimination transcript can require substantial
memory, so share a budget when working with larger degrees.

Choose `method="modular"` on `sparse.gcd` or `sparse.gcd_checked` to use sparse
modular reconstruction over integers or rationals. The algorithm discovers and
interpolates the occurring monomials, avoiding a dense box of parameter degrees.
Both methods compute the same normalized GCD. The checked form also retains
maximality evidence. `method="auto"` keeps the existing core dispatch; it does
not automatically choose this dedicated sparse algorithm. Other method names
raise `DOMAIN`; requesting modular reconstruction over `QQi` or a prime field
raises `UNSUPPORTED`.

## Choose the order for sparse division

`divmod`, `divide_many`, `divide_checked`, and `divide_many_checked` accept
`order="lex"`, `"grlex"`, or `"grevlex"`, plus an optional `axes=` permutation in
descending priority. For example, `axes=[1, 0]` makes the second variable more
significant. Every axis must appear exactly once; an empty list uses tuple order.
The default remains lexicographic order with axis zero first.

In the example, lexicographic division reduces $x$ by the leading $x$ of
$x+y^2$, leaving $-y^2$. Graded lexicographic order makes $y^2$ the leading
monomial, which cannot divide $x$, so the quotient is zero and the remainder
is $x$. Reversing the axis priority under lexicographic order gives the same
remainder in this example.

The chosen order affects the calculation, while output terms always return in
canonical increasing lexicographic order. You can pass them straight into other
array operations. Checked results retain the computational order in their
evidence. Their checker establishes the recomposition identity; it does not
independently certify reducedness or a unique remainder. Divisor list order also
matters. Use a Gröbner basis when you need ideal normal forms.
