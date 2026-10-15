# Sparse polynomial arrays

Use `poly.sparse` for `(power, coefficient)` terms. Operands share a power
shape: `Int` for one variable, or a tuple of `Int` powers for several variables.
Axis numbers start at zero, and output terms use increasing lexicographic order.

Every computational function accepts `max_steps=` instead of the explicit
shared budget shown below. See [array conventions](array-api.md) for inference
and errors, or [exact array algebra](../tutorial/array-exact.md) for worked examples.

## Construction and interoperability

<!-- api: src/poly/sparse.mojo Array -->

<!-- api: src/poly/sparse.mojo asarray -->

<!-- api: src/poly/sparse.mojo astype -->

<!-- api: src/poly/sparse.mojo to_core -->

<!-- api: src/poly/sparse.mojo from_core -->

<!-- api: src/poly/sparse.mojo format -->

<!-- api: src/poly/sparse.mojo to_kernel_bytes -->

<!-- api: src/poly/sparse.mojo from_kernel_bytes -->

## Arithmetic and calculus

<!-- api: src/poly/sparse.mojo add -->

<!-- api: src/poly/sparse.mojo subtract -->

<!-- api: src/poly/sparse.mojo negative -->

<!-- api: src/poly/sparse.mojo multiply -->

<!-- api: src/poly/sparse.mojo divide -->

<!-- api: src/poly/sparse.mojo power -->

<!-- api: src/poly/sparse.mojo degree -->

<!-- api: src/poly/sparse.mojo leading_term -->

<!-- api: src/poly/sparse.mojo variables -->

<!-- api: src/poly/sparse.mojo derivative -->

<!-- api: src/poly/sparse.mojo integral -->

<!-- api: src/poly/sparse.mojo definite_integral -->

<!-- api: src/poly/sparse.mojo compose -->

<!-- api: src/poly/sparse.mojo substitute -->

<!-- api: src/poly/sparse.mojo affine_transform -->

<!-- api: src/poly/sparse.mojo taylor_shift -->

## Exact algebra

Use `method="modular"` on `gcd` or `gcd_checked` to select sparse modular
reconstruction over `ZZ` or `QQ`. The default `"auto"` preserves core dispatch.

Exact `divmod`, `divide_many`, and their checked counterparts accept `order=`
and `axes=`. These select the reduction order; returned arrays still use
canonical increasing lexicographic storage. See the [worked
examples](../tutorial/array-checked.md#choose-the-order-for-sparse-division).

<!-- api: src/poly/sparse.mojo divmod -->

<!-- api: src/poly/sparse.mojo divide_exact -->

<!-- api: src/poly/sparse.mojo pseudo_divmod -->

<!-- api: src/poly/sparse.mojo divide_many -->

<!-- api: src/poly/sparse.mojo content -->

<!-- api: src/poly/sparse.mojo primitive_part -->

<!-- api: src/poly/sparse.mojo gcd -->

<!-- api: src/poly/sparse.mojo extended_gcd -->

<!-- api: src/poly/sparse.mojo resultant -->

<!-- api: src/poly/sparse.mojo subresultants -->

Factorization keeps the inferred domain and supports several variables over
`AA` and `QQbar`. Tuple powers retain unused axes. Univariate factors have degree
at most two over `AA` and one over `QQbar`; multivariate irreducibles can have
higher degree in either field. See
[algebraic factorization](../tutorial/algebraic-coefficients.md#factor-over-real-or-complex-algebraic-numbers)
for a complete example.

<!-- api: src/poly/sparse.mojo factor -->

<!-- api: src/poly/sparse.mojo squarefree -->

## Evaluation

<!-- api: src/poly/sparse.mojo evaluate -->

<!-- api: src/poly/sparse.mojo evaluate_many -->

## Roots and polynomial systems

See [root workflows](../tutorial/array-roots.md) and
[equations and inequalities](../tutorial/array-systems.md) for complete examples.

<!-- api: src/poly/sparse.mojo real_roots -->

<!-- api: src/poly/sparse.mojo complex_roots -->

<!-- api: src/poly/sparse.mojo root_between -->

<!-- api: src/poly/sparse.mojo count_roots -->

<!-- api: src/poly/sparse.mojo sturm_sequence -->

<!-- api: src/poly/sparse.mojo sign_on_interval -->

<!-- api: src/poly/sparse.mojo real_root_cover -->

<!-- api: src/poly/sparse.mojo complex_root_cover -->

<!-- api: src/poly/sparse.mojo real_root_cover_to -->

<!-- api: src/poly/sparse.mojo complex_root_cover_to -->

<!-- api: src/poly/sparse.mojo real_root_checker -->

<!-- api: src/poly/sparse.mojo complex_root_checker -->

<!-- api: src/poly/sparse.mojo Basis -->

<!-- api: src/poly/sparse.mojo groebner -->

<!-- api: src/poly/sparse.mojo reduce -->

<!-- api: src/poly/sparse.mojo contains -->

<!-- api: src/poly/sparse.mojo Membership -->

<!-- api: src/poly/sparse.mojo membership -->

<!-- api: src/poly/sparse.mojo check_membership -->

<!-- api: src/poly/sparse.mojo check_groebner -->

<!-- api: src/poly/sparse.mojo dumps -->

<!-- api: src/poly/sparse.mojo loads -->

<!-- api: src/poly/sparse.mojo dump -->

<!-- api: src/poly/sparse.mojo load -->

<!-- api: src/poly/sparse.mojo loads_numeric -->

<!-- api: src/poly/sparse.mojo degree_bound -->

<!-- api: src/poly/sparse.mojo certified_degree -->

<!-- api: src/poly/sparse.mojo specialize -->

<!-- api: src/poly/sparse.mojo evaluate_algebraic -->

<!-- api: src/poly/sparse.mojo enclose -->

<!-- api: src/poly/sparse.mojo adaptive_evaluate -->

<!-- api: src/poly/sparse.mojo jet_evaluate -->

<!-- api: src/poly/sparse.mojo roots_at -->

## Retained exact evidence

See [checked exact results](array-exact-results.md) for functional accessors and the precise claim each checker establishes.

<!-- api: src/poly/sparse.mojo gcd_checked -->

<!-- api: src/poly/sparse.mojo coprimality -->

<!-- api: src/poly/sparse.mojo resultant_checked -->

<!-- api: src/poly/sparse.mojo resultant_certified -->

<!-- api: src/poly/sparse.mojo divide_checked -->

<!-- api: src/poly/sparse.mojo divide_many_checked -->

<!-- api: src/poly/sparse.mojo factor_checked -->

<!-- api: src/poly/sparse.mojo squarefree_checked -->

## Incremental enclosure and numerical transport

<!-- api: src/poly/sparse.mojo enclose_stream -->

<!-- api: src/poly/sparse.mojo enclose_checker -->

<!-- api: src/poly/sparse.mojo load_numeric -->

<!-- api: src/poly/sparse.mojo complex_roots_to -->

<!-- api: src/poly/sparse.mojo complex_roots_checker -->
