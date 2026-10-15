# Dense polynomial arrays

Use `poly.dense` for ascending coefficients.

Every computational function accepts `max_steps=` instead of the explicit
shared budget shown below. See [array conventions](array-api.md) for inference
and errors, or [exact array algebra](../tutorial/array-exact.md) for worked examples.

## Construction and interoperability

<!-- api: src/poly/dense.mojo Array -->

<!-- api: src/poly/dense.mojo asarray -->

<!-- api: src/poly/dense.mojo astype -->

<!-- api: src/poly/dense.mojo to_core -->

<!-- api: src/poly/dense.mojo from_core -->

<!-- api: src/poly/dense.mojo from_roots -->

<!-- api: src/poly/dense.mojo format -->

<!-- api: src/poly/dense.mojo to_kernel_bytes -->

<!-- api: src/poly/dense.mojo from_kernel_bytes -->

## Arithmetic and calculus

<!-- api: src/poly/dense.mojo add -->

<!-- api: src/poly/dense.mojo subtract -->

<!-- api: src/poly/dense.mojo negative -->

<!-- api: src/poly/dense.mojo multiply -->

<!-- api: src/poly/dense.mojo divide -->

<!-- api: src/poly/dense.mojo power -->

<!-- api: src/poly/dense.mojo degree -->

<!-- api: src/poly/dense.mojo derivative -->

<!-- api: src/poly/dense.mojo integral -->

<!-- api: src/poly/dense.mojo definite_integral -->

<!-- api: src/poly/dense.mojo compose -->

<!-- api: src/poly/dense.mojo affine_transform -->

<!-- api: src/poly/dense.mojo taylor_shift -->

<!-- api: src/poly/dense.mojo jet_evaluate -->

## Exact algebra

<!-- api: src/poly/dense.mojo divmod -->

<!-- api: src/poly/dense.mojo divide_exact -->

<!-- api: src/poly/dense.mojo pseudo_divmod -->

<!-- api: src/poly/dense.mojo divide_many -->

<!-- api: src/poly/dense.mojo content -->

<!-- api: src/poly/dense.mojo primitive_part -->

<!-- api: src/poly/dense.mojo gcd -->

<!-- api: src/poly/dense.mojo extended_gcd -->

<!-- api: src/poly/dense.mojo resultant -->

<!-- api: src/poly/dense.mojo subresultants -->

Factorization keeps the inferred domain. `AA` gives linear factors and
irreducible quadratics; `QQbar` gives linear factors. Both algebraic fields
require a univariate input; sparse tuple powers may include unused axes. See
[algebraic factorization](../tutorial/algebraic-coefficients.md#factor-over-real-or-complex-algebraic-numbers)
for a complete example.

<!-- api: src/poly/dense.mojo factor -->

<!-- api: src/poly/dense.mojo squarefree -->

## Evaluation and interpolation

<!-- api: src/poly/dense.mojo evaluate -->

<!-- api: src/poly/dense.mojo evaluate_many -->

<!-- api: src/poly/dense.mojo interpolate -->

<!-- api: src/poly/dense.mojo hermite_interpolate -->

<!-- api: src/poly/dense.mojo prepare_points -->

<!-- api: src/poly/dense.mojo PointPlan -->

## Roots and polynomial systems

See [root workflows](../tutorial/array-roots.md) and
[equations and inequalities](../tutorial/array-systems.md) for complete examples.

<!-- api: src/poly/dense.mojo real_roots -->

<!-- api: src/poly/dense.mojo complex_roots -->

<!-- api: src/poly/dense.mojo root_between -->

<!-- api: src/poly/dense.mojo count_roots -->

<!-- api: src/poly/dense.mojo sturm_sequence -->

<!-- api: src/poly/dense.mojo sign_on_interval -->

<!-- api: src/poly/dense.mojo real_root_cover -->

<!-- api: src/poly/dense.mojo complex_root_cover -->

<!-- api: src/poly/dense.mojo real_root_cover_to -->

<!-- api: src/poly/dense.mojo complex_root_cover_to -->

<!-- api: src/poly/dense.mojo real_root_checker -->

<!-- api: src/poly/dense.mojo complex_root_checker -->

<!-- api: src/poly/dense.mojo dumps -->

<!-- api: src/poly/dense.mojo loads -->

<!-- api: src/poly/dense.mojo dump -->

<!-- api: src/poly/dense.mojo load -->

<!-- api: src/poly/dense.mojo loads_numeric -->

<!-- api: src/poly/dense.mojo degree_bound -->

<!-- api: src/poly/dense.mojo certified_degree -->

<!-- api: src/poly/dense.mojo enclose -->

<!-- api: src/poly/dense.mojo adaptive_evaluate -->

<!-- api: src/poly/dense.mojo NumericDivision -->

<!-- api: src/poly/dense.mojo divide_checked -->

<!-- api: src/poly/dense.mojo check_division -->

## Retained exact evidence

See [checked exact results](array-exact-results.md) for functional accessors and the precise claim each checker establishes.

<!-- api: src/poly/dense.mojo gcd_checked -->

<!-- api: src/poly/dense.mojo coprimality -->

<!-- api: src/poly/dense.mojo resultant_checked -->

<!-- api: src/poly/dense.mojo resultant_certified -->

<!-- api: src/poly/dense.mojo divide_many_checked -->

<!-- api: src/poly/dense.mojo factor_checked -->

<!-- api: src/poly/dense.mojo squarefree_checked -->

<!-- api: src/poly/dense.mojo interpolate_checked -->

<!-- api: src/poly/dense.mojo hermite_interpolate_checked -->

## Incremental enclosure and numerical transport

<!-- api: src/poly/dense.mojo enclose_stream -->

<!-- api: src/poly/dense.mojo enclose_checker -->

<!-- api: src/poly/dense.mojo load_numeric -->

<!-- api: src/poly/dense.mojo complex_roots_to -->

<!-- api: src/poly/dense.mojo complex_roots_checker -->
