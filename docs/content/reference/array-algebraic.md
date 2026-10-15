# Exact algebraic values

Use `poly.algebraic` for exact real and complex algebraic arithmetic. Its scalar
types also serve as coefficients in the inferred `AA` and `QQbar` fields. Start
with [algebraic coefficients](../tutorial/algebraic-coefficients.md) for construction
and polynomial workflows; approximation stays explicit.

Every computational operation has a shared `Budget` overload and a convenience
overload with `max_steps=`. Borrowed evidence is immutable; call `.copy()` on
it when you need independent mutable data for a core checker or serialization.

Save and restore these values with [poly.io](array-io.md).

<!-- api: src/poly/_algebraic_value.mojo Real -->

<!-- api: src/poly/_algebraic_value.mojo Complex -->

<!-- api: src/poly/algebraic.mojo ascomplex -->

<!-- api: src/poly/algebraic.mojo conjugate -->

<!-- api: src/poly/algebraic.mojo asreal -->

<!-- api: src/poly/algebraic.mojo add -->

<!-- api: src/poly/algebraic.mojo subtract -->

<!-- api: src/poly/algebraic.mojo multiply -->

<!-- api: src/poly/algebraic.mojo negative -->

<!-- api: src/poly/algebraic.mojo reciprocal -->

<!-- api: src/poly/algebraic.mojo divide -->

<!-- api: src/poly/algebraic.mojo compare -->

<!-- api: src/poly/algebraic.mojo check -->

<!-- api: src/poly/algebraic.mojo bounds -->

<!-- api: src/poly/algebraic.mojo refine -->

<!-- api: src/poly/algebraic.mojo approximate -->

<!-- api: src/poly/algebraic.mojo evaluate -->

<!-- api: src/poly/algebraic.mojo Specialization -->

<!-- api: src/poly/algebraic.mojo coefficients -->

## Retained ordering evidence

`compare_checked` and `sign_checked` retain their operands and an ordering certificate. Read `relation()` for the claimed sign and use `algebraic.check` to verify it.

<!-- api: src/poly/algebraic.mojo Comparison -->

<!-- api: src/poly/algebraic.mojo compare_checked -->

<!-- api: src/poly/algebraic.mojo sign_checked -->

For `AA` coefficients and real assignments, `sparse.specialize` returns an
iterable `AlgebraicSpecialization`. Its checker binds the original polynomial
and assignments to each returned coefficient. `sparse.real_roots` accepts this
result; `sparse.roots_at` retains it together with complete root evidence.

<!-- api: src/poly/_specialization.mojo AlgebraicSpecialization -->
