# Root results and refinement

Root search starts in `poly.dense` or `poly.sparse`. Use `poly.roots` to inspect, verify, and refine the returned evidence.

Every computational operation has a shared `Budget` overload and a convenience
overload with `max_steps=`. Borrowed evidence is immutable; call `.copy()` on
it when you need independent mutable data for a core checker or serialization.

<!-- api: src/poly/roots.mojo ExactRoots -->

<!-- api: src/poly/_complex_roots.mojo ComplexRoots -->

<!-- api: src/poly/_complex_roots.mojo ComplexRootsCertificate -->

<!-- api: src/poly/roots.mojo check_complex -->

<!-- api: src/poly/roots.mojo RootCount -->

<!-- api: src/poly/roots.mojo AlgebraicCount -->

<!-- api: src/poly/roots.mojo check -->

<!-- api: src/poly/roots.mojo RealCover -->

<!-- api: src/poly/roots.mojo refine -->

<!-- api: src/poly/roots.mojo ComplexCover -->

<!-- api: src/poly/roots.mojo intervals -->

<!-- api: src/poly/roots.mojo rectangles -->

<!-- api: src/poly/roots.mojo NumericCover -->

<!-- api: src/poly/roots.mojo Multiplicities -->

<!-- api: src/poly/roots.mojo multiplicities -->

<!-- api: src/poly/roots.mojo Diagnostics -->

<!-- api: src/poly/roots.mojo diagnostics -->

<!-- api: src/poly/roots.mojo Polished -->

<!-- api: src/poly/roots.mojo polish -->

<!-- api: src/poly/roots.mojo SpecializedRoots -->

Interval queries accept real algebraic coefficients and endpoints. For `AA` or
`QQbar` coefficients, `roots.diagnostics` requires `context=BallContext(...)`;
`roots.polish` uses its explicit complex context. Both retain and check the
coefficient-conversion witnesses.

Use `poly.io` to archive complete root results, accepted covers, multiplicities,
diagnostics, and successful polishing evidence. An incomplete cover or an
unverified polishing result is rejected by the archive writer and reader.

Algebraic covers bind the numerical regions to their exact `AA` or `QQbar`
coefficients. Conversion uses `max_precision`; refinement encloses the original
coefficients again under the new options. Streaming checkers prepare and verify
the same coefficient family without retaining completed local proofs.

<!-- api: src/poly/_algebraic_covers.mojo AlgebraicCover -->

<!-- api: src/poly/_algebraic_covers.mojo AlgebraicRealRootChecker -->

<!-- api: src/poly/_algebraic_covers.mojo AlgebraicComplexRootChecker -->

<!-- api: src/poly/_complex_roots.mojo ExactComplexRootEvent -->

<!-- api: src/poly/_complex_roots.mojo ExactComplexRootReceipt -->

<!-- api: src/poly/_complex_roots.mojo ExactComplexRootSink -->

<!-- api: src/poly/_complex_roots.mojo ExactComplexRootChecker -->

## Prepared Sturm queries

Build a rational plan with `dense.sturm_sequence` or `sparse.sturm_sequence`.
Repeated counts share its chain, while `check` verifies the source binding,
quotient identities, and endpoint variations. See [reusable Sturm
chains](../tutorial/array-roots.md#reuse-a-sturm-chain) for a worked example.

<!-- api: src/poly/_sturm.mojo SturmPlan -->

<!-- api: src/poly/_sturm.mojo SturmCount -->

<!-- api: src/poly/roots.mojo sign_changes -->

<!-- api: src/poly/roots.mojo count_roots -->

<!-- api: src/poly/roots.mojo polynomials -->
