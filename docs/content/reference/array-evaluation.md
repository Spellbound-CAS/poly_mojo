# Certified array evaluation

`dense.enclose` and `sparse.enclose` produce a value and an independently checkable
bound trace. Algebraic coefficients and points additionally retain checked
root-refinement and containment witnesses in the complete result. The borrowed
`evidence()` is the arithmetic trace only. Their `adaptive_evaluate` counterparts additionally try to meet
explicit accuracy options. The [tutorial](../tutorial/array-certified.md) explains
why valid evidence and successful accuracy are separate questions.

<!-- api: src/poly/evaluation.mojo Enclosure -->

<!-- api: src/poly/evaluation.mojo check -->

## Incremental witnesses

`dense.enclose_stream` and `sparse.enclose_stream` accept the same
fixed-precision inputs as `enclose`. Construct `enclose_checker` with the same
inputs and context, forward every event, then check the final receipt. See the
[streaming tutorial](../tutorial/array-certified.md#check-one-witness-step-at-a-time)
for a complete example.

Algebraic coefficients or points select `AlgebraicEnclosureStream`. Its producer
refines one input per event; its checker verifies that input's root identity,
interval containment, and outward numerical bounds. Coefficient events precede
point events, then arithmetic events. The receipt remains a core
`StreamedEnclosure`, with `steps` counting the complete facade event sequence.
Pass it to the matching facade checker.

Without algebraic inputs, the return types remain the core `EnclosureStream`
and `EnclosureStreamChecker`, and events remain `NumericBox` values. Type
inference lets the same forwarding loop work with either family.

<!-- api: src/poly/evaluation.mojo AlgebraicEnclosureStream -->

<!-- api: src/poly/evaluation.mojo EnclosureStep -->

<!-- api: src/poly/evaluation.mojo RefinedInput -->
