# Numerical coefficient arrays

Pass `context=` to `dense.asarray` or `sparse.asarray` to select numerical
arithmetic. The [tutorial](../tutorial/array-numeric.md) explains conversion,
rounding, enclosure, and degree uncertainty. Computational calls accept a
shared `Budget` or `max_steps=`; contexts remain explicit in arithmetic.

The dense and sparse references include numerical overloads beside their
exact counterparts. Numerical array storage is exposed as
`dense.NumericArray` and `sparse.NumericArray`.

<!-- api: src/poly/_numeric.mojo DenseArray -->

<!-- api: src/poly/_numeric.mojo SparseArray -->

For method signatures, see [dense arrays](array-dense.md) and
[sparse arrays](array-sparse.md). Root covers, diagnostics, and polishing
are described in the [root reference](array-roots.md).

<!-- api: src/poly/_numeric.mojo Values -->

For checked fixed and adaptive evaluation, see [evaluation results](array-evaluation.md).
