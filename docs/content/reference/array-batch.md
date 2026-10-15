# Batch maps, broadcasting, and reductions

`poly.batch` applies budget-aware functions to shaped collections of values.
The [tutorial](../tutorial/array-batching.md) explains shape rules, traversal
order, callback signatures, and empty inputs.

`reduce` accepts `axis=` as an integer, a tuple or finite iterable of integers,
`Optional[Int]`, or `None` for all axes. `keepdims=True` retains reduced axes
with size 1. Negative axes count from the end; duplicate and out-of-range axes
raise `DOMAIN`. Each output is a row-major left fold from its own `initial`.

Without either shape keyword, `reduce` returns one value. Supplying `axis` or
`keepdims` returns a batch, including a zero-dimensional batch when all axes
are removed. `Array.item()` extracts the only value from a one-item result.
An empty axis selection applies `function(initial, item)` independently to each
input; an empty reduced group returns `initial` without calling the function.

Save and restore these values with [poly.io](array-io.md).

<!-- api: src/poly/batch.mojo Array -->

<!-- api: src/poly/batch.mojo asarray -->

<!-- api: src/poly/batch.mojo scalar -->

<!-- api: src/poly/batch.mojo broadcast_to -->

<!-- api: src/poly/batch.mojo vmap -->

<!-- api: src/poly/batch.mojo outer -->

<!-- api: src/poly/batch.mojo reduce -->

<!-- api: src/poly/batch.mojo accumulate -->

`accumulate` accepts one integer axis, including negative axes, or
`Optional[Int]`. It starts each line from `initial` and returns the original
shape. `axis=None` and the call without `axis` use one flattened scan. Invalid
axes raise `DOMAIN`; Boolean, floating-point, tuple, and iterable axes are
rejected at compilation.
