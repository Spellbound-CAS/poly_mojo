# Batch polynomial calculations

Use `poly.batch` to apply a function to a shaped collection of polynomials,
points, or results. Each item remains a complete value. A batch of three
polynomials has shape `(3,)` even when their coefficient counts differ.

## Reuse a scalar function

<!-- example: docs/examples/array_batching.mojo -->

The two input batches have shapes `(2, 1)` and `(1, 3)`. The binary map returns
a `(2, 3)` grid: each row uses one polynomial from the first batch, and each
column uses one from the second. The function itself handles one pair and
passes the shared budget into `dense.add`.

`batch.asarray(items)` makes a one-dimensional batch from a raw iterable.
Passing a prepared batch without `shape=` shares its values and preserves its shape. Set `shape=` to interpret
a flat iterable in row-major order. The product of the dimensions must equal
the number of items. `batch.scalar(value)` makes a zero-dimensional batch,
which broadcasts to any shape. Polynomial coefficients are never interpreted
as extra batch axes.

## Broadcast compatible shapes

Shapes align from the last axis. Corresponding dimensions must match or one
must be 1; missing leading axes act as dimensions of size 1. For example,
`(4, 1, 2)` and `(3, 2)` produce `(4, 3, 2)`. A dimension of size 1 can also
broadcast to an empty dimension of size 0. Incompatible shapes raise `DOMAIN`.

`broadcast_to(values, shape)` creates a shared view. It records zero strides
for repeated axes instead of copying polynomial coefficients. Indexing and
iteration visit logical items in row-major order; `shape()` returns an
immutable view of the dimensions.

## Maps, pairwise calculations, and folds

| Function | Result |
|---|---|
| `vmap(function, values)` | One output per item, preserving shape |
| `vmap(function, a, b)` | Binary map over the broadcast shape |
| `outer(function, a, b)` | Every pair; shape is `a.shape + b.shape` |
| `reduce(function, values, initial)` | One row-major left fold, returning a value |
| `reduce(function, values, initial, axis=..., keepdims=...)` | One fold per group of batch items, returning a batch |
| `accumulate(function, values, initial)` | Every running fold result, in the original shape |

A fold always has an explicit initial value. Reducing an empty batch returns
that value; accumulating an empty batch returns an empty batch. No operation
reassociates a fold, so numerical rounding follows the stated traversal order.

## Reduce selected batch dimensions

Use `axis=1` to combine the columns within each row of a two-dimensional batch.
With `keepdims=True`, each reduced dimension stays in the result with size 1,
so the result can broadcast against the original batch in a later binary map.

<!-- example: docs/examples/array_batch_reductions.mojo -->

The two row sums are the polynomials $6$ and $15+x$. Their batch has shape
`(2, 1)`; omitting `keepdims=True` would give `(2,)`. The final reduction combines
all six polynomials into $21+x$. None of these operations reduces a polynomial's
coefficient array or changes its variables.

`axis` accepts an integer, a tuple or finite iterable of integers, or `None`:

| Argument | Meaning |
|---|---|
| `axis=0` | Reduce the first batch dimension |
| `axis=-1` | Reduce the last batch dimension |
| `axis=(0, 2)` | Reduce two dimensions together |
| `axis=None` | Reduce every batch dimension |
| `axis=()` | Keep every dimension; apply the callback to `initial` and each item |

You can pass `Optional[Int]` when an axis is chosen at runtime. Duplicate axes,
including aliases such as `0` and `-2` for a two-dimensional batch, raise
`DOMAIN`. Out-of-range axes also raise `DOMAIN`. A zero-dimensional batch has
no numbered axes; use `None` or an empty tuple for it.

The call without `axis` or `keepdims` keeps its original scalar return type.
Supplying either keyword selects a batch result, even for `axis=None` or
`keepdims=False`. This keeps the return type fixed when the options change at
runtime. Use `result.item()` to extract the only value from a result of any
shape containing exactly one item; empty or larger results raise `DOMAIN`.

Each group starts with its own copy of `initial`. If a reduced dimension is
empty, the group's result is that initial value. For example, reducing axis 1
of shape `(2, 0, 3)` produces shape `(2, 3)` containing six initial values.
If a dimension that remains in the result is empty, the result is empty too.

Output groups follow the retained axes in row-major order. Within a group,
items follow the reduced axes in their original row-major order; writing
`axis=(2, 0)` instead of `(0, 2)` does not change the fold. This order also holds
for broadcast views and callbacks whose result depends on evaluation order.

## Keep a running result along an axis

`batch.accumulate(function, values, initial, axis=...)` returns every running
fold in the original shape. Each line along the selected axis starts with its
own `initial`. Negative axes and `Optional[Int]` work as they do for reductions.
A scan accepts one axis; tuples and axis lists are not scan arguments.

<!-- example: docs/examples/array_workflow_accumulate.mojo -->

The first row contains `1, 3, 6`; the second starts again from zero and contains
`4, 9, 15`. `axis=0` would accumulate down the columns. Callbacks visit items in
row-major order, and broadcast views use their logical values. Empty batches
stay empty. Omitting `axis`, or passing `None`, scans the flattened batch while
preserving its original shape. A scalar batch accepts `None` and has no numbered
axes. Accumulation preserves every dimension, so it does not take `keepdims`.

## Use an array library for layout

Keep transpose, stack, and concatenate operations in your multidimensional
array library. The boundary with `poly.batch` is an iterable of complete items
and an explicit shape. For APN numerical batches, pass the shape when preparing
the items, then construct an APN batch from the result when you need its array
operations again.

<!-- example: docs/examples/array_batch_interop.mojo -->

APN transposes the original `(2, 3)` numerical array. The polynomial package
then folds each row of the `(3, 2)` view and returns the sums to an APN vector.
Conversion consumes logical values in row-major order; it materializes new
container storage and does not preserve view strides. Scalar value buffers may
still be shared by their value types. The `poly.batch` construction and reduction
share their usual budgets; APN's construction has its own resource contract.

This bridge applies to scalar types supported by APN batches. A whole polynomial
or a result with evidence is one `poly.batch` item, but is not an APN numerical
scalar. Supporting those objects in another array library requires that
library's object-storage interface.

## Pass work through the callback

Callbacks are named Mojo functions with a final `mut budget: Budget` argument
and `raises PolyError`. Unary and binary maps can return another value type,
including a polynomial, scalar, tuple, or result with evidence. The pinned
Mojo toolchain uses thin function values here; captured closures are not
accepted. To pass runtime data, include it in batch items or broadcast a
`batch.scalar` as the second argument.

Maps and folds run sequentially and eagerly. They share one work allowance
across shape preparation, allocation, and every callback. A callback must
charge its own computation, usually by forwarding the budget into a library
operation. Output storage is allocated under that allowance; failed calls do
not change input arrays. See the [batch reference](../reference/array-batch.md)
for signatures.
