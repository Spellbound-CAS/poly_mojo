# Solve polynomial sign problems with CAD

For coefficient-array inputs, start with the [array workflow](array-systems.md).

Cylindrical algebraic decomposition (CAD) divides real space into cells on
which your input polynomials have constant signs. It can describe a whole
solution set or find a witness for a conjunction of polynomial inequalities.
The result carries evidence for projection, root selection, and coverage.

CAD can be expensive as the number of variables and degrees grow. Choose the
smallest operation that answers your question.

| Task | Operation |
|---|---|
| Decompose the whole space into sign-invariant cells | `cad` |
| Locate a supplied point in its cell | `cad_cell` |
| Check a conjunction using a complete decomposition | `cad_covering` |
| Search selectively for a witness or a covering proof of failure | `cylindrical_covering` |

## Describe a circle and its interior

<!-- example: docs/examples/cad.mojo -->

The polynomial $x^2+y^2-1$ is zero on the unit circle, negative inside it,
and positive outside. A complete decomposition produces 13 cells in the
plane for the chosen variable order. The selective covering call finds an
interior witness using two stacks, without constructing the full decomposition.

The variable list is ordered **from base to top**. With `[0, 1]`, the algorithm
first decomposes the $x$ axis, then builds stacks in $y$ above those base cells.
Every input variable must be included. Changing the order can change both the
shape of the decomposition and the amount of work.

## Express the signs you want

A covering request supplies one sign mask per input polynomial. The bits are
1 for negative, 2 for zero, and 4 for positive; combine them to allow several
signs.

| Constraint on a polynomial $p$ | Mask |
|---|---|
| $p<0$ | `1` |
| $p=0$ | `2` |
| $p>0$ | `4` |
| $p\leq0$ | `3` |
| $p\geq0$ | `6` |
| $p\ne0$ | `5` |

Thus `[UInt8(1)]` requests the circle's interior. Several masks describe a
conjunction, one constraint for each input polynomial.

## Read cells geometrically

Each non-root cell records its parent, level, and position in the parent's
stack. Even positions are sectors between roots; odd positions are sections
on roots. A section's boundary is a root function over the base cell, not a
constant equal to its value at the sample point.

The sample point is useful for evaluating and naming a cell. Its signs alone
do not establish signs throughout the cell: the projection and lifting
certificate provides that guarantee. Use `check_cad` or the corresponding
covering checker before relying on the geometric claim.

## Projection choice and limits

The example selects Lazard projection with `lazard=True`. Its lifting handles
vanishing factors through certified derivative information. McCallum projection
can report `UNSUPPORTED` on nullification outside its supported cases; Lazard
is the alternative for those inputs.

A selective covering is not a partial CAD that can be treated as a full
partition. It proves the requested conjunction satisfiable through a witness,
or unsatisfiable through its covering evidence. Keep that narrower claim when
using or storing the result. See the [CAD reference](../reference/cad.md) for
cell fields and checker contracts, and [incremental transport](../guides/interchange.md)
when retaining the entire certificate is costly.
