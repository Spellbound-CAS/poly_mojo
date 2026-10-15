# Evaluation at many points and interpolation

Interpolation constructs a polynomial from prescribed values. Over exact
rational coefficients, both the input data and the resulting polynomial can
remain exact. A reusable plan is useful when several calculations share the
same sample points.

For coefficient and term inputs through `poly`, see
[exact array algebra](array-exact.md). This walkthrough uses the core types
and their explicit contracts.

## Reconstruct a quadratic

The values $1,4,9$ at $0,1,2$ determine a unique polynomial of degree below
three: $(x+1)^2$.

<!-- example: docs/examples/tutorial_interpolation.mojo -->

The example builds a plan, interpolates the supplied values, and checks the
result against the known quadratic. `check_interpolation` independently
checks the values at the nodes. The final evaluation reuses the same plan.

## Choose a plan

`MultipointPlan.prepare` builds a product tree for repeated evaluation. It
also supports interpolation, computing interpolation weights when needed.
`prepare_interpolation` additionally retains the weights, which is useful
when many value lists share the same distinct nodes.

Both plans are immutable. Preparing a plan and using it consume work from
the supplied budget. A plan retains its nodes and product tree; its memory
cost is part of the choice to reuse it.

Evaluation can repeat a node. Interpolation requires distinct nodes and one
value for each node. Repeated-node data need the Hermite interface, where
derivatives supply the additional conditions. A one-shot `interpolate` call
is also available when you do not need to keep a plan.

## Include derivative information

`hermite_interpolate` accepts values and derivative constraints at the
specified nodes. By default, its data are raw derivatives. With
`derivatives=False`, they are Taylor/Hasse coefficients instead: the
coefficient of $(x-a)^k$ at a node $a$.

Those conventions differ by a factor of $k!$ in characteristic zero. The
Taylor/Hasse form avoids assuming that these factorials are invertible,
which is useful in finite characteristic. Use `check_hermite_interpolation`
to verify the supplied conditions under the same convention.

## Interpolate uncertain data

`numeric_interpolate` and `numeric_hermite_interpolate` work with numerical
domains and explicit contexts. Ball interpolation requires separated node
intervals so its denominators can be certified nonzero. Raising the working
precision cannot separate overlapping input intervals.

For exact sample data that will later be evaluated numerically, interpolate
exactly first and then choose an explicit numerical conversion. See the
[interpolation reference](../reference/toolkit.md#interpolation) and
[domain guide](../guides/domains.md).
