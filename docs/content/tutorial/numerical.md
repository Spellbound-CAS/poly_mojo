# Evaluate with arbitrary precision and bounds

Use numerical point arithmetic when you need a value at a chosen precision.
Use enclosure evaluation when you need a bound on the answer. Both workflows
support real and complex values, but they make different promises.

This example evaluates $x^2-2$ at $x=1.5$ using 128-bit arithmetic, then
bounds the values of the exact polynomial over a small interval around 1.5.

<!-- example: docs/examples/tutorial_numerical.mojo -->

## Read the point calculation

The example constructs a `NumericDensePoly[FloatDomain]` directly from point
coefficients. To convert an existing exact polynomial, use `numeric_from_dense`.
`ArithmeticContext` controls the precision and rounding used by
`FloatDomain`. A stored floating-point coefficient is a point value: its
original exact value, if different, is no longer part of the polynomial.

The example's calculation is exactly representable and gives $0.25$.
In general, each arithmetic step may round. Choosing a precision does not
promise a correctly rounded result for the entire polynomial expression.
For complex points, use `ComplexDomain` with `ComplexContext`.

## Read the enclosure calculation

`enclose_dense` keeps the original exact coefficients and evaluates over the
supplied ball. Its result includes an enclosure and a trace that an independent
checker can verify. `check_enclose_dense` reconstructs the computation from
the original polynomial and input region, then checks rational bounds for the
steps. It does not trust a midpoint or rerun the floating-point producer.

An accepted enclosure covers **every point in the input region**. It need not
be the narrowest possible bound. Dependency between repeated occurrences of
$x$ can widen interval evaluation even when the arithmetic is accurate.

For uncertain coefficients, use `NumericDensePoly[BallDomain]` or the sparse
counterpart and the numerical enclosure entry points. Their bounds cover the
represented coefficient family as well as the evaluation region. Complex
balls represent rectangles and use `BallContext` too.

## Ask for an accuracy target

Adaptive evaluation increases working precision up to `max_precision` in
`NumericOptions`. It tests the enclosure against the requested absolute and
relative width targets. Set `absolute_bits` for a useful target near zero;
a relative-only target cannot accept a positive-width interval around zero.

An enclosure can be valid without meeting the requested accuracy. Check the
enclosure certificate, inspect the outcome, and use `check_numeric_accuracy`
when your application relies on the width claim. At a precision ceiling, the
last valid enclosure can still be useful.

Increasing precision reduces arithmetic error. To reduce uncertainty already
present in coefficients or arguments, you need better input data or smaller
input regions. See [numerical root covers](numerical-roots.md) for the same
distinction in root finding and the [numerical reference](../spec/numerical.md)
for conversion and accuracy contracts.
