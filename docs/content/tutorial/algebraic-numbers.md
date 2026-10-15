# Work with exact real algebraic numbers

For coefficient-array inputs, start with the [array workflow](array-roots.md).

An exact algebraic value lets you calculate with numbers such as $\sqrt{2}$
without replacing them by decimal approximations. The library represents the
value by a rational polynomial and enough root information to select one real
root of that polynomial.

## Select a root and calculate with it

<!-- example: docs/examples/algebraic_numbers.mojo -->

The interval $(1,2)$ selects the positive root of $x^2-2$. The negative root
satisfies the same polynomial, so the interval is essential to the value's
identity. Squaring the selected root gives the rational number 2 exactly.

Arithmetic operations return certificates whose `.value` is the result. The
checker verifies both the polynomial relation and the selection of the result
root. Checking only that an answer satisfies an annihilating polynomial would
leave its conjugate ambiguous.

Factories normalize rational results to linear defining polynomials and equal
rational endpoints. A certificate checker can also accept a valid nonminimal
defining polynomial; acceptance does not prove that a supplied polynomial is
irreducible.

## Compare, evaluate, and specialize

Use `algebraic_sign` for a certified comparison with zero. Algebraic arithmetic
supports exact addition, subtraction, multiplication, negation, and inversion.
Inversion requires a nonzero value; combine it with multiplication for division.

The specialized evaluation APIs substitute real algebraic values into rational
sparse polynomials. `algebraic_roots` accepts a list of real algebraic
coefficients and returns all distinct real roots. The coefficients may be
algebraically dependent: the final certificate checks the original polynomial,
not just a larger elimination polynomial used during the search.

When specializing a sparse polynomial before root finding, assign every
variable except the root variable. Specialization may lower the degree. An
identically zero specialization raises `DOMAIN`, because its root set is not
a finite list.

## Choose the right representation

`RealAlgebraic` is a mutable core root witness. The array API wraps witnesses in
immutable `algebraic.Real` coefficients; `algebraic.Complex` combines two real
components. Both implement `CoefficientField` for ordinary dense and sparse
arithmetic. See [algebraic coefficients](algebraic-coefficients.md) for construction,
promotion, and operation coverage. For approximate real or complex answers with
bounds, use [numerical root covers](numerical-roots.md).

See the [algebraic reference](../reference/algebraic.md) for operation signatures
and the [extension theory](../theory/extensions.md) for elimination and root
selection. [CAD](cad.md) builds on these exact values to describe cells in
several real variables.
