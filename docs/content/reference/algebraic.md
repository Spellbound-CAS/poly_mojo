# Exact real algebraic numbers

For the array workflow, see [algebraic coefficients](../tutorial/algebraic-coefficients.md).
Values are exact real roots of monic rational polynomials with a unique-root
Sturm witness. Factories choose irreducible definitions and represent rationals
by a linear polynomial and equal endpoints. Checkers also accept nonminimal
definitions; they do not assert irreducibility.

Arithmetic returns certificates whose `value` is the answer. Binary operations
retain a resultant ideal identity, a nonzero annihilator, operand refinements,
and a root-count witness selecting the intended conjugate. Unary operations
verify the direct polynomial transformation. Independent checkers use those
identities and interval bounds without factoring, computing gcds or isolating
roots. Negation and inversion pass `None` as the checker's second operand.
`algebraic_sign` is a certified comparison against `real_rational(0)`.

`algebraic_roots` accepts coefficients in ascending degree and returns every
distinct real root. It eliminates independent coefficient variables to obtain
a nonzero rational annihilator, accounts for **all** its real roots, and
filters each candidate by certified evaluation of the original polynomial.
Algebraically dependent coefficients are valid. Constants have no roots;
the zero polynomial raises `DOMAIN` because it has no finite root list.

Sparse evaluation and partial specialization retain a fixed arithmetic
schedule. Specialization keeps even coefficients which become zero so its
checker can verify that none was omitted. `algebraic_roots_at` specializes all
variables except the designated root variable and checks both specialization
and completeness. A vanished leading coefficient is allowed; an identically
zero specialization raises `DOMAIN`.

All operations share the supplied budget and leave borrowed inputs unchanged
on failure. See the [extension mathematics](../theory/extensions.md).

<!-- example: docs/examples/algebraic_numbers.mojo -->

<!-- api: src/polynomials/algebraic_arithmetic.mojo AlgebraicArithmeticCertificate -->

<!-- api: src/polynomials/algebraic_types.mojo AlgebraicComparisonCertificate -->

<!-- api: src/polynomials/algebraic_evaluation.mojo AlgebraicEvaluationCertificate -->

<!-- api: src/polynomials/algebraic_specialization.mojo AlgebraicPolynomialEvaluation -->

<!-- api: src/polynomials/algebraic_specialization.mojo AlgebraicRootsAtCertificate -->

<!-- api: src/polynomials/algebraic_roots.mojo AlgebraicRootsCertificate -->

<!-- api: src/polynomials/algebraic_specialization.mojo AlgebraicSpecializationCertificate -->

<!-- api: src/polynomials/algebraic_types.mojo RealAlgebraic -->

<!-- api: src/polynomials/algebraic.mojo algebraic_add -->

<!-- api: src/polynomials/algebraic.mojo algebraic_compare -->

<!-- api: src/polynomials/algebraic_evaluation.mojo algebraic_evaluate -->

<!-- api: src/polynomials/algebraic_specialization.mojo algebraic_evaluate_sparse -->

<!-- api: src/polynomials/algebraic.mojo algebraic_inverse -->

<!-- api: src/polynomials/algebraic.mojo algebraic_multiply -->

<!-- api: src/polynomials/algebraic.mojo algebraic_negate -->

<!-- api: src/polynomials/algebraic_roots.mojo algebraic_roots -->

<!-- api: src/polynomials/algebraic_specialization.mojo algebraic_roots_at -->

<!-- api: src/polynomials/algebraic.mojo algebraic_sign -->

<!-- api: src/polynomials/algebraic_specialization.mojo algebraic_specialize -->

<!-- api: src/polynomials/algebraic.mojo algebraic_subtract -->

<!-- api: src/polynomials/algebraic_arithmetic.mojo check_algebraic_arithmetic -->

<!-- api: src/polynomials/algebraic.mojo check_algebraic_comparison -->

<!-- api: src/polynomials/algebraic_specialization.mojo check_algebraic_evaluate_sparse -->

<!-- api: src/polynomials/algebraic_evaluation.mojo check_algebraic_evaluation -->

<!-- api: src/polynomials/algebraic_roots.mojo check_algebraic_roots -->

<!-- api: src/polynomials/algebraic_specialization.mojo check_algebraic_roots_at -->

<!-- api: src/polynomials/algebraic_specialization.mojo check_algebraic_specialization -->

<!-- api: src/polynomials/algebraic.mojo check_real_algebraic -->

<!-- api: src/polynomials/algebraic.mojo real_rational -->

<!-- api: src/polynomials/algebraic.mojo real_root -->
