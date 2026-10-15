# Theory

These pages state the mathematics this package implements: definitions,
algorithms, and the theorems that make them correct, with proofs or, for
classical results, citations. They also define each certificate and prove
that checking it establishes the claim it supports. The following labels distinguish definitions,
proofs, algorithms, and implementation requirements:

| Label | Meaning |
|---|---|
| **Definition** | Fixes a notion; the code implements exactly this. |
| **Theorem**, **Lemma**, **Proposition** | Proved here. |
| **Theorem (cited)** | A classical result, with a reference. |
| **Algorithm** | Numbered steps, with the theorem that proves it correct. |
| **Certificate** | Data, the relation a checker verifies, and the claim that relation implies. |
| **Requirement** | A property of the implementation that the proofs assume, with its tests. |

Results are numbered per page: R for rings, V for division, G for gcds,
E for resultants, F for rational functions, Z for factorization, B for
Gröbner bases, T for real roots, A for the [algebraic/CAD extensions](extensions.md),
and N for [numerical enclosures and roots](numerical.md).

## What justifies what

| Component | Results | Page |
|---|---|---|
| Representations, orders, canonical forms, the kernel codec | R-1 to R-6 | [Polynomial rings](rings.md) |
| Division with remainder, pseudo-division, multivariate division | V-1 to V-5 | [Division](division.md) |
| Content, primitive parts, gcds, modular gcds, Bézout and coprimality certificates | G-1 to G-9 | [Greatest common divisors](gcd.md) |
| Resultants and subresultants | E-1 to E-5 | [Resultants](resultants.md) |
| Normalized fractions and cancellation | F-1 to F-3 | [Rational functions](rational-functions.md) |
| Squarefree and univariate factorization | Z-1 to Z-8 | [Factorization](factorization.md) |
| Gröbner bases and membership certificates | B-1 to B-7 | [Gröbner bases](groebner.md) |
| Sturm sequences and real root isolation | T-1 to T-6 | [Real roots](roots.md) |
| Numerical evaluation and uniform real/complex root covers | N-1 to N-6 | [Numerical guarantees](numerical.md) |

## What the kernel checks

None of these algorithms is trusted. The numerical extension has independent
library checkers and does not extend the kernel codec. For exact claims supported by the integration, the kernel checks
always an identity between polynomials, such as $p = g p'$ or
$f = \sum_i h_i g_i$, by the `FORMAL_POLYNOMIAL` kind's arithmetic,
together with, for real roots, finitely many sign evaluations at rational
points. Each certificate on these pages is designed so that its check is
of that form, and each check is proved to imply its claim. That is why the
algorithms can be as clever as speed requires.
