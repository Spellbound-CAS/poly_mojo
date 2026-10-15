# Rational functions

**Definition F-1.** A *rational function* over $\mathbb Q$ is an element of
the fraction field $\mathbb Q(X)$ of $\mathbb Q[X]$. Its *normalized form*
is the pair $(p, q)$ with $p/q$ the element, $q \neq 0$,
$\gcd(p, q) = 1$, and $q$ with leading coefficient $1$ in the chosen
monomial order (and $q = 1$ when $p = 0$).

**Proposition F-1.** Every rational function has exactly one normalized
form.

*Proof.* If $p/q = p'/q'$ in normalized forms, $pq' = p'q$; since
$\gcd(p, q) = 1$, $q \mid q'$, and symmetrically $q' \mid q$, so they are
associates; both being monic, $q = q'$, and then $p = p'$. Existence: divide
any representation by the gcd of numerator and denominator and by the
denominator's leading coefficient. $\square$

**Theorem F-2 (cancellation).** If $p = g p'$, $q = g q'$ (Certificate
G-C1) and $q \neq 0$, then $p/q = p'/q'$ in $\mathbb Q(X)$. Pointwise, for
every point $a$ with $q(a) \neq 0$, also $q'(a) \neq 0$ and
$p(a)/q(a) = p'(a)/q'(a)$.

*Proof.* In the field, $p/q = (g p')/(g q') = p'/q'$ since $g \neq 0$.
At a point, $q(a) = g(a) q'(a) \neq 0$ forces both factors nonzero. $\square$

The pointwise statement is what a CAS answer needs: in the `scalar`
profile it is the transfer principle (Theorem X-6 of the kernel's theory
pages) with the guard $q \neq 0$ and the finiteness of the variables.
Points where $q' \ne 0$ but $q = 0$ (the removed singularities, such as
$x = 1$ in $(x^2 - 1)/(x - 1)$) are exactly where the original and the
cancelled forms disagree as functions on the carrier, which is why
`polynorm` states the guard $q \neq 0$, not $q' \neq 0$.

**Proposition F-3 (arithmetic).** For normalized $p_1/q_1$ and
$p_2/q_2$, the normalized forms of their sum and product are obtained by
forming $p_1 q_2 + p_2 q_1$ over $q_1 q_2$, respectively $p_1 p_2$ over
$q_1 q_2$, and normalizing; the gcd can be computed more cheaply from
$\gcd(q_1, q_2)$ (Henrici's method) without changing the result, by
Proposition F-1. $\square$

**Certificate F-C1 (normalization).** Data: the cofactors of G-C1 for the
cancellation, and optionally G-C2 or G-C3 showing that the result is fully
reduced. The claim needed for an equation is Theorem F-2's; full
reduction is a quality claim.


The array facade requires G-C3 for every normalized function, so
`poly.rational.check(f)` proves full reduction as well as the cancellation
identity and monic denominator. The same argument applies over QQi, AA, and
QQbar. The [rational-function reference](../reference/array-rational.md) lists
the multivariate coprimality and partial-fraction irreducibility checks.
