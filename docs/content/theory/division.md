# Division

## Univariate division

**Theorem V-1 (division with remainder).** Let $K$ be a field and
$f, g \in K[x]$ with $g \neq 0$. There are unique $q, r \in K[x]$ with
$f = qg + r$ and $\deg r < \deg g$.

*Proof.* Existence by induction on $\deg f$: if $\deg f < \deg g$ take
$q = 0$, $r = f$; otherwise $f' = f - (\operatorname{lc}(f)/\operatorname{lc}(g))\, x^{\deg f - \deg g} g$
has smaller degree. Uniqueness: from $qg + r = q'g + r'$,
$(q - q')g = r' - r$, whose degree is less than $\deg g$, so $q = q'$.
$\square$

**Algorithm V-1.** The schoolbook loop of the proof, on dense
coefficient vectors, with $O(\deg g \cdot (\deg f - \deg g + 1))$ field
operations.

**Theorem V-2 (pseudo-division).** Let $R$ be an integral domain,
$f, g \in R[x]$, $g \neq 0$, $\delta = \max(\deg f - \deg g + 1, 0)$. There
are unique $q, r \in R[x]$ with
$\operatorname{lc}(g)^{\delta} f = qg + r$ and $\deg r < \deg g$.

*Proof.* Apply Theorem V-1 over the fraction field to
$\operatorname{lc}(g)^\delta f$: every step divides only by
$\operatorname{lc}(g)$, and the factor $\operatorname{lc}(g)^\delta$
supplies one such factor per step, so $q$ and $r$ lie in $R[x]$.
Uniqueness as in V-1. $\square$

Pseudo-division keeps computations in $\mathbb Z[x]$ when coefficients
must stay integral, as in the subresultant algorithm ([gcd](gcd.md)).

**Proposition V-3 (exact division over the integers).** For
$f, g \in \mathbb Z[x]$, $g \neq 0$: $g$ divides $f$ in $\mathbb Z[x]$ if
and only if division over $\mathbb Q$ gives remainder $0$ and a quotient
with integer coefficients.

*Proof.* If $f = qg$ with $q \in \mathbb Z[x]$, uniqueness in $\mathbb Q[x]$
gives that quotient and remainder $0$. The converse is the definition.
$\square$

**Certificate V-C1 (division).** Data: $q$, $r$. Check: $f = qg + r$
and $\deg r < \deg g$. Claim: $r$ is the remainder of $f$ by $g$; in
particular $g \mid f$ when $r = 0$. The check is one multiplication and one
addition; by uniqueness (V-1) it determines the claim.

## Division by several polynomials

**Theorem V-4 (multivariate division, cited).** Fix a monomial order
$\prec$ and $g_1, \dots, g_s \in K[X]$, $K$ a field, all nonzero. Every
$f \in K[X]$ can be written $f = a_1 g_1 + \cdots + a_s g_s + r$ where no
monomial of $r$ is divisible by any $\operatorname{LT}(g_i)$, and
$\operatorname{LM}(a_i g_i) \preceq \operatorname{LM}(f)$ whenever
$a_i g_i \neq 0$. The standard algorithm computes such $a_i, r$ [Cox,
Little and O'Shea 2015, ch. 2, §3, Theorem 3].

**Proposition V-5 (termination).** The multivariate division algorithm
terminates, because each step strictly decreases the leading monomial of
the remaining dividend, and a monomial order is a well-order.

The quotients and remainder are not unique in general: they depend on the
order of the $g_i$. Uniqueness of the remainder returns when the $g_i$
form a Gröbner basis ([Gröbner bases](groebner.md)).

**Implementation invariant.** Multiplication by a fixed monomial
preserves an admissible order. Thus the terms of a divisor times one quotient
term can be merged directly with the pending dividend, cancelling equal
monomials without constructing a temporary product polynomial. By V-5, each
divisor's quotient terms and the extracted remainder terms arrive in strictly
decreasing order. Reversing each list yields its canonical increasing order;
there are no duplicate monomials to combine. KERNEL order is not multiplicative
and uses the established canonical arithmetic path when the private fused
helper is called directly.

**Certificate V-C2 (division by a list).** Data: $a_1, \dots, a_s$,
$r$. Check: $f = \sum_i a_i g_i + r$. Claim: $f - r$ lies in the ideal
$\langle g_1, \dots, g_s \rangle$. The degree conditions of V-4 are not
part of the claim and need not be checked.
