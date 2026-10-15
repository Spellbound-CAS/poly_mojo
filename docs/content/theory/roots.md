# Real roots

Real root isolation serves real-algebraic certificates: MetiTarski's
univariate decisions, sign conditions in
`linarith`'s nonlinear extensions, and the root counts that the CAS
reports. All computations are exact, in $\mathbb Q$.

**Definition T-1 (Sturm sequence).** For $f \in \mathbb Q[x]$ of
positive degree, its *Sturm sequence* is $f_0 = f$, $f_1 = f'$, and
$f_{i+1} = -\operatorname{rem}(f_{i-1}, f_i)$ while $f_i \neq 0$, ending
with $f_\ell \neq 0$, $f_{\ell+1} = 0$. For a point $a$, $V(a)$ is the
number of sign changes in $f_0(a), \dots, f_\ell(a)$, zeros omitted.

**Theorem T-2 (Sturm's theorem, cited).** For $a < b$ with
$f(a) f(b) \neq 0$, the number of distinct real roots of $f$ in $(a, b)$
is $V(a) - V(b)$. This holds whether or not $f$ is squarefree: $f_\ell$ is
then $\gcd(f, f')$ up to a constant, and dividing the sequence by it does
not change the sign changes at points where $f$ does not vanish [Sturm
1829; Basu, Pollack and Roy 2006, ch. 2].

**Proposition T-3 (Cauchy's bound).** Every complex root $z$ of
$f = \sum_{i=0}^d c_i x^i$, $c_d \neq 0$, satisfies
$\lvert z \rvert < 1 + \max_{i < d} \lvert c_i / c_d \rvert$.

*Proof.* If $\lvert z \rvert \ge 1 + M$ with $M$ the maximum, then
$\lvert \sum_{i<d} c_i z^i \rvert \le M \lvert c_d\rvert (\lvert z\rvert^d - 1)/(\lvert z\rvert - 1) < \lvert c_d z^d \rvert$,
so $f(z) \neq 0$. $\square$

**Algorithm T-4 (isolation by bisection).** Start from $(-R, R)$ with $R$
the bound of T-3, a power of 2. Count roots with Sturm's theorem; split
intervals with more than one root at their midpoint (if the midpoint is a
root, record it as an exact root and continue on both halves); keep
intervals with exactly one. The result is a list of disjoint open
intervals with rational (dyadic) endpoints, each containing exactly one
root, plus the exact rational roots found. Descartes' rule with the
Vincent–Collins–Akritas subdivision is a faster alternative with the same
output [Collins and Akritas 1976].

**Theorem T-5.** Algorithm T-4 terminates and isolates every real root.

*Proof.* The distinct real roots are finitely many, at positive distance
from each other; after enough bisections every interval contains at most
one, and every split preserves the counts by Theorem T-2. $\square$

**Certificate T-C1 (root count).** Data: the Sturm sequence
$f_0, \dots, f_\ell$, the quotients $q_i$ with
$f_{i-1} = q_i f_i - f_{i+1}$, and the points $a < b$. Check: $f_0 = f$,
$f_1 = f'$, each relation $f_{i-1} = q_i f_i - f_{i+1}$ for
$1 \le i \le \ell$ with $f_{\ell+1} = 0$ (polynomial identities; the last
says that $f_\ell$ divides $f_{\ell-1}$, hence every $f_i$, so $f_\ell$ is
$\gcd(f, f')$ up to a constant), $f(a) f(b) \neq 0$, and the values
$f_i(a)$, $f_i(b)$ by exact rational evaluation. Claim: $f$ has exactly
$V(a) - V(b)$ distinct real roots in $(a, b)$. The relations, not the
degrees, carry the proof: at a root $c$ of $f_i$ with $0 < i < \ell$ they
give $f_{i-1}(c) = -f_{i+1}(c)$, the sign alternation on which Sturm's
theorem rests.

The kernel can check the identities and the evaluations; the step from
them to the root count is Sturm's theorem, which must be a theorem of the
proof checker's library of real-algebraic results. Until it is, root counts
enter proofs as hypotheses with provenance.

**Theorem T-6 (sign on an interval).** If $f$ has no root in $[a, b]$
(T-C1 with count $0$ on a slightly larger interval, or $f(a) f(b) \neq 0$
and count $0$ on $(a, b)$) and $f(c) > 0$ for some $c \in [a, b]$, then
$f > 0$ on $[a, b]$.

*Proof.* The intermediate value theorem for polynomials. $\square$
