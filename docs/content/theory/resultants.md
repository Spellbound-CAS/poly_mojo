# Resultants

**Definition E-1.** For $f = \sum_{i=0}^{m} a_i x^i$ and
$g = \sum_{j=0}^{n} b_j x^j$ over a commutative ring $R$, with $m, n \ge 1$
the formal degrees, the *Sylvester matrix* $\operatorname{Syl}(f, g)$ is
the $(m + n) \times (m + n)$ matrix whose first $n$ rows hold the
coefficients of $x^{n-1} f, \dots, xf, f$ and whose last $m$ rows hold
those of $x^{m-1} g, \dots, xg, g$, each row listing the coefficients of
$x^{m+n-1}, \dots, x, 1$ in that order; the *resultant* is
$\operatorname{res}(f, g) = \det \operatorname{Syl}(f, g) \in R$. For
multivariate polynomials, $\operatorname{res}_{x}$ takes coefficients in
the ring of the other variables.

This is the usual convention: over a field containing the roots, with
$f = a_m \prod_i (x - \alpha_i)$ and $g = b_n \prod_j (x - \beta_j)$,
$\operatorname{res}(f, g) = a_m^n b_n^m \prod_{i, j} (\alpha_i - \beta_j)$
[von zur Gathen and Gerhard 2013, ch. 6]. Listing the rows of each block in
the opposite order multiplies the determinant by
$(-1)^{n(n-1)/2 + m(m-1)/2}$, which is $-1$ already for $m = 1$, $n = 2$.

**Theorem E-2 (common factors, cited).** Over a UFD $R$, for $f, g$ with
$a_m b_n \neq 0$: $\operatorname{res}(f, g) = 0$ if and only if $f$ and $g$
have a common factor of positive degree in $x$ [Cox, Little and O'Shea
2015, ch. 3, §6, Proposition 8; von zur Gathen and Gerhard 2013, §6.3].

**Theorem E-3 (the resultant is in the ideal).** There are
$u, v \in R[x]$ with $\deg u < n$, $\deg v < m$, and
$u f + v g = \operatorname{res}(f, g)$.

*Proof.* The columns of $S = \operatorname{Syl}(f, g)$ run by decreasing
power of $x$, so the last column holds the constant terms. For $u$ of
degree less than $n$ and $v$ of degree less than $m$, the row vector of
their coefficients (each from the highest power down), multiplied by $S$,
is the coefficient vector of $uf + vg$. The adjugate satisfies $\operatorname{adj}(S)\, S = \det(S)\, I$,
so the last row $w$ of $\operatorname{adj}(S)$ satisfies
$w S = (0, \dots, 0, \det S)$. Reading $w$ as the coefficients of $(u, v)$
gives $uf + vg = \det S = \operatorname{res}(f, g)$, with every entry of $w$
a minor of $S$, hence in $R$ [Cox, Little and O'Shea 2015, ch. 3, §6,
Proposition 9]. $\square$

**Theorem E-4 (subresultants, cited).** The subresultants
$S_j(f, g)$, for $0 \le j < \min(m, n)$, are polynomials of degree at most
$j$ whose coefficients are minors of $\operatorname{Syl}(f, g)$;
$S_0 = \operatorname{res}(f, g)$, and $\deg \gcd(f, g)$ is the least $j$
with $S_j \neq 0$; the subresultant PRS computes them [Collins 1967;
Brown and Traub 1971].

**Algorithm E-5.** The resultant is computed by the subresultant PRS over
$R$ (no division outside $R$), with the cofactors $u$, $v$ of Theorem E-3
carried along the sequence, as in the extended Euclidean algorithm.

**Certificate E-C1 (resultant identity).** Data: $u$, $v$, $s$. Check:
$u f + v g = s$. Claim: $s$ lies in the ideal $\langle f, g \rangle$; if
$s$ is nonzero and free of $x$, it serves in the coprimality certificate
G-C3. That $s$ *is* the resultant is not claimed and not needed.
