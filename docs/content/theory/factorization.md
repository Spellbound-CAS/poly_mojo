# Factorization

Factorization returns a scalar and irreducible factors with their exact
multiplicities over the integers, rationals, Gaussian rationals or a prime
field. Dense inputs are univariate; sparse inputs can contain several
variables. The functional API also factors univariate and multivariate polynomials over
`AA` and `QQbar`. Core product checkers verify reconstruction independently; the
algebraic-field facade check additionally establishes irreducibility.

## Squarefree factorization

**Definition Z-1.** A polynomial is *squarefree* when no irreducible
factor divides it twice. The *squarefree decomposition* of
$f \in K[x]$, $K$ of characteristic $0$, is $f = c\, a_1 a_2^2 \cdots a_k^k$
with $c \in K$, each $a_i$ squarefree and monic, and the $a_i$ pairwise
coprime; it is unique.

**Theorem Z-2 (Yun's algorithm, cited).** For $K$ of characteristic $0$,
with $b_1 = f/\gcd(f, f')$, $c_1 = f'/\gcd(f, f')$, $d_1 = c_1 - b_1'$, and,
for $i \ge 1$, $a_i = \gcd(b_i, d_i)$, $b_{i+1} = b_i/a_i$,
$c_{i+1} = d_i/a_i$, $d_{i+1} = c_{i+1} - b_{i+1}'$, until $b_i = 1$, the
$a_i$ are the factors of the squarefree decomposition [Yun 1976; von zur
Gathen and Gerhard 2013, §14.6].

**Certificate Z-C1 (squarefree decomposition).** Data: $c$, $a_1, \dots, a_k$.
Check: $f = c \prod_i a_i^i$. Claim: this product equation. That each
$a_i$ is squarefree and that they are pairwise coprime are separate claims,
certified when needed by G-C2 or G-C3 for $(a_i, a_i')$ and $(a_i, a_j)$.

## Factorization over a finite field

**Theorem Z-3 (cited).** A squarefree $f \in \mathbb F_p[x]$ factors into
irreducibles by distinct-degree factorization (gcds with
$x^{p^d} - x$) followed by equal-degree factorization, either by
Cantor and Zassenhaus's randomized splitting or by Berlekamp's algorithm
[Cantor and Zassenhaus 1981; Berlekamp 1967; von zur Gathen and Gerhard
2013, ch. 14]. The randomized algorithm, run with a fixed seed, is
deterministic in its output: the factorization is unique, so only
the running time depends on the seed.

## Factorization over the integers

**Theorem Z-4 (Hensel lifting, cited).** Let $f \in \mathbb Z[x]$ be
primitive and squarefree, $p$ a prime not dividing $\operatorname{lc}(f)$
with $f \bmod p$ squarefree, and $f \equiv \operatorname{lc}(f)\, g_1 \cdots g_r \pmod p$
with the $g_i$ monic and pairwise coprime modulo $p$. For every $k \ge 1$
there are unique monic $G_i \in \mathbb Z[x]$, defined modulo $p^k$, with
$G_i \equiv g_i \pmod p$ and $f \equiv \operatorname{lc}(f)\, G_1 \cdots G_r \pmod{p^k}$
[Zassenhaus 1969; von zur Gathen and Gerhard 2013, §15.4].

**Theorem Z-5 (Landau–Mignotte bound, cited).** If $g \in \mathbb Z[x]$
divides $f \in \mathbb Z[x]$ and $\deg g = d$, then every coefficient of
$g$ is at most $2^d \lVert f \rVert_2 \cdot \lvert \operatorname{lc}(g)/\operatorname{lc}(f)\rvert$
in absolute value, $\lVert f \rVert_2$ being the Euclidean norm of $f$'s
coefficient vector [Mignotte 1974; von zur Gathen and Gerhard 2013, §6.6].

**Algorithm Z-6 (Zassenhaus's algorithm).**

1. Write $f = c \cdot \operatorname{pp}(f)$ and factor the primitive part's
   squarefree decomposition (Theorem Z-2); factor each squarefree factor
   $h$ as follows.
2. Choose the first prime $p$ in the fixed sequence with $p \nmid \operatorname{lc}(h)$
   and $h \bmod p$ squarefree; factor $h \bmod p$ (Theorem Z-3).
3. Lift the factorization to $p^k$ with $p^k > 2B\,\lvert\operatorname{lc}(h)\rvert$,
   $B$ the bound of Theorem Z-5 for factors of degree at most
   $\deg h / 2$ (Theorem Z-4).
4. For subsets $S$ of the lifted factors, in increasing size, form
   $\operatorname{lc}(h) \prod_{i \in S} G_i$ reduced to symmetric
   coefficients modulo $p^k$, take its primitive part, and test whether it
   divides $h$ (Proposition V-3); on success, record it, remove $S$, and
   divide.
5. What remains is irreducible.

**Theorem Z-7 (correctness).** Algorithm Z-6 returns the factorization of
$f$ into its content and irreducible primitive factors with positive
leading coefficients, with multiplicities.

*Proof.* Every irreducible factor $q$ of $h$ over $\mathbb Z$ reduces modulo
$p$ to a product of some of the $g_i$ (unique factorization in
$\mathbb F_p[x]$), hence, by uniqueness in Theorem Z-4, to the product of the
corresponding $G_i$ modulo $p^k$ up to the leading coefficient. Since the
coefficients of $\operatorname{lc}(h) q / \operatorname{lc}(q)$ are bounded
by $B \lvert \operatorname{lc}(h)\rvert < p^k / 2$, the symmetric
representative of $\operatorname{lc}(h) \prod_{i \in S} G_i$ equals it
exactly, and its primitive part is $q$. Trying subsets in increasing size
finds the irreducible factors in order, and the trial division of step 4
admits only true divisors. Gauss's lemma (Theorem G-1) relates the
factorizations over $\mathbb Z$ and $\mathbb Q$. $\square$

The subset search is exponential in the worst case (Swinnerton-Dyer
polynomials); van Hoeij's lattice-based recombination removes that cost
and is available for large splitting sets. It does not change the certificates.

**Certificate Z-C2 (factorization).** Data: $c$ and the pairs
$(f_i, e_i)$. Check: $f = c \prod_i f_i^{e_i}$. Claim: this product
equation. The core checker does not certify irreducibility of the $f_i$: an
answer is correct as an *equation* even when a supplied factor splits further.
The `AA` and `QQbar` facade check below uses the same record with additional
degree and sign conditions.

**Proposition Z-8 (over the rationals).** For $f \in \mathbb Q[x]$, write
$f = (a/b) F$ with $F \in \mathbb Z[x]$ primitive; the factorization of $f$
is $a/b$ times that of $F$, with the factors made monic or primitive as the
caller asks. $\square$


## Positive characteristic and sparse multiplicity blocks

For monic $f$, set $c=\gcd(f,f')$ and $w=f/c$. At multiplicity $m=1,2,\ldots$,
compute $y=\gcd(w,c)$, emit $w/y$ when it is nonconstant, and replace
$(w,c)$ by $(y,c/y)$. Over a prime field, factors whose multiplicity is divisible
by $p$ remain in $c$. All their exponents are divisible by $p$; divide the
exponents by $p$, repeat the decomposition, and multiply the resulting
multiplicities by $p$. Coefficients stay unchanged because Frobenius is the
identity on $\mathbb F_p$. An explicit loop carries this multiplicity scale
through arbitrarily many derivative-free stages. See the squarefree reduction
in [von zur Gathen and Panario's survey](https://people.csail.mit.edu/dmoshkov/courses/codes/poly-factorization.pdf).

For several variables, replace the initial GCD by
$c=\gcd(f,\partial_{x_1}f,\ldots,\partial_{x_n}f)$. A factor occurring with
multiplicity $m$ contributes $m-1$ to this GCD unless the characteristic
divides $m$. The same recurrence therefore applies. In characteristic $p$,
inverse Frobenius divides every variable exponent by $p$. The sparse routine
avoids a whole-polynomial Kronecker image; its internal GCDs may still allocate
dense univariate or coefficient-tower views. All blocks are monic in the input's
stored order and are returned in increasing multiplicity. `KERNEL` storage
order is not multiplicative, so decomposition uses lex internally, then
renormalizes returned blocks and adjusts the scalar by each block's full power.

## Gaussian-rational factorization

For a monic squarefree block $f\in\mathbb Q(i)[x]$, try integer shifts
$g(x)=f(x+s i)$, starting at $s=0$, until the norm
$N(g)=g\overline g$ is squarefree over $\mathbb Q$. Writing $g=a+i b$ gives
$N(g)=a^2+b^2$, so the implementation uses two rational polynomial squares
and a sum. The norm is factored over $\mathbb Q$. Each rational irreducible
$h$ supplies $\gcd(g,h)$ over $\mathbb Q(i)$; nonconstant GCDs are shifted
back by $-s i$ to recover the factors of $f$.

The squarefree-norm condition separates conjugate factors, making these GCDs
irreducible over the extension. This is Trager's reduction, described in the
[SymPy polynomial internals](https://docs.sympy.org/latest/modules/polys/internals.html#sympy.polys.factortools.dup_ext_factor).
Squarefree blocks supply the original multiplicities. Every shift, rational
factorization and GCD shares the caller's budget.

## Sparse field factorization

The QQi and prime-field sparse routines first separate squarefree blocks,
then apply a complete Kronecker reduction to each block. This removes repeated
factors before they inflate the image and subset search. For a block $f$ with
occurring variables $x_0,\ldots,x_{n-1}$, let
$b_j=1+\deg_{x_j}f$ and $w_j=\prod_{k<j}b_k$. Substituting $x_j=t^{w_j}$
encodes each monomial exponent vector as its mixed-radix integer. The image
uses at most $\prod_j b_j$ dense coefficient slots. Checked size arithmetic
and budget charges precede that allocation.

Every divisor has per-variable degrees bounded by those of $f$, so its image
decodes uniquely inside the same degree box. Factor the image, enumerate
subsets of its irreducible factor occurrences, decode their products, and
accept only exact sparse divisors. Trying smaller subsets first finds
irreducible factors; if a remaining polynomial were reducible, one of its
proper factors would use at most half the remaining image-factor occurrences
and would have been tried. Repeated image factors are retained as occurrences,
then equal sparse factors are combined into multiplicities.

A split image does not imply a split input: decoded candidates must pass
exact division. Image storage can grow rapidly with the number of variables,
and recombination is exponential in the worst case. The algorithm offers
complete exact coverage within the work allowance; it does not promise
efficient factorization for every sparse input.


## Real and complex algebraic fields

Every nonconstant univariate polynomial over $\overline{\mathbb Q}$ splits into linear
factors: its roots are algebraic over a finite extension of $\mathbb Q$, hence
algebraic over $\mathbb Q$. This is the `QQbar` factorization used by the
functional API. For real algebraic coefficients (`AA`), non-real roots occur
in conjugate pairs with matching multiplicity. If $a+bi$ is such a root, the
pair contributes

$$(x-a-bi)(x-a+bi)=x^2-2ax+(a^2+b^2).$$

Its discriminant is $-4b^2<0$, so it is irreducible over the real field. Real
roots contribute linear factors. These exhaust the possibilities; an
irreducible univariate polynomial over `AA` has degree one or two.

Consequently, a checker can establish complete irreducible factorization
without repeating the producer's root search. Validate the coefficient
witnesses, require monic distinct factors, check the degree of each factor and
the negative discriminant of each real quadratic, and verify the powered
product with positive multiplicities. Unique factorization then gives the
reported multiplicities and completeness. The implementation also checks the
weighted degree sum before expansion to bound malformed certificates.

See [the tutorial](../tutorial/algebraic-coefficients.md#factor-over-real-or-complex-algebraic-numbers)
for choosing a field through coefficient values.


### Completeness in several variables

Multivariate polynomials over an algebraically closed field need not split into
linear factors: $x^2+y$ is irreducible over `QQbar`. The algebraic sparse path
uses the same reduction to univariate factorization as Kronecker's algorithm;
see section 5 of [Kaltofen's factorization survey](https://users.cs.duke.edu/~elk27/bibliography/82/Ka82_survey.pdf).

Let $f=gh$ over either field. Degrees in each variable add, so both divisors lie
in the box $0\leq\deg_{x_j}g,\deg_{x_j}h\leq\deg_{x_j}f$. Mixed-radix encoding
is injective on that box and is a ring substitution. Consequently, the image
of each divisor is a scalar times a subproduct of the image's univariate
irreducible factors. One side has image degree at most half the total. Trying
every such subproduct, decoding it, and testing exact sparse divisibility must
find a proper divisor whenever one exists. Failed decoding candidates do not
weaken this completeness argument.

The checker receives the image factorizations rather than recomputing them.
Their irreducibility follows from the univariate degree and discriminant tests
above; their product establishes completeness. An independent exhaustive search
then establishes multivariate irreducibility. Exact reconstruction from distinct
irreducible factors establishes the output multiplicities. This evidence avoids
repeating the image-root searches, but it does not make verification polynomial-time.
