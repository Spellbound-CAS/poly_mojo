# Polynomial rings

## Rings and coefficients

**Definition R-1.** A *coefficient ring* $K$ is one of $\mathbb Z$,
$\mathbb Q$, $\mathbb Q(i)$, or $\mathbb F_p$ for a prime $p < 2^{62}$. For a finite set
of variable indices $V \subset \mathbb Z$, $K[X_V]$ is the ring of
polynomials in the indeterminates $X_v$, $v \in V$: finitely supported
functions from monomials to $K$, a *monomial* being a function
$m : V \to \mathbb N$, written $X^m = \prod_v X_v^{m(v)}$, with pointwise
addition and the convolution product. The *support* of $f$ is
$\operatorname{supp}(f) = \{ m \mid f(m) \neq 0 \}$, and its *total degree*
is $\max_{m \in \operatorname{supp} f} \sum_v m(v)$ ($-\infty$ for $0$).

$\mathbb Z[X_V]$ and $\mathbb Q[X_V]$ are integral domains and unique
factorization domains; $\mathbb F_p[X_V]$ and $\mathbb Q[X_V]$ are
polynomial rings over fields [von zur Gathen and Gerhard 2013, ch. 2 and
6]. For one variable $x$, $K[x]$ over a field is a Euclidean domain.

In the code a coefficient ring is a value, its *context*:
nothing for $\mathbb Z$, $\mathbb Q$ and $\mathbb Q(i)$, the prime $p$ for
$\mathbb F_p$. Every polynomial carries the context of its coefficients,
so that $0$ and $1$ of $\mathbb F_p$ are always built with the right $p$.

**Definition R-1a (Gaussian rationals).**
$\mathbb Q(i)=\{a+bi:a,b\in\mathbb Q\}$, where $i^2=-1$.
Addition is component-wise and
$(a+bi)(c+di)=(ac-bd)+(ad+bc)i$. The representation is unique; zero means
both components are zero. A nonzero element has the exact inverse
$(a-bi)/(a^2+b^2)$, since $a^2+b^2>0$. Thus these values form a field,
and their polynomial rings have the usual field-based arithmetic and
Gröbner properties. QQi stores the pair using apn_mojo ExactComplex with
canonical Rational components and the unique NoContext context.

## Monomial orders

**Definition R-2.** A *monomial order* is a total order $\prec$ on
monomials that is a well-order, with $1 \preceq X^m$ for all $m$, and
$X^a \prec X^b \Rightarrow X^{a+c} \prec X^{b+c}$. Given a list of
variables $x_1 > x_2 > \cdots > x_n$:

- `LEX`: compare exponent vectors $(m(x_1), \dots, m(x_n))$
  lexicographically;
- `GRLEX`: total degree first, then `LEX`;
- `GREVLEX`: total degree first; on a tie, the monomial with the smaller
  exponent in the *last* variable where they differ is larger.

These are monomial orders [Cox, Little and O'Shea 2015, ch. 2, §2].

**Definition R-3 (the kernel's order).** `KERNEL` orders monomials, each
written as the list of its pairs $(v, m(v))$ with $m(v) > 0$ and $v$
increasing, lexicographically pair by pair, a proper prefix first; pairs
compare by index, then exponent. It is the order of the kernel's
`FORMAL_POLYNOMIAL` encoding (Definition P-2 of the kernel's theory).

**Proposition R-1.** `KERNEL` is a strict total order on monomials in
which $1$ is least, but it is not a monomial order: it is neither
compatible with multiplication nor a well-order. It is used for canonical
storage only, never for division or Gröbner bases.

*Proof.* Lexicographic order on finite words over a totally ordered
alphabet, with a proper prefix first, is a strict total order, and the
empty word ($1$) is a prefix of every word. Multiplication: $X_0 \prec X_1$
(first pairs $(0,1) \prec (1,1)$), but $X_0 \cdot X_0 = X_0^2 \succ X_0 X_1 = X_1 \cdot X_0$
(first pairs $(0,2) \succ (0,1)$). Well-order: indices range over
$\mathbb Z$, and $X_{-1} \succ X_{-2} \succ X_{-3} \succ \cdots$ descends
forever. (Even over non-negative indices,
$m_k = X_0 X_1 \cdots X_{k-1} X_k^2$ satisfies $m_k \succ m_{k+1}$, by the
pairs at position $k$: $(k, 2) \succ (k, 1)$.) $\square$

## Canonical forms

**Definition R-4.** The *canonical form* of $f$ in order $\prec$ is the
list of pairs $(m, f(m))$ for $m \in \operatorname{supp}(f)$, sorted by
$\prec$, each coefficient in its canonical form (an `Integer`, a reduced
`Rational` with positive denominator, a `ModInt` in $[0, p)$).

**Proposition R-2.** For a fixed strict total order, every polynomial has
exactly one canonical form, and two polynomials are equal exactly when
their canonical forms are.

*Proof.* The support is unique, so is its sorted listing, and so are the
coefficients' canonical forms. For QQi the canonical coefficient is its
ordered pair of canonical rationals. $\square$

**Requirement R-R1 (representation invariant).** A `SparsePoly` always
holds a canonical form: strictly increasing monomials in its order, no
zero coefficient, every monomial's indices strictly increasing and
exponents positive. A `DensePoly` holds coefficients $c_0, \dots, c_d$
with $c_d \neq 0$, or the empty list for $0$. Every operation preserves
these invariants; `check_invariants` verifies them in tests.

## Arithmetic

**Proposition R-3 (correctness of sparse arithmetic).** Merging two
canonical forms, adding coefficients of equal monomials and dropping
zeros, computes the canonical form of the sum; forming all products of
terms, adding exponent maps, multiplying coefficients, then sorting and
combining equal monomials, computes the canonical form of the product.

*Proof.* Immediate from Definition R-1 and Proposition R-2. $\square$

**Heap-merge variant.** Under a multiplicative monomial order, the products
of one fixed term with a sorted polynomial form a sorted row. A heap containing
the smallest remaining product of each row emits all nonzero term products in
monomial order. Combining adjacent equal monomials gives the same canonical
product as R-3. Choosing the shorter input for the rows bounds pending products
by $\min(m,n)$; the coefficient products still number at most $mn$. Ties use row
indices, and no hash iteration is involved. This argument does not apply to
KERNEL storage order, which retains sorting of the full product list.

**Theorem R-4 (Karatsuba, cited).** For dense univariate polynomials of
degree less than $n = 2k$, writing $f = f_1 x^k + f_0$ and
$g = g_1 x^k + g_0$,
$fg = f_1 g_1 x^{2k} + \big((f_1 + f_0)(g_1 + g_0) - f_1 g_1 - f_0 g_0\big) x^k + f_0 g_0$,
which recursively takes $O(n^{\log_2 3})$ coefficient operations [von zur
Gathen and Gerhard 2013, §8.1].

## Evaluation and substitution

**Proposition R-5.** For $f \in K[X_V]$ and a point $a \in K^V$, Horner
evaluation in each variable computes $f(a)$, the image of $f$ under the
unique $K$-algebra homomorphism $X_v \mapsto a_v$. Substituting
polynomials for variables is the homomorphism $X_v \mapsto g_v$. The
formal derivative $\partial f / \partial X_v$ satisfies the product rule.

*Proof.* The universal property of polynomial rings; the product rule is a
direct computation on monomials. $\square$

For substitution in one variable, the implementation groups only the occupied exponents
$e_0 < \cdots < e_k$. Starting with $r=c_{e_k}$, each step replaces $r$ by
$r g^{e_{i+1}-e_i}+c_{e_i}$, then the result is multiplied by $g^{e_0}$.
Expanding this sparse Horner identity gives $\sum_i c_{e_i}g^{e_i}$, including
when $g$ contains the substituted variable or is zero. Storage for coefficient
groups depends on occupied terms, not the largest exponent.

## The kernel codec

**Theorem R-6 (the codec is a bijection onto canonical encodings).** The
map from a `SparsePoly[Rational]` in `KERNEL` order to the byte encoding
of the kernel's `FORMAL_POLYNOMIAL` kind,
writing the number of terms, then each term's pairs as `i64` and `u32` and
its coefficient's rational encoding, is a bijection between polynomials
other than $0$ and $1$ and the encodings the kernel's decoder accepts.
Decoding is its inverse.

*Proof.* By Proposition R-2 a polynomial has one canonical form in
`KERNEL` order; the encoding writes it field by field and is injective;
the kernel's decoder accepts exactly the field sequences of canonical
forms (it rejects order violations, zero exponents and coefficients,
non-canonical rationals, trailing bytes, and the values $0$ and $1$).
$\square$

The codec lets the kernel's arithmetic and this package's be compared
byte for byte (contract P), and lets `polynorm` hand certificates to the
kernel as literals.
