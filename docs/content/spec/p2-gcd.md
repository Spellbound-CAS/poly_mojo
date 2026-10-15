# P2 Division, gcd, resultants, rational functions

Implemented and tested; see the [API reference](../reference/p2.md).

Algorithm provenance and independent regression cases: [reference material](../contributing/reference.md#p2-division-gcd-resultants-rational-functions).

## P2.1 Division

```text
struct DivisionResult[C]:  quotient: DensePoly[C], remainder: DensePoly[C]
def divmod[C: CoefficientField](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> DivisionResult[C]
def pseudo_divmod[C: CoefficientRing](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> DivisionResult[C]
def divide_exact[C](f: SparsePoly[C], g: SparsePoly[C], mut budget: Budget) raises PolyError -> Optional[SparsePoly[C]]
struct ListDivision[C]:    quotients: List[SparsePoly[C]], remainder: SparsePoly[C]
def divide_by_list[C: CoefficientField](f: SparsePoly[C], divisors: List[SparsePoly[C]], mut budget: Budget) raises PolyError -> ListDivision[C]
```

`divmod` and `pseudo_divmod` follow Theorems V-1 and V-2 (`DOMAIN` when
`g` is zero); `divide_exact` returns the quotient when `g` divides `f`
(over $\mathbb Z$, Proposition V-3); `divide_by_list` is Theorem V-4's
algorithm in the divisors' common order.

**Required checks:** the identities $f = qg + r$ (with the degree bound) and
$\operatorname{lc}(g)^\delta f = qg + r$ hold on generated inputs;
`divide_exact` finds the quotient of every generated product and refuses
non-multiples; `divide_by_list` satisfies Certificate V-C2's identity.

## P2.2 Content, primitive part, gcd

```text
def content[C](p: SparsePoly[C], main_variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]
def primitive_part[C](p: SparsePoly[C], main_variable: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]
struct GcdResult[C]:
    gcd: SparsePoly[C]                  # normalized (Definition G-1)
    cofactors: (SparsePoly[C], SparsePoly[C])   # Certificate G-C1
    bezout: Optional[(SparsePoly[C], SparsePoly[C])]   # G-C2: univariate over a field
    coprimality: Optional[CoprimalityCertificate[C]]  # G-C3 for the cofactors, on request
def gcd[C](f: SparsePoly[C], h: SparsePoly[C], mut budget: Budget, certify_maximal: Bool = False) raises PolyError -> GcdResult[C]
def extended_gcd[C: CoefficientField](f: DensePoly[C], h: DensePoly[C], mut budget: Budget) raises PolyError -> GcdResult[C]
```

Univariate gcds over $\mathbb Z$ and $\mathbb Q$ use the modular algorithm
(Algorithm G-5) above a degree threshold and the subresultant PRS below
it (the measured threshold is degree 8); multivariate gcds use Algorithm G-6. `certify_maximal` asks for G-C2
or G-C3, which cost more than the cofactors.

**Required checks:** on generated $f = g a$, $h = g b$ with known $g$ and
coprime $a, b$, the result equals $g$ normalized; the cofactors satisfy
G-C1; Bézout coefficients satisfy G-C2; the modular and subresultant
algorithms agree.

## P2.3 Primes and the Chinese remainder theorem

```text
def prime(i: Int, mut budget: Budget) raises PolyError -> UInt64  # the i-th prime of the fixed sequence, descending below 2^62
def crt(residues: List[UInt64], primes: List[UInt64], mut budget: Budget) raises PolyError -> Integer   # symmetric representative
```

The sequence is the primes below $2^{62}$ in decreasing order, generated
by a deterministic Miller–Rabin test with a fixed witness set proven
correct for 64-bit integers, and cached.

**Required checks:** the first primes of the sequence match a list computed by an
independent method (the reference script); CRT reconstructs generated
integers.

## P2.4 Unlucky primes and early termination

Algorithm G-5's steps 3 and 5 are implemented exactly as specified; a
test forces an unlucky prime (by choosing inputs whose subresultant is
divisible by the first prime) and checks that the result is still
correct.

## P2.5 Resultants and coprimality certificates

```text
struct ResultantResult[C]:  resultant: SparsePoly[C], cofactors: (SparsePoly[C], SparsePoly[C])   # E-C1
def resultant[C](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> ResultantResult[C]
def subresultants[C](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> List[SparsePoly[C]]
struct ContentGcdStep[C]:
    gcd: SparsePoly[C]
    cofactors: (SparsePoly[C], SparsePoly[C])
    coprimality_node: Optional[Int]
struct CoprimalityNode[C]:
    variable: Int64
    u: SparsePoly[C], v: SparsePoly[C], s: SparsePoly[C]
    left_content_steps: List[ContentGcdStep[C]]
    right_content_steps: List[ContentGcdStep[C]]
    contents_node: Optional[Int]
struct CoprimalityCertificate[C]:
    nodes: List[CoprimalityNode[C]]  # referenced nodes precede their users
    root: Int
def coprimality_certificate[C](f: SparsePoly[C], h: SparsePoly[C], mut budget: Budget) raises PolyError -> Optional[CoprimalityCertificate[C]]
```

The content folds and terminal unit identities are explicit content
witnesses. Nodes reference only earlier nodes, content proofs descend in
variables, and the checker verifies only algebraic identities and structural
conditions.

**Required checks:** resultants vanish exactly on generated pairs with a common
factor (Theorem E-2); the cofactor identity holds; resultants agree with
the Sylvester determinant computed by fraction-free Gaussian elimination
on small inputs; coprimality certificates check (Theorem G-8) for coprime
pairs and are refused for non-coprime ones.

## P2.6 Rational functions

```text
struct RationalFunction(Copyable, Movable, Equatable, Writable):
    numerator: SparsePoly[QQ], denominator: SparsePoly[QQ]     # normalized (Definition F-1)
    @staticmethod
    def normalize(p: SparsePoly[QQ], q: SparsePoly[QQ], mut budget: Budget) raises PolyError -> (RationalFunction, CancellationCertificate)
    def add, sub, mul, div, neg  # each takes mut budget: Budget and raises PolyError
struct CancellationCertificate:  common: SparsePoly[QQ], cofactors: (SparsePoly[QQ], SparsePoly[QQ])
```

**Required checks:** normalized forms are unique (Proposition F-1) on generated
equal fractions; cancellation certificates check; the field axioms hold.

## Certificates and checkers

```text
def check_division(f, g, result, mut budget: Budget) raises PolyError -> Bool                 # V-C1
def check_list_division(f, divisors, result, mut budget: Budget) raises PolyError -> Bool     # V-C2
def check_common_divisor(f, h, result, mut budget: Budget) raises PolyError -> Bool           # G-C1
def check_bezout(f, h, result, mut budget: Budget) raises PolyError -> Bool                   # G-C2
def check_coprimality(f, h, certificate, mut budget: Budget) raises PolyError -> Bool         # G-C3, recursively
def check_resultant_identity(f, g, result, mut budget: Budget) raises PolyError -> Bool       # E-C1
def check_cancellation(p, q, certificate, mut budget: Budget) raises PolyError -> Bool        # F-C1
```

Each checker recomputes only the relation its certificate states
(PG2), and has a negative test with a corrupted certificate.
