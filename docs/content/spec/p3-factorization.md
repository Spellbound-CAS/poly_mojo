# P3 Factorization

Algorithm provenance and independent regression cases: [reference material](../contributing/reference.md#p3-factorization).

## P3.1 Squarefree decomposition

```text
struct SquarefreeDecomposition[C]:  unit: C, factors: List[(DensePoly[C], Int)]   # a_i with exponent i, a_i ≠ 1, by increasing i
def squarefree[C: CoefficientField](f: DensePoly[C], mut budget: Budget) raises PolyError -> SquarefreeDecomposition[C]
def check_squarefree_product(f, decomposition, mut budget: Budget) raises PolyError -> Bool          # Z-C1
```

Yun's algorithm (Theorem Z-2) handles characteristic-zero fields. Over
`ModInt`, derivative GCDs separate ordinary multiplicities, and repeated
inverse Frobenius steps recover factors whose multiplicities are divisible
by the characteristic.

**Required checks:** the product identity holds; each $a_i$ is coprime to its
derivative and to every other $a_j$ (checked by `gcd`); generated inputs
with planted multiplicities are recovered exactly.

## P3.2 Factorization over a prime field

```text
def factor_mod_p(f: DensePoly[ModInt], mut budget: Budget, seed: UInt64) raises PolyError -> List[DensePoly[ModInt]]
                                        # f squarefree and monic; irreducible monic factors, sorted
```

Distinct-degree factorization, then Cantor–Zassenhaus equal-degree
splitting with SplitMix64 seeded by `seed` (Theorem Z-3).
The output is sorted by degree, then by the coefficients from the leading
one down, compared as integers in $[0, p)$ (SymPy's order),
so it does not depend on the seed.

**Required checks:** products of generated irreducibles (found by exhaustive
search for small degrees) are factored back exactly; the output is the
same for several seeds.

## P3.3 Factorization over the integers and rationals

```text
struct Factorization[C]:  unit: C, factors: List[(DensePoly[C], Int)]
def factor(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> Factorization[QQ]
def factor_integer(f: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Factorization[ZZ]
def check_factorization(f, factorization, mut budget: Budget) raises PolyError -> Bool                   # Z-C2
```

Algorithm Z-6, with Hensel lifting (Theorem Z-4) to the bound of Theorem
Z-5, and subset recombination by trial division. Factors over
$\mathbb Q$ are monic and the unit carries the leading coefficient; over
$\mathbb Z$ they are primitive with positive leading coefficient, and the
unit is the content with its sign. The factors are sorted by degree, then
by the coefficients from the leading one down, compared as rationals over
$\mathbb Q$ and as integers over $\mathbb Z$.

**Required checks:** classical examples factor correctly (for instance
$x^4 + 1$ is irreducible over $\mathbb Q$ though reducible modulo every
prime; $x^{12} - 1$ has the cyclotomic factors); products of generated
irreducible factors are recovered; Swinnerton-Dyer polynomials of small
order factor (slowly, within a measured budget); every result satisfies
`check_factorization`.


## Additional exact fields and sparse multiplicity blocks

The original `factor_mod_p` contract remains the monic-squarefree splitting
primitive. `factor_finite_field` accepts arbitrary prime-field input and returns
`Factorization[ModInt]`, including the leading scalar and multiplicities.
`squarefree` also accepts ModInt, extracting inverse Frobenius roots when the
derivative vanishes. Zero and constants retain the existing empty-factor-list
convention.

`factor_gaussian` returns `Factorization[QQi]`. Its factors are monic, ordered by
degree and leading-to-constant coefficients with real components compared first.
The corresponding sparse producers are `factor_multivariate_gaussian` and
`factor_multivariate_finite_field`; they preserve the input order and use the
existing `MultivariateFactorization` certificate.

`squarefree_multivariate` accepts a field-valued sparse polynomial and returns
`SparseSquarefreeDecomposition`: a scalar and monic multiplicity blocks sorted
by multiplicity. `check_sparse_squarefree_product` verifies positive
multiplicities and exact reconstruction. It does not claim squarefreeness or
pairwise coprimality. All existing certificate fields and public declarations
remain unchanged. See the [algorithm notes](../theory/factorization.md).
