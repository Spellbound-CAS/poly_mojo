# P5 Real roots

Algorithm provenance and independent regression cases: [reference material](../contributing/reference.md#p5-real-roots).

## P5.1 Sturm sequences and root counts

```text
struct SturmSequence:  polynomials: List[DensePoly[QQ]], quotients: List[DensePoly[QQ]]   # Definition T-1
def sturm_sequence(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> SturmSequence
def sign_changes(sequence: SturmSequence, a: QQ, mut budget: Budget) raises PolyError -> Int
struct RootCountCertificate:  sequence: SturmSequence, a: QQ, b: QQ, count: Int
def count_roots(f: DensePoly[QQ], a: QQ, b: QQ, mut budget: Budget) raises PolyError -> RootCountCertificate
                                        # DOMAIN if f(a) f(b) = 0 or a ≥ b
def check_root_count(f, certificate, mut budget: Budget) raises PolyError -> Bool                              # T-C1
```

**Required checks:** counts agree with planted roots of generated products of
linear and irreducible quadratic factors, including repeated roots;
`check_root_count` rejects certificates with a corrupted relation or a
wrong count.

## P5.2 Isolation

```text
struct RootInterval:  low: QQ, high: QQ           # open, exactly one root inside
struct RootIsolation:  exact_roots: List[QQ], intervals: List[RootInterval]   # sorted, disjoint
def isolate_real_roots(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> RootIsolation
def refine(f: DensePoly[QQ], interval: RootInterval, width: QQ, mut budget: Budget) raises PolyError -> RootInterval
def sign_on_interval(f: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> Optional[Int]
                                        # +1 or −1 when f has no root in [low, high] (Theorem T-6)
```

Algorithm T-4 with dyadic endpoints (Proposition T-3 for the initial
bound). `refine` bisects an isolating interval until it is narrower than
`width`.

**Required checks:** every interval contains exactly one root (checked with
`count_roots`), intervals are disjoint and sorted, exact rational roots are
found exactly, and the number of roots equals the number of distinct real
roots of generated inputs with known roots.
