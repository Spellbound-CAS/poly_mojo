# P4 Gröbner bases

Algorithm provenance and independent regression cases: [reference material](../contributing/reference.md#p4-grobner-bases).

## P4.1 Bases and reduction

The built-in exact fields are QQ, QQi ($\mathbb Q(i)$), and prime ModInt.
QQi uses apn_mojo ExactComplex values with rational real and imaginary parts.
No approximate zero test or conversion from numerical Complex is performed.

```text
struct GroebnerBasis[C]:
    order: MonomialOrder
    basis: List[SparsePoly[C]]                    # reduced (Theorem B-4)
    cofactors: Optional[List[List[SparsePoly[C]]]] # per basis element, its combination of the inputs (Proposition B-6)
def groebner[C: CoefficientField](generators: List[SparsePoly[C]], order: MonomialOrder,
                                  mut budget: Budget, track_cofactors: Bool = False) raises PolyError -> GroebnerBasis[C]
def reduce[C: CoefficientField](f: SparsePoly[C], basis: GroebnerBasis[C], mut budget: Budget) raises PolyError -> ListDivision[C]
```

Buchberger's algorithm with the Gebauer–Möller criteria (Theorem B-3),
selecting S-pairs by the normal strategy (least LCM of leading monomials
first, ties by the pair's indices), then inter-reduction to the reduced
basis. All choices are deterministic (G1).

**Required checks:** every returned basis passes Buchberger's criterion (Theorem
B-2), checked independently by recomputing all S-polynomial remainders;
the reduced basis does not depend on the order of the generators; with
cofactor tracking, each basis element equals its combination of the
inputs.

## P4.2 Membership certificates

```text
struct MembershipCertificate[C]:  cofactors: List[SparsePoly[C]]     # f = Σ h_k f_k
def membership[C: CoefficientField](f: SparsePoly[C], generators: List[SparsePoly[C]], order: MonomialOrder,
                                    mut budget: Budget) raises PolyError -> Optional[MembershipCertificate[C]]
def check_membership(f, generators, certificate, mut budget: Budget) raises PolyError -> Bool                  # B-C1
```

`membership` computes a basis with cofactor tracking, reduces $f$, and
returns the cofactors of Proposition B-6 when the remainder is $0$.

**Required checks:** for generated $f = \sum h_k f_k$, a certificate is found
and checks; for $f$ outside the ideal (remainder nonzero), none is
returned; the infeasibility case $1 \in I$ (Theorem B-7) is a test.
