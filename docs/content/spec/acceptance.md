# Validation and compatibility

Validation combines algorithm tests, independent expected values, corrupted
certificates, reproducibility, compiled consumers, and compatibility snapshots.
The [maintainer guide](../contributing/handoff.md) gives current commands.

## Tests by phase

| Phase | Tests |
|---|---|
| P1 | ring axioms and invariants for every coefficient type and representation; orders; Karatsuba against schoolbook; evaluation as a homomorphism; the codec's vectors and rejections |
| P2 | division identities; gcd recovery of planted common factors; modular and subresultant gcds agree; an unlucky prime; resultants and the Sylvester determinant; coprimality certificates; normalized fractions |
| P3 | squarefree decompositions with planted multiplicities; factorization over $\mathbb F_p$ independent of the seed; classical factorizations over $\mathbb Q$ ($x^4 + 1$, cyclotomic polynomials, small Swinnerton-Dyer polynomials) |
| P4 | Buchberger's criterion on every returned basis; independence of generator order; membership certificates found and refused; $1 \in I$ |
| P5 | root counts and isolation against planted roots, repeated roots included |
| all | every `check_*` function accepts its algorithm's certificates on generated inputs and rejects corrupted ones |

## Cross-checks

| Check | How |
|---|---|
| Independent arithmetic | P1's sums, products, negations, and powers agree with the kernel's `poly_reference.py` on the vectors of `tests/reference/kernel/poly_ops.txt` (generated from a copy of the script, `tests/reference/kernel/`) |
| Determinism | every algorithm, run twice and with the polynomials' terms supplied in different orders, gives identical results and certificates |
| The reference corpus | every case of `tests/reference/corpus.txt` for the operations of a phase passes in that phase's tests, as the [reference material](../contributing/reference.md#using-it-in-tests) page describes; `pixi run corpus-check` keeps the corpus itself consistent |
