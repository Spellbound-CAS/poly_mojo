# Ground rules

These rules preserve the package's correctness and compatibility. The G rules
are general engineering requirements; the PG rules are specific to
polynomials.

| Rule | Requirement |
|---|---|
| G1 | Deterministic results and certificates |
| G2 | Explicit worklists instead of recursion over input-sized structures |
| G3–G4 | Work budgets and checked capacities |
| G5 | No views retained across growth of their backing storage |
| G9 | Typed errors and unchanged borrowed inputs on failure |
| G11–G12 | Positive and negative tests with independently derived expectations |
| G13–G14 | Public docstrings and documentation updated with behavior |
| G16 | Unresolved computations cannot be reported as proved answers |

The package adds these polynomial-specific requirements:

| ID | Rule |
|---|---|
| PG1 | **Canonical forms always.** Every `SparsePoly` and `DensePoly` satisfies Requirement R-R1 at every public boundary; `check_invariants` runs in tests after every operation. |
| PG2 | **Certificates are complete.** A certificate holds everything its `check_*` function needs; checking never calls the algorithm that produced it. |
| PG3 | **Exact symbolic core; explicit numerical domains.** Exact algorithms use exact coefficients only. apn_mojo Float/Complex/Ball/ComplexBall appear only in dedicated numerical modules with explicit precision contexts. Numerical certificate checkers use exact rational bounds and distinguish point values from interval families. |
| PG4 | **Fixed primes and seeds.** Modular algorithms take primes from `primes.mojo`'s fixed sequence; randomized steps take a SplitMix64 seeded by the call. |
| PG5 | **Recursion over variables is bounded.** Multivariate algorithms that recurse over variables (gcd, resultants) recurse at most once per variable, which is bounded by the input; recursion over input-sized structures (terms, degrees) is iterative (G2). |
| PG6 | **The kernel's encoding is the contract.** `to_kernel_bytes` and `from_kernel_bytes` follow the kernel's canonical byte format exactly, and their vectors come from the kernel's `poly_reference.py` (copied, with its vectors, in `tests/reference/kernel/`). |
| PG7 | **Pure Mojo.** The package is Mojo over `std` and apn_mojo: no `external_call`, no dynamic loading, nothing from `std.ffi` but `_get_global`, no C or C++ library; ported code comes only from MIT, BSD, ISC, or Apache-2.0 sources and keeps their notices, never from GPL or LGPL code (`pixi run layering` checks the first part). |
