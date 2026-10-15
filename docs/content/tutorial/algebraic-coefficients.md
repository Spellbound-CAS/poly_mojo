# Calculate with algebraic coefficients

Use exact roots as ordinary polynomial coefficients. A polynomial containing
$\sqrt{2}$ stays exact through multiplication, division, interpolation, and
Gröbner reduction. You don't need to choose a number field or a common primitive
element before combining values.

## Build a polynomial from roots

`dense.root_between([-2, 0, 1], 1, 2)` selects $\sqrt{2}$. Put the returned
`algebraic.Real` in a coefficient sequence to infer `AA`, the field of real
algebraic numbers. Mixed scalar types use `poly.coefficients((...))`; sparse
construction accepts mixed `(power, value)` tuples directly.

<!-- example: docs/examples/array_algebraic_coefficients.mojo -->

The first multiplication produces $x^2-2$ exactly. Its domain remains `AA`,
even though its resulting coefficients happen to be integers. Domain inference
uses input types; arithmetic does not silently narrow the result's scalar type.
Empty arrays and zero polynomials also retain their domain.

## Include exact complex values

`algebraic.ascomplex(real, imag)` combines two exact real values. Either part can
be an integer, a rational, or a `Real`, including values from different algebraic
extensions. The resulting `algebraic.Complex` infers `QQbar`, the field of complex
algebraic numbers. This type is separate from APN's rounded `Complex`.

Use `algebraic.add`, `subtract`, `multiply`, `divide`, `negative`, `reciprocal`,
and `conjugate` for scalar arithmetic. `value.real()` and `value.imag()` return
exact `Real` components. Equality compares exact values: different root witnesses
for the same number compare equal. Complex values have no ordering.

| Inputs combined in a calculation | Inferred domain |
|---|---|
| `Real` with integers or rationals | `AA` |
| `Real` with APN `ExactComplex` | `QQbar` |
| `algebraic.Complex` with any characteristic-zero exact scalar | `QQbar` |

Finite-field values have no implicit algebraic conversion. Approximation also
stays explicit: pass a precision `context=` to `dense.asarray` or
`sparse.asarray` to convert a whole array. Scalar `bounds`, `refine`, and
`approximate` remain available when you need rational intervals or values.
[Certified evaluation](array-certified.md) retains the conversion witnesses.

## Reuse the polynomial operations

The ordinary dense and sparse functions accept these coefficients: arithmetic,
calculus, evaluation, composition and substitution, division, GCDs, resultants,
squarefree decomposition, and interpolation use the inferred field. Sparse
Gröbner bases and membership use it too. Checked variants retain their usual
certificates and claims; for example, `exact.check(dense.divide_checked(...))`
verifies the division identity. Formal series and alternative bases use the
same field adapters.

`sparse.specialize` accepts `(axis, value)` assignments. With algebraic or complex
coefficient arrays it returns another sparse array, retaining the original power
shape. With integer or rational coefficients and real assignments, it preserves
the existing `algebraic.Specialization` result and its group evaluation proofs.
Pass that result to `sparse.asarray` to continue ordinary polynomial arithmetic,
or use `algebraic.coefficients` to extract one remaining coefficient axis.

`dense.real_roots` accepts real algebraic coefficient arrays and retains its
completeness certificate. `sparse.roots_at` retains its integer/rational-input contract. For algebraic
coefficients, use `algebraic.coefficients(sparse.specialize(...), axis=...)`
and pass that sequence to `dense.real_roots`. `dense.complex_roots` and `sparse.complex_roots` accept both `AA` and `QQbar`
coefficients and return exact complex values with multiplicities. See the
[complex-root workflow](array-roots.md#find-all-exact-complex-roots) for ordering,
verification, and sparse axis selection. Complex root-cover APIs return
certified numerical regions.

## Factor over real or complex algebraic numbers

`dense.factor(p)` and `sparse.factor(p)` return the familiar `(scalar, factors)`
pair. Each factor has a positive multiplicity, and their powered product times
the scalar equals `p`. For a univariate polynomial, `AA` factors are linear or
quadratic, and `QQbar` factors are linear. Sparse multivariate polynomials can
have irreducible factors of higher degree in either field.

The input field determines the answer. Integer coefficients for $x^2-2$ keep
it irreducible over `ZZ`. Including `algebraic.asreal(1)` requests `AA`, where
it splits into $(x-\sqrt2)(x+\sqrt2)$. To split $x^2+1$, choose `QQbar` by
including `algebraic.ascomplex(1)`.

<!-- example: docs/examples/array_algebraic_factor.mojo -->

Univariate factors are monic; the scalar is the input's leading coefficient. They sort by
increasing degree, then by coefficients read from leading to constant. Complex
coefficients use real part first, then imaginary part for this deterministic
order. Zero and constants have no factors, with their value returned as the
scalar. Neither factors nor zero results narrow their coefficient domain.

`factor_checked` retains the original polynomial and factorization. Over these
two fields, `exact.check` verifies the entire factorization: coefficient
witnesses, monicity, distinct sorted factors, positive multiplicities, and the
product. For univariate inputs, it checks irreducibility directly: a linear
factor is irreducible, and an `AA` quadratic must have strictly negative
discriminant. Checking does not solve for the polynomial's roots again.

### Factor several variables

Use tuple powers for multivariate input. The same functions infer the field,
retain unused axes, and return sparse factors. For example, $x^2+y^2$ stays
irreducible over `AA` but splits into $(x-iy)(x+iy)$ over `QQbar`:

<!-- example: docs/examples/array_algebraic_multivariate_factor.mojo -->

Multivariate factors sort by total degree, then their stored monomials and
coefficients. Each factor's greatest stored coefficient is one; the returned
scalar accounts for that normalization. Use the returned scalar when
reconstructing the input. A complex coefficient field does not force every
multivariate factor to be linear: $x^2+y$ is irreducible over `QQbar`.

The multivariate checker uses retained univariate image factorizations to
verify irreducibility. It checks those certificates and rules out every proper
divisor combination without rerunning the factorization producer. Keep the
complete checked result for `exact.check`: `result.evidence()` exposes the
core product certificate, whose own checker establishes reconstruction only.

Set `max_steps` for large inputs. The algorithm removes coordinate factors and
repeated factors first, but general recombination builds a dense image whose
storage grows with the product of the per-variable degree ranges. The search
and its independent check can take exponential work. Budget exhaustion raises
without returning partial factors. Root isolation and algebraic arithmetic
can also dominate work for small inputs.

## Costs, evidence, and boundaries

Every computation accepts a shared `Budget` or `max_steps=`. Coefficient-field
arithmetic also verifies each scalar certificate before using its result in
polynomial operations or checkers. A small coefficient
array can still require substantial work: algebraic arithmetic constructs and
factors annihilating polynomials, then selects the intended root. Degrees and
coefficient sizes can grow when combining unrelated values. Budget exhaustion
raises `LIMIT` and leaves borrowed inputs unchanged.

`algebraic.check` verifies a real witness and its retained arithmetic or evaluation
step; for a complex value it checks both components. It does not certify an entire
history of scalar calls. Polynomial checked results retain their separate,
operation-specific evidence.

[Rational functions](array-rational.md) accept AA and QQbar coefficients.
[CAD workflows](array-cad-workflows.md) accept AA coefficients directly and keep
user coordinates and signs in their original order. To save exact algebraic
arrays and their witnesses, see [serialization](array-serialization.md).
See the [algebraic reference](../reference/array-algebraic.md) for scalar helpers.


For ratios of these polynomials, use [`poly.rational`](array-rational.md).
Construction, arithmetic, substitution, evaluation, partial fractions, Laurent
expansions, and residues retain the inferred exact field. Algebraic evaluation
points and expansion centers can also promote an existing QQ function.
