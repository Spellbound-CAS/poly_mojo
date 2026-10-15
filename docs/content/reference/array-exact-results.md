# Checked exact results

Use `dense.gcd_checked`, `sparse.divide_checked`, and the other checked producers
when you need both a polynomial answer and evidence for an independent check.
They accept the same coefficient sequences and sparse terms as ordinary
operations. Domains are inferred, and integer division and interpolation promote
to rational coefficients.

Import `exact` from `poly` to inspect and verify a result. `exact.value(result)`
returns a GCD, resultant, or interpolant. Division has `quotient`, `quotients`,
and `remainder`; factorization has `unit` and `factors`. Every polynomial
accessor returns independent storage in the original dense or sparse shape.
A dense resultant is a constant array, with the canonical empty array for zero.
`resultant_certified` retains the determinant transcript for a stronger check;
`resultant_checked` retains the existing ideal identity. Sparse checked GCDs
accept `method="modular"`; checked division accepts `order=` and `axes=`.

`factor_checked` supports integers, rationals, Gaussian rationals and prime-field
coefficients, including repeated factors. It also supports univariate and
multivariate `AA` and `QQbar` inputs. `squarefree_checked`
supports all these domains, promoting integers to rationals. Dense inputs are univariate;
sparse inputs can have any number of axes. One-axis sparse squarefree results
retain the existing dense evidence format; multiaxis results retain sparse
blocks. Functional extraction and checking work the same way for both.

## What a successful check establishes

| Producer | Claim verified by `exact.check` |
|---|---|
| `gcd_checked` | A common divisor with coprime cofactors, hence a GCD up to a unit; also handles two zero inputs |
| `coprimality` | The two inputs have no common nonunit factor |
| Dense `divide_checked` | The quotient/remainder identity and degree bound |
| Sparse `divide_checked`, `divide_many_checked` | The ordered quotient/remainder identity; no claim of a unique or reduced remainder |
| `resultant_checked` | An identity expressing the returned polynomial as a combination of the inputs; not the determinant definition of the resultant |
| `resultant_certified` | The full Sylvester determinant, over QQ, QQi, GF(p), AA, and QQbar; integer inputs promote to QQ |
| `factor_checked` over `AA` or `QQbar` | Complete normalized irreducible factorization, including positive multiplicities and strict factor ordering |
| `factor_checked` over other domains | Unit times the powered factors equals the input; no irreducibility claim |
| `squarefree_checked` | The same product identity; no squarefreeness claim |
| `interpolate_checked`, `hermite_interpolate_checked` | Distinct nodes, the degree bound, and every supplied value or derivative constraint |

`coprimality` returns `None` when the inputs have a common factor. A returned
certificate still needs checking. `algebraic.compare_checked` and
`algebraic.sign_checked` use `algebraic.check` for their exact ordering evidence.

For algebraic factorization, the facade checker validates coefficient witnesses,
monicity, and reconstruction. Univariate irreducibility follows from degree and
negative real quadratic discriminants. Multivariate results retain auxiliary
image factorizations, which the checker validates before exhausting candidate
proper divisors. It never calls the factorization producer, but the exhaustive
check can require exponential work. Keep the complete `Result` for this claim:
`evidence()` exposes the core product certificate. The underlying core
`check_factorization` and `check_multivariate_factorization` keep their
product-only contracts.

Results retain their original operands and share storage when copied.
`evidence()` borrows an immutable core witness. An explicit certificate copy is
useful for interoperability; its allocation is outside the facade's work budget.
Result construction, functional accessors, and checks accept either a shared
`Budget` or `max_steps=`. A limit raises an error, never an unchecked answer.
See the [walkthrough](../tutorial/array-checked.md) for a complete example.

<!-- api: src/poly/exact.mojo Result -->

<!-- api: src/poly/exact.mojo check -->

<!-- api: src/poly/exact.mojo value -->

<!-- api: src/poly/exact.mojo quotient -->

<!-- api: src/poly/exact.mojo quotients -->

<!-- api: src/poly/exact.mojo remainder -->

<!-- api: src/poly/exact.mojo cofactors -->

<!-- api: src/poly/exact.mojo unit -->

<!-- api: src/poly/exact.mojo factors -->

Save a complete result with `io.dumps(result)` and restore it with
`io.loads[type_of(result)](text)`. Loading verifies the retained operands and
certificate before returning. See [checked archives](../tutorial/array-serialization.md#save-a-checked-result).

<!-- api: src/poly/_determinant.mojo Certificate -->

<!-- api: src/poly/_determinant.mojo Step -->
