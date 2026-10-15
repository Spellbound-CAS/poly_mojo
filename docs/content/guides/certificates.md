# Understand what a certificate proves

A certificate is evidence that a separate checker can verify. The producer
may use expensive search, modular reconstruction, or adaptive numerical work;
the checker verifies a stated relation against the original inputs.

A successful check proves that relation. It does not automatically prove every
property suggested by the producer's name.

## Match the check to your claim

| Evidence | Established claim | Additional claim it does not establish by itself |
|---|---|---|
| Division result | The division identity and the checker's remainder conditions | Properties of a different divisor list or order |
| Common-divisor cofactors | The proposed divisor divides both inputs | That it is a greatest common divisor |
| Bézout or coprimality evidence | The corresponding identity or coprimality claim | A relation for inputs other than those checked |
| Core cancellation certificate | The original guarded fraction identity | Equality at removed poles or maximal cancellation |
| `rational.check(function)` | Cancellation identity, coprimality, and monic denominator over QQ, QQi, AA, or QQbar | Equality at removed poles |
| `rational.check(partial_fractions)` | Normalized source, irreducible factors with multiplicities, proper terms, and recomposition | An analytic bound for a separate Laurent truncation |
| Resultant identity | The result lies in the input ideal | Equality to the Sylvester determinant |
| Core factorization product | Unit and factor powers multiply back to the input | Irreducibility of every factor |
| `exact.check` for `AA`/`QQbar` factorization | Complete irreducible factorization with normalized factors and multiplicities | The core product certificate alone does not include retained multivariate image proofs |
| Membership certificate | The target is a sum of multiples of the original generators | That a supplied list is a Gröbner basis |
| Sturm/root evidence | The specified root counts and isolating conditions | Multiplicity unless the API explicitly includes it |
| Algebraic arithmetic | The arithmetic relation and selection of the result root | Minimality of a supplied defining polynomial |
| Enclosure trace | A bound valid for the represented inputs | The narrowest enclosure or a requested accuracy target |
| Complete numerical root cover | Global coverage and the reported count/width/isolation claims | Isolation or accuracy when their result fields are false |
| CAD evidence | The decomposition or covering claim of that API | A full decomposition from a selective covering alone |

For example, a gcd producer supplies a candidate $g$ and cofactors with
$a=ga'$ and $b=gb'$. These two products establish divisibility. To conclude
maximality, also use the applicable Bézout or coprimality evidence. The
[division and gcd tutorial](../tutorial/division.md) follows this distinction.

## Keep the inputs with the evidence

Pass the original polynomial, domain context, order, and relevant options to
the checker specified by the API. Generator order matters: membership
multipliers refer to the original list, including zero generators. Numerical
options matter too, because width and isolation requirements form part of
the claim.

Treat a false checker result as rejected evidence. A checker can also raise
when it cannot complete within the budget; that is not acceptance. For
numerical results, check the evidence and then require the outcome fields
needed by your application. A valid enclosure can be too wide for your task.

## Independence and trust

Checkers avoid repeating the producer's search. Exact identity checks use
polynomial arithmetic. Numerical checkers reconstruct the evaluation graph
and verify rational endpoint bounds independently of the floating-point
producer. CAD projection checking verifies its determinant obligations,
which are stronger than a generic resultant ideal-membership identity.

Independent verification need not be cheap. Rational partial-fraction checks
can use a separate exhaustive irreducibility search when prime reduction
cannot settle the claim. Keep the supplied budget explicit when checking
untrusted or large results; see the [rational-function tutorial](../tutorial/array-rational.md#read-the-decomposition-and-its-evidence).

The library's checkers are useful within ordinary applications. A separate
proof system still needs its own acceptance path; a library certificate does
not automatically become a kernel proof. The exact kernel codec is a specific
integration interface and does not encode numerical enclosures.

## Incremental certificates

Streaming checkers can consume evidence as it arrives. Acceptance of an
individual event does not establish that the stream is complete. Call `finish`
with the producer's receipt and accept the result only if final validation
succeeds. See [interchange and streaming](interchange.md) for a complete example.
