# Isolate real roots exactly

For coefficient-array inputs, start with the [array workflow](array-roots.md).

Root isolation locates each distinct real root in an interval without
replacing it by an unchecked floating-point approximation. With rational
coefficients, the library can also recover rational roots exactly.

## Count and isolate

For $p(x)=(x^2-2)(3x-1)$, the real roots are $-\sqrt{2}$, $1/3$, and
$\sqrt{2}$. The following program first checks a root count, then separates
the rational root from the two irrational-root intervals.

<!-- example: docs/examples/tutorial_roots.mojo -->

`count_roots` takes ordered rational endpoints that are not roots. Its
certificate contains the evidence for the count; `check_root_count` verifies
the remainder identities and sign variations without repeating isolation.

`isolate_real_roots` returns `exact_roots` and `intervals`. The first list
contains rational roots, including non-dyadic values such as $1/3$. The
second list gives isolating intervals for the remaining real roots. Read
both lists when collecting the answer.

## Understand multiplicity and boundaries

The exact Sturm interface counts **distinct** real roots. A repeated root
counts once. Nonzero constants have no roots; the zero polynomial is rejected
because it does not have a finite root set.

Isolation uses rational arithmetic and returns intervals with nonroot dyadic
endpoints. The default width is below $2^{-40}$; rational-root reconstruction
can require further narrowing. `refine` takes a positive rational target
width to narrow an interval further.

`sign_on_interval` returns no sign if any root lies in the closed interval.
A missing sign is a mathematical outcome, distinct from a work-limit error.

## Choose the next interface

Use these exact APIs when rational roots and exact isolating intervals are
the result you need. Use [real algebraic values](algebraic-numbers.md) when
you want to compute with a selected root as a number.

For configurable accuracy, complex roots, uncertain coefficients, or a root
cover with explicit completeness and isolation flags, use the
[numerical root workflow](numerical-roots.md). It also supports refining
previous evidence and adding multiplicity certificates for exact real roots.

See the [root reference](../reference/p5.md) and
[Sturm theory](../theory/roots.md) for the exact endpoint conventions.
