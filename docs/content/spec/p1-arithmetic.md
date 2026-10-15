# P1 Representations and arithmetic

Algorithm provenance and independent regression cases: [reference material](../contributing/reference.md#p1-representations-and-arithmetic).

All computational calls pass a shared explicit budget.
Fixed-signature equality/copy/hash/format protocols are exempt; algorithms
precharge their use. Budget exhaustion does not modify polynomial operands.

## P1.1 Errors and budgets

```text
struct PolyErrorKind:  LIMIT | CAPACITY | DOMAIN | UNSUPPORTED | ENCODING     # ImplicitlyCopyable
struct PolyError(kind: PolyErrorKind, message: String)                        # Writable
struct Budget:
    def __init__(out self, steps: UInt64 = UInt64.MAX)
    def charge(mut self, n: UInt64) raises PolyError        # LIMIT when it would go negative
    def remaining(self) -> UInt64
```

**Required checks:** a function that exhausts its budget raises `LIMIT` and
leaves its inputs unchanged; apn_mojo errors are converted at each call
site (a compile-fail test, as in the kernel, shows that a
`raises PolyError` function cannot call apn_mojo unconverted).

## P1.2 Coefficient rings

```text
trait CoefficientRing(Copyable, Equatable, Writable, Deinitable):
    comptime Context: Copyable & Equatable & Writable & Deinitable   # the ring's description
    @staticmethod
    def zero(ctx: Self.Context) -> Self
    @staticmethod
    def one(ctx: Self.Context) -> Self
    def context(self) -> Self.Context      # the ring this element belongs to
    def is_zero(self) -> Bool
    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self
    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self
    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self
    def neg(self, mut budget: Budget) raises PolyError -> Self
    def size(self) -> Int                  # bytes of the value, for budgets

trait CoefficientField(CoefficientRing):
    def inverse(self, mut budget: Budget) raises PolyError -> Self     # DOMAIN on zero

struct NoContext(...):  (empty; all values equal)            # the context of ZZ, QQ and QQi
struct Modulus(...):    value: UInt64                          # the context of ModInt: a prime < 2^62

struct ZZ(CoefficientRing):      value: Integer;                  Context = NoContext
struct QQ(CoefficientField):     value: Rational;                 Context = NoContext
struct QQi(CoefficientField):    value: ExactComplex;             Context = NoContext
struct ModInt(CoefficientField): value: UInt64, modulus: UInt64;  Context = Modulus
```

The constants take the ring's context: `ModInt.zero(ctx)` and
`ModInt.one(ctx)` get their modulus from `ctx`, never from global or
thread-local state. `ModInt` operations on different moduli raise `DOMAIN`.
Products modulo $p$ use 128-bit intermediate arithmetic or Montgomery
multiplication.

`QQi` lives in `exact_complex.mojo`, backed by apn_mojo `ExactComplex`.
Each component is a canonical rational. Arithmetic delegates its rational
steps to the existing budgeted QQ methods, including intermediate products
and the squared norm used in inversion. The sum of component `size()` values
accounts for copies and comparisons. Inverting zero raises `DOMAIN` before
arithmetic; exhaustion during any component operation leaves both inputs
unchanged. Its principal consumer is the existing generic Gröbner API.

**Required checks:** the ring axioms hold on generated values for each type
(property tests); `ModInt` agrees with `Integer` arithmetic reduced modulo
$p$; inverses satisfy $a \cdot a^{-1} = 1$; `zero(ctx)`, `one(ctx)`, and
`context()` agree (`x.context() == ctx` for constants made from `ctx`);
values of two moduli, interleaved in one computation, each stay correct,
and mixing them raises `DOMAIN`.

## P1.3 Monomials and orders

```text
struct Monomial(Copyable, Movable, Equatable, Hashable, Writable):
    pairs: List[(Int64, UInt32)]        # indices strictly increasing, exponents ≥ 1
    def degree(self, mut budget: Budget) raises PolyError -> UInt64          # total degree
    def mul(self, other: Monomial, mut budget: Budget) raises PolyError -> Monomial     # CAPACITY on exponent overflow
    def divides(self, other: Monomial, mut budget: Budget) raises PolyError -> Bool
    def div(self, other: Monomial, mut budget: Budget) raises PolyError -> Monomial     # DOMAIN unless divides
    def lcm(self, other: Monomial, mut budget: Budget) raises PolyError -> Monomial

struct MonomialOrder(Copyable, Movable, Equatable):
    kind: KERNEL | LEX | GRLEX | GREVLEX
    variables: List[Int64]              # empty for KERNEL; the order x1 > x2 > … otherwise
    def compare(self, a: Monomial, b: Monomial, mut budget: Budget) raises PolyError -> Int   # −1, 0, 1
```

Comparing under `LEX`, `GRLEX`, or `GREVLEX` a monomial with a variable
not in the list raises `DOMAIN`.

**Required checks:** `KERNEL` matches Definition R-3 on the kernel reference's
vectors; the three monomial orders satisfy Definition R-2 on generated
monomials (totality, multiplicativity, well-foundedness on bounded
degree); `KERNEL`'s counterexamples of Proposition R-1 are tests.

## P1.4 Sparse polynomials

```text
struct SparsePoly[C: CoefficientRing](Copyable, Movable, Equatable, Writable):
    order: MonomialOrder
    ring: C.Context                     # the coefficients' ring
    terms: List[Term[C]]                # canonical (R-R1), increasing in order

    @staticmethod
    def zero(order: MonomialOrder, ring: C.Context, mut budget: Budget) raises PolyError -> Self
    @staticmethod
    def constant(c: C, order: MonomialOrder, mut budget: Budget) raises PolyError -> Self  # ring from c.context()
    @staticmethod
    def variable(index: Int64, order: MonomialOrder, ring: C.Context, mut budget: Budget) raises PolyError -> Self
    @staticmethod
    def from_terms(terms: List[Term[C]], order: MonomialOrder, ring: C.Context, mut budget: Budget) raises PolyError -> Self
                                        # sorts, combines equal monomials, drops zeros;
                                        # DOMAIN if a coefficient's context is not ring
    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self
    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self
    def neg(self, mut budget: Budget) raises PolyError -> Self
    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self
    def scale(self, c: C, mut budget: Budget) raises PolyError -> Self
    def pow(self, n: UInt64, mut budget: Budget) raises PolyError -> Self      # repeated squaring
    def is_zero(self) -> Bool
    def leading(self, mut budget: Budget) raises PolyError -> Term[C]                              # DOMAIN on zero
    def total_degree(self, mut budget: Budget) raises PolyError -> Int                                               # −1 for zero
    def degree_in(self, index: Int64, mut budget: Budget) raises PolyError -> Int
    def variables(self, mut budget: Budget) raises PolyError -> List[Int64]
    def reorder(self, order: MonomialOrder, mut budget: Budget) raises PolyError -> Self
    def check_invariants(self, mut budget: Budget) raises PolyError
```

Operands of `add`, `sub`, and `mul` must have the same order and the same
ring (`DOMAIN` otherwise); results keep both, and `pow(0)` is `one(ring)`
even for the zero polynomial. For LEX, GRLEX and GREVLEX, multiplication merges
sorted product rows with a binary heap at 64 or more term pairs, retaining
at most one pending product per term of the shorter input. Smaller products
retain the direct term-product and sorting path. Single-term multiplication
preserves its ordered fast path.

Eligible KERNEL products also stream from 256 pairs: a common single
variable, or fixed monomial supports with an overlapping variable. Those cases
preserve row order under translation. KERNEL itself is still not a general
multiplicative order. Large mixed-support products may instead use an internal
LEX order and restore KERNEL order afterward; disjoint products keep the direct
path. Every path charges arithmetic, traversal and allocation to the budget.

**Required checks:** the commutative-ring axioms hold on generated polynomials;
`pow` agrees with repeated `mul`; every operation preserves R-R1; results
agree with the `add`, `mul`, `neg`, and `pow` vectors of
`tests/reference/kernel/poly_ops.txt`, which come from the kernel's
`poly_reference.py`.

## P1.5 Dense univariate polynomials

```text
struct DensePoly[C: CoefficientRing](Copyable, Movable, Equatable, Writable):
    variable: Int64
    ring: C.Context                     # the coefficients' ring
    coefficients: CoefficientBuffer[C]  # c_0 … c_d, c_d ≠ 0; empty for zero
    @staticmethod
    def from_coefficients(variable: Int64, ring: C.Context, var coefficients: List[C], mut budget: Budget) raises PolyError -> Self
    @staticmethod
    def from_coefficients(variable: Int64, ring: C.Context, var coefficients: CoefficientBuffer[C], mut budget: Budget) raises PolyError -> Self
    def degree(self) -> Int
    def add, sub, neg, scale, mul, pow  # as for SparsePoly; mul uses Karatsuba above a threshold
    def to_sparse(self, order: MonomialOrder, mut budget: Budget) raises PolyError -> SparsePoly[C]
    @staticmethod
    def from_sparse(p: SparsePoly[C], mut budget: Budget, *, variable: Optional[Int64] = None) raises PolyError -> Self   # DOMAIN if p has another variable; keeps p.ring
```

The optional `variable` preserves the variable of zero and constants during
a round trip; without it, a nonconstant polynomial supplies its sole variable
and constants use index zero. More than one variable is `DOMAIN`.
`from_coefficients` validates coefficient contexts and removes trailing zeros.
It adopts either owned storage form without flattening. `CoefficientBuffer`
keeps contiguous lists for known-size results and fixed-size chunks for
unknown-length construction. Indexing borrows entries; internal arithmetic reads
both layouts without an adapter copy. An explicit `to_list(budget)` copies into
contiguous storage. See [dense coefficient storage](../reference/contract-p.md#dense-coefficient-storage).

The measured Karatsuba threshold is 32 coefficients for `ZZ` and `QQ`,
64 for `ModInt`; `DensePoly[C].karatsuba_threshold()` reports it.
`mul_schoolbook` exposes the direct convolution for cross-checks and
measurements.

**Required checks:** Karatsuba agrees with schoolbook multiplication on generated
inputs of every degree up to twice the threshold; conversions round-trip.

## P1.6 Evaluation, substitution, derivatives

```text
def evaluate[C](p: SparsePoly[C], point: List[(Int64, C)], mut budget: Budget) raises PolyError -> C
def evaluate_dense[C](p: DensePoly[C], x: C, mut budget: Budget) raises PolyError -> C      # Horner
def substitute[C](p: SparsePoly[C], index: Int64, q: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]
def derivative[C](p: SparsePoly[C], index: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]
```

`evaluate` raises `DOMAIN` if a variable of `p` has no value.

**Required checks:** evaluation is a ring homomorphism on generated inputs
(Proposition R-5); the product rule holds for `derivative`.

## P1.7 The kernel codec

```text
def to_kernel_bytes(p: SparsePoly[QQ], mut budget: Budget) raises PolyError -> List[UInt8]
                                       # ENCODING for 0 and 1; re-sorts into KERNEL order if needed
def from_kernel_bytes(data: Span[UInt8, _], mut budget: Budget) raises PolyError -> SparsePoly[QQ]
                                       # in KERNEL order; ENCODING on anything not canonical
```

The kernel ceiling of 1,048,576 bytes also applies (`CAPACITY`); oversized
counts and truncated fields are rejected before input-sized allocation.

**Required checks:** both reproduce the kernel reference's vectors byte for
byte (Theorem R-6: the `enc` lines of `tests/reference/kernel/poly_ops.txt`
and the table in `poly_vectors.txt`), reject with `ENCODING` every `bad`
line of `poly_ops.txt`, and round-trip generated polynomials.
