"""Exact coefficient rings with explicit contexts and budgeted arithmetic."""
from apn_mojo import Integer, Rational
from .errors import Budget, PolyError, PolyErrorKind, _sum, _product, _charge_size

@fieldwise_init
struct NoContext(ImplicitlyCopyable, Equatable, Writable):
    """The unique context of integers, rationals, and Gaussian rationals."""

    def __eq__(self, other: Self) -> Bool:
        """Compare the unique contexts.

        Args:
            other: The other context.

        Returns:
            Always True.
        """
        return True

trait CoefficientRing(Copyable, Equatable, Writable, Deinitable):
    """An exact commutative ring with explicit context and charged arithmetic."""
    comptime Context: Copyable & Equatable & Writable & Deinitable
    """The type describing the coefficient ring."""

    @staticmethod
    def zero(ctx: Self.Context) -> Self:
        """Construct the additive identity.

        Args:
            ctx: The ring context.

        Returns:
            Zero in that ring.
        """
        ...

    @staticmethod
    def one(ctx: Self.Context) -> Self:
        """Construct the multiplicative identity.

        Args:
            ctx: The ring context.

        Returns:
            One in that ring.
        """
        ...

    def context(self) -> Self.Context:
        """Report the ring context.

        Returns:
            This value's ring.
        """
        ...

    def is_zero(self) -> Bool:
        """Test the additive identity in constant time.

        Returns:
            Whether the value is zero.
        """
        ...

    def size(self) -> Int:
        """Report the magnitude storage size in bytes.

        Returns:
            A positive size used for precharging arithmetic.
        """
        ...

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact sum.

        Args:
            other: The other coefficient in the same ring.
            budget: The shared work allowance.

        Returns:
            The result in the same ring.

        Raises:
            PolyError: DOMAIN for different rings, LIMIT or CAPACITY for work limits.
        """
        ...

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact difference.

        Args:
            other: The other coefficient in the same ring.
            budget: The shared work allowance.

        Returns:
            The result in the same ring.

        Raises:
            PolyError: DOMAIN for different rings, LIMIT or CAPACITY for work limits.
        """
        ...

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact product.

        Args:
            other: The other coefficient in the same ring.
            budget: The shared work allowance.

        Returns:
            The result in the same ring.

        Raises:
            PolyError: DOMAIN for different rings, LIMIT or CAPACITY for work limits.
        """
        ...

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Compute the additive inverse.

        Args:
            budget: The shared work allowance.

        Returns:
            The negated coefficient.

        Raises:
            PolyError: LIMIT or CAPACITY for work limits.
        """
        ...

trait CoefficientField(CoefficientRing):
    """A coefficient ring with inverses for nonzero elements."""

    def inverse(self, mut budget: Budget) raises PolyError -> Self:
        """Compute the multiplicative inverse.

        Args:
            budget: The shared work allowance.

        Returns:
            The inverse in this field.

        Raises:
            PolyError: DOMAIN for zero, LIMIT or CAPACITY for work limits.
        """
        ...

def _work(a: Int, b: Int, quadratic: Bool, mut budget: Budget) raises PolyError:
    var n = _sum(a, b)
    if quadratic:
        n = _product(n, n)
    _charge_size(n, budget)

def _submul(mut r: ZZ, c: ZZ, h: ZZ, mut budget: Budget) raises PolyError:
    """r -= c h in place, without a temporary product, charged as r.sub(c.mul(h)) with the product's size taken as the factors' total."""
    _work(c.size(), h.size(), True, budget)
    _work(r.size(), _sum(c.size(), h.size()), False, budget)
    try:
        r.value._submul(c.value, h.value)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


@fieldwise_init
struct ZZ(CoefficientRing):
    """An exact integer coefficient."""
    var value: Integer
    """The canonical arbitrary-precision value."""
    comptime Context = NoContext
    """The unique ring context."""

    @staticmethod
    def zero(ctx: NoContext) -> Self:
        """Construct zero.

        Args:
            ctx: The ring context.

        Returns:
            The additive identity.
        """
        return Self(Integer(0))

    @staticmethod
    def one(ctx: NoContext) -> Self:
        """Construct one.

        Args:
            ctx: The ring context.

        Returns:
            The multiplicative identity.
        """
        return Self(Integer(1))

    def context(self) -> NoContext:
        """Report the unique context.

        Returns:
            The ring's context.
        """
        return NoContext()

    def is_zero(self) -> Bool:
        """Test whether this coefficient is zero.

        Returns:
            Whether the value is zero.
        """
        return self.value == Integer(0)

    def size(self) -> Int:
        """Report a positive magnitude size in bytes.

        Returns:
            The size used by the work model.
        """
        return self.value.magnitude_bit_length() // 8 + 1

    def __eq__(self, other: Self) -> Bool:
        """Compare exact values (the fixed-signature protocol exception).

        Args:
            other: The value to compare.

        Returns:
            Whether the canonical coefficients agree.
        """
        return self.value == other.value

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact sum after charging its size.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact result.

        Raises:
            PolyError: LIMIT or CAPACITY if the computation cannot proceed.
        """
        _work(self.size(), other.size(), False, budget)
        try:
            return Self(self.value + other.value)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact difference after charging its size.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact result.

        Raises:
            PolyError: LIMIT or CAPACITY if the computation cannot proceed.
        """
        _work(self.size(), other.size(), False, budget)
        try:
            return Self(self.value - other.value)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact product after charging its size.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact result.

        Raises:
            PolyError: LIMIT or CAPACITY if the computation cannot proceed.
        """
        _work(self.size(), other.size(), True, budget)
        try:
            return Self(self.value * other.value)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Negate the coefficient after charging its size.

        Args:
            budget: The shared work allowance.

        Returns:
            The additive inverse.

        Raises:
            PolyError: LIMIT if the allowance is insufficient.
        """
        _charge_size(self.size(), budget)
        return Self(-self.value)

@fieldwise_init
struct QQ(CoefficientField):
    """An exact rational coefficient."""
    var value: Rational
    """The canonical arbitrary-precision value."""
    comptime Context = NoContext
    """The unique ring context."""

    @staticmethod
    def zero(ctx: NoContext) -> Self:
        """Construct zero.

        Args:
            ctx: The ring context.

        Returns:
            The additive identity.
        """
        return Self(Rational(0))

    @staticmethod
    def one(ctx: NoContext) -> Self:
        """Construct one.

        Args:
            ctx: The ring context.

        Returns:
            The multiplicative identity.
        """
        return Self(Rational(1))

    def context(self) -> NoContext:
        """Report the unique context.

        Returns:
            The ring's context.
        """
        return NoContext()

    def is_zero(self) -> Bool:
        """Test whether this coefficient is zero.

        Returns:
            Whether the value is zero.
        """
        return self.value == Rational(0)

    def size(self) -> Int:
        """Report a positive magnitude size in bytes.

        Returns:
            The size used by the work model.
        """
        return self.value.numerator().magnitude_bit_length() // 8 + self.value.denominator().magnitude_bit_length() // 8 + 2

    def __eq__(self, other: Self) -> Bool:
        """Compare exact values (the fixed-signature protocol exception).

        Args:
            other: The value to compare.

        Returns:
            Whether the canonical coefficients agree.
        """
        return self.value == other.value

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact sum after charging its size.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact result.

        Raises:
            PolyError: LIMIT or CAPACITY if the computation cannot proceed.
        """
        _work(self.size(), other.size(), True, budget)
        try:
            return Self(self.value + other.value)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact difference after charging its size.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact result.

        Raises:
            PolyError: LIMIT or CAPACITY if the computation cannot proceed.
        """
        _work(self.size(), other.size(), True, budget)
        try:
            return Self(self.value - other.value)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the exact product after charging its size.

        Args:
            other: The other coefficient.
            budget: The shared work allowance.

        Returns:
            The exact result.

        Raises:
            PolyError: LIMIT or CAPACITY if the computation cannot proceed.
        """
        _work(self.size(), other.size(), True, budget)
        try:
            return Self(self.value * other.value)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Negate the coefficient after charging its size.

        Args:
            budget: The shared work allowance.

        Returns:
            The additive inverse.

        Raises:
            PolyError: LIMIT if the allowance is insufficient.
        """
        _charge_size(self.size(), budget)
        return Self(-self.value)

    def inverse(self, mut budget: Budget) raises PolyError -> Self:
        """Invert a nonzero rational.

        Args:
            budget: The shared work allowance.

        Returns:
            The reciprocal.

        Raises:
            PolyError: DOMAIN on zero, LIMIT or CAPACITY on work limits.
        """
        if self.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "inverse of zero")
        _work(self.size(), 0, True, budget)
        try:
            return Self(Rational(1) / self.value)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

def _mod_mul(a: UInt64, b: UInt64, p: UInt64) -> UInt64:
    # Factors below 2^32 multiply within a word, whose remainder is a native
    # division instead of a 128-bit one; the residue is the same.
    if (a | b) >> 32 == 0:
        return (a * b) % p
    return UInt64((UInt128(a) * UInt128(b)) % UInt128(p))

def _mod_pow(a: UInt64, n: UInt64, p: UInt64, mut budget: Budget) raises PolyError -> UInt64:
    var base = a % p
    var exponent = n
    var result = UInt64(1)
    while exponent:
        budget.charge(1)
        if exponent & 1:
            result = _mod_mul(result, base, p)
        exponent >>= 1
        if exponent:
            base = _mod_mul(base, base, p)
    return result

def _prime(p: UInt64, mut budget: Budget) raises PolyError -> Bool:
    if p < 2 or p >= UInt64(1) << 62:
        return False
    if p == 2:
        return True
    if p % 2 == 0:
        return False
    var d = p - 1
    var s = 0
    while d % 2 == 0:
        budget.charge(1)
        d >>= 1
        s += 1
    # Deterministic Miller-Rabin witnesses for every unsigned 64-bit input.
    var witnesses: List[UInt64] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    for a in witnesses:
        budget.charge(1)
        if a % p == 0:
            continue
        var x = _mod_pow(a, d, p, budget)
        if x == 1 or x == p - 1:
            continue
        var accepted = False
        for _ in range(1, s):
            budget.charge(1)
            x = _mod_mul(x, x, p)
            if x == p - 1:
                accepted = True
                break
        if not accepted:
            return False
    return True

struct Modulus(ImplicitlyCopyable, Equatable, Writable):
    """A validated prime strictly below two to the power 62."""
    var value: UInt64
    """The prime; callers must preserve the validated value."""

    def __init__(out self, value: UInt64, mut budget: Budget) raises PolyError:
        """Validate a prime field context using deterministic Miller-Rabin.

        Args:
            value: The proposed prime.
            budget: The shared work allowance.

        Raises:
            PolyError: DOMAIN for a nonprime or out-of-range modulus, LIMIT on exhaustion.
        """
        if not _prime(value, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "modulus must be a prime below 2^62")
        self.value = value

    def __init__(out self, *, _prime_value: UInt64):
        """Recover an already validated internal prime context.

        Args:
            _prime_value: A prime previously checked by the public constructor.
        """
        self.value = _prime_value

    def __eq__(self, other: Self) -> Bool:
        """Compare prime field contexts.

        Args:
            other: The other field context.

        Returns:
            Whether the primes agree.
        """
        return self.value == other.value

struct ModInt(CoefficientField):
    """A canonical residue carrying its prime modulus."""
    var value: UInt64
    """The representative in the range zero through modulus minus one."""
    var modulus: UInt64
    """The validated prime modulus."""
    comptime Context = Modulus
    """The field context type."""

    def __init__(out self, value: UInt64, ctx: Modulus):
        """Reduce a native unsigned integer in a validated field.

        Args:
            value: The integer to reduce.
            ctx: A validated prime context.
        """
        self.value = value % ctx.value
        self.modulus = ctx.value

    @staticmethod
    def zero(ctx: Modulus) -> Self:
        """Construct field zero.

        Args:
            ctx: The prime field.

        Returns:
            Its additive identity.
        """
        return Self(0, ctx)

    @staticmethod
    def one(ctx: Modulus) -> Self:
        """Construct field one.

        Args:
            ctx: The prime field.

        Returns:
            Its multiplicative identity.
        """
        return Self(1, ctx)

    def context(self) -> Modulus:
        """Recover the validated context.

        Returns:
            The field of this residue.
        """
        return Modulus(_prime_value=self.modulus)

    def size(self) -> Int:
        """Report the native residue storage size.

        Returns:
            Eight bytes.
        """
        return 8

    def is_zero(self) -> Bool:
        """Test for the additive identity.

        Returns:
            Whether the residue is zero.
        """
        return self.value == 0

    def __eq__(self, other: Self) -> Bool:
        """Compare residues and their fields.

        Args:
            other: The other residue.

        Returns:
            Whether both value and modulus agree.
        """
        return self.value == other.value and self.modulus == other.modulus

    def _same(self, other: Self) raises PolyError:
        if self.modulus != other.modulus:
            raise PolyError(PolyErrorKind.DOMAIN, "different prime fields")

    def _with_residue(self, value: UInt64) -> Self:
        # Internal arithmetic proves 0 <= value < modulus. Keep the public
        # constructor's reduction for arbitrary native inputs.
        var result = self.copy()
        result.value = value
        return result^

    def add(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the modular sum.

        Args:
            other: The other residue in the same field.
            budget: The shared work allowance.

        Returns:
            The reduced result.

        Raises:
            PolyError: DOMAIN for mixed fields or LIMIT on exhaustion.
        """
        self._same(other)
        budget.charge(1)
        var sum = self.value + other.value  # < 2p < 2^63; cannot overflow.
        return self._with_residue(sum - self.modulus if sum >= self.modulus else sum)

    def sub(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the modular difference.

        Args:
            other: The other residue in the same field.
            budget: The shared work allowance.

        Returns:
            The reduced result.

        Raises:
            PolyError: DOMAIN for mixed fields or LIMIT on exhaustion.
        """
        self._same(other)
        budget.charge(1)
        var difference = self.value - other.value if self.value >= other.value else self.modulus - (other.value - self.value)
        return self._with_residue(difference)

    def mul(self, other: Self, mut budget: Budget) raises PolyError -> Self:
        """Compute the modular product.

        Args:
            other: The other residue in the same field.
            budget: The shared work allowance.

        Returns:
            The reduced result.

        Raises:
            PolyError: DOMAIN for mixed fields or LIMIT on exhaustion.
        """
        self._same(other)
        budget.charge(1)
        return self._with_residue(_mod_mul(self.value, other.value, self.modulus))

    def neg(self, mut budget: Budget) raises PolyError -> Self:
        """Compute the additive inverse in this field.

        Args:
            budget: The shared work allowance.

        Returns:
            The canonical additive inverse.

        Raises:
            PolyError: LIMIT on exhaustion.
        """
        budget.charge(1)
        return self._with_residue(0 if self.is_zero() else self.modulus - self.value)

    def inverse(self, mut budget: Budget) raises PolyError -> Self:
        """Invert by Fermat's theorem using budgeted repeated squaring.

        Args:
            budget: The shared work allowance.

        Returns:
            The multiplicative inverse.

        Raises:
            PolyError: DOMAIN for zero or LIMIT on exhaustion.
        """
        if self.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "inverse of zero")
        return self._with_residue(_mod_pow(self.value, self.modulus - 2, self.modulus, budget))

def _pow[C: CoefficientRing](a: C, n: UInt64, mut budget: Budget) raises PolyError -> C:
    var result = C.one(a.context())
    _charge_size(a.size(), budget)
    var base = a.copy()
    var exponent = n
    while exponent:
        budget.charge(1)
        if exponent & 1:
            result = result.mul(base, budget)
        exponent >>= 1
        if exponent:
            base = base.mul(base, budget)
    return result^
