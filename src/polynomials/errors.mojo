"""Typed failures, deterministic work budgets, and checked size arithmetic."""
from std.sys import size_of

@fieldwise_init
struct PolyErrorKind(ImplicitlyCopyable, Equatable, Writable):
    """Classify a polynomial operation's failure."""
    var code: UInt8
    """The stable kind code."""
    comptime LIMIT = Self(1)
    """The work budget is exhausted."""
    comptime CAPACITY = Self(2)
    """A size or work counter cannot be represented."""
    comptime DOMAIN = Self(3)
    """The operands do not belong to the operation's domain."""
    comptime UNSUPPORTED = Self(4)
    """The operation is not supported for these operands."""
    comptime ENCODING = Self(5)
    """The bytes do not encode a canonical polynomial."""

    def __eq__(self, other: Self) -> Bool:
        """Compare error kinds.

        Args:
            other: The other kind.

        Returns:
            Whether both kinds agree.
        """
        return self.code == other.code

@fieldwise_init
struct PolyError(Copyable, Movable, Writable):
    """A typed failure, with no implicit conversion from plain Error."""
    var kind: PolyErrorKind
    """The machine-readable failure kind."""
    var message: String
    """The human-readable diagnostic."""

struct Budget(ImplicitlyCopyable, Writable):
    """A deterministic abstract-work allowance shared by nested calls."""
    var _remaining: UInt64

    def __init__(out self, steps: UInt64 = UInt64.MAX):
        """Create a budget.

        Args:
            steps: The initial work allowance.
        """
        self._remaining = steps

    def charge(mut self, n: UInt64) raises PolyError:
        """Charge before performing work; a rejected charge changes nothing.

        Args:
            n: The work to charge.

        Raises:
            PolyError: LIMIT if the allowance is insufficient.
        """
        if n > self._remaining:
            raise PolyError(PolyErrorKind.LIMIT, "work budget exhausted")
        self._remaining -= n

    def remaining(self) -> UInt64:
        """Return the unspent allowance.

        Returns:
            The remaining abstract work units.
        """
        return self._remaining

def _sum(a: Int, b: Int) raises PolyError -> Int:
    if a < 0 or b < 0 or a > Int.MAX - b:
        raise PolyError(PolyErrorKind.CAPACITY, "size sum overflow")
    return a + b

def _product(a: Int, b: Int) raises PolyError -> Int:
    if a < 0 or b < 0 or (b != 0 and a > Int.MAX // b):
        raise PolyError(PolyErrorKind.CAPACITY, "size product overflow")
    return a * b

def _charge_size(n: Int, mut budget: Budget) raises PolyError:
    if n < 0:
        raise PolyError(PolyErrorKind.CAPACITY, "negative size")
    budget.charge(UInt64(n))


def _allocation[T: AnyType](count: Int, mut budget: Budget) raises PolyError:
    # Check byte multiplication before allowing a count-sized allocation.
    _ = _product(count, size_of[T]())
    _charge_size(count, budget)
