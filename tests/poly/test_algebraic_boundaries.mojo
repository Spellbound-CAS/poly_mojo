"""Algebraic witness normalization, immutable ownership, and single-pass input budgets."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ
from poly import dense, sparse, algebraic
from std.iter import IterableOwned, Iterator, StopIteration
from std.memory import ArcPointer
from std.testing import assert_equal, assert_true


@fieldwise_init
struct Counts(Movable):
    var reads: Int


@fieldwise_init
struct Once(IterableOwned, Iterator):
    comptime Element = algebraic.Real
    comptime IteratorOwnedType: Iterator = Self
    var values: List[algebraic.Real]
    var index: Int
    var counts: ArcPointer[Counts]

    def __iter__(var self) -> Self:
        return self^

    def __next__(mut self) raises StopIteration -> algebraic.Real:
        self.counts[].reads += 1
        if self.index == len(self.values):
            raise StopIteration()
        var value = self.values[self.index].copy()
        self.index += 1
        return value^


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var witness = a.witness().copy()
    var budget = Budget()
    var multiplier = dense.to_core([Rational(-3), Rational(1)])
    witness.polynomial = witness.polynomial.mul(multiplier, budget)
    var restored = algebraic.asreal(witness)
    assert_true(restored == a)
    assert_true(algebraic.check(restored))
    witness.polynomial.coefficients[0] = QQ(Rational(100))
    assert_true(restored == a)
    assert_true(algebraic.check(restored))
    var malformed = False
    try:
        _ = algebraic.asreal(witness)
    except e:
        malformed = e.kind == PolyErrorKind.DOMAIN
    assert_true(malformed)
    var counts = ArcPointer(Counts(0))
    var once = Once([a.copy(), algebraic.asreal(1)], 0, counts)
    var p = dense.asarray(once^)
    assert_equal(p.domain(), "AA")
    assert_equal(counts[].reads, 3)
    counts = ArcPointer(Counts(0))
    once = Once([a.copy()], 0, counts)
    var limited = False
    try:
        _ = dense.asarray(once^, max_steps=0)
    except e:
        limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_equal(counts[].reads, 0)
    for steps in [0, 1, 20, 100]:
        var work = Budget(UInt64(steps))
        limited = False
        try:
            _ = dense.multiply(p, p, work)
        except e:
            limited = e.kind == PolyErrorKind.LIMIT
        assert_true(limited)
        assert_true(work.remaining() <= UInt64(steps))
        assert_true(p[0] == a)
    var empty = sparse.asarray(List[Tuple[Tuple[Int, Int], algebraic.Complex]]())
    assert_equal(empty.domain(), "QQbar")
    assert_equal(algebraic.coefficients(empty, axis=1).domain(), "QQbar")
    assert_equal(len(algebraic.coefficients(empty, axis=1)), 0)
    print("array algebraic boundary checks passed")
