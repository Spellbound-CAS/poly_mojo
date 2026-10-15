"""Behavioral checks for the public polynomial input boundary."""
import poly
from apn_mojo import Integer, Rational, ExactComplex, Batch
from polynomials import ModInt, Budget, PolyError, PolyErrorKind
from poly._common import _allocation
from std.collections import Array
from std.iter import IterableOwned, Iterator, StopIteration, iter
from std.memory import ArcPointer
from std.testing import assert_equal, assert_true

@fieldwise_init
struct Counters(Movable):
    var iterations: Int
    var reads: Int

struct Once[T: Copyable & Deinitable](IterableOwned, Iterator):
    comptime Element = Self.T
    comptime IteratorOwnedType: Iterator = Self
    var values: List[Self.T]
    var index: Int
    var counters: ArcPointer[Counters]

    def __init__(out self, var values: List[Self.T], counters: ArcPointer[Counters]):
        self.values = values^
        self.index = 0
        self.counters = counters

    def __iter__(var self) -> Self:
        self.counters[].iterations += 1
        return self^

    def __next__(mut self) raises StopIteration -> Self.T:
        self.counters[].reads += 1
        if self.index == len(self.values):
            raise StopIteration()
        var result = self.values[self.index].copy()
        self.index += 1
        return result^

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        return 0, None


def test_dense() raises:
    var original = [1, 2, 0]
    var p = poly.dense.asarray(original)
    assert_equal(p.domain(), "ZZ")
    assert_equal(p.to_list(), [Integer(1), Integer(2)])
    assert_equal(original, [1, 2, 0])
    assert_equal(poly.dense.evaluate([1, 2, 3], 2), Integer(17))
    var q = poly.dense.add(p, [-1, -2])
    assert_equal(len(q), 0)
    assert_equal(q.domain(), "ZZ")
    assert_equal(poly.dense.multiply([1, 2], [3, 4]).to_list(), [Integer(3), Integer(10), Integer(8)])
    assert_equal(poly.dense.derivative([1, 2, 3]).to_list(), [Integer(2), Integer(6)])
    assert_equal(len(poly.dense.asarray(List[Int]())), 0)
    var fixed: Array[Int, 3] = [1, 2, 0]
    assert_equal(poly.dense.asarray(fixed).to_list(), p.to_list())
    assert_equal(poly.dense.asarray(Span(fixed)).to_list(), p.to_list())
    var batch = Batch[Integer]([1, 2, 0])
    assert_equal(poly.dense.asarray(batch).to_list(), p.to_list())
    assert_equal(poly.dense.asarray(Batch[Rational]([1, 2])).domain(), "QQ")
    assert_equal(poly.dense.asarray([UInt64.MAX])[0], Integer(UInt64.MAX))
    assert_equal(poly.dense.asarray([Int64.MIN])[0], Integer(Int64.MIN))
    var list_copy = p.to_list()
    list_copy[0] = Integer(99)
    assert_equal(p[0], Integer(1))
    var total = Integer(0)
    for value in p:
        total += value
    assert_equal(total, Integer(3))


def test_promotion() raises:
    var rational_values = [1, Rational(1, 2)]
    var a = poly.dense.asarray(rational_values)
    var b = poly.dense.asarray([Rational(1, 2), 1])
    assert_equal(a.domain(), "QQ")
    assert_equal(b.domain(), "QQ")
    assert_equal(a[1], Rational(1, 2))
    var n = Int(3)
    var r = Rational(1, 2)
    var mixed = poly.coefficients((n, r, ExactComplex(0, 1)))
    var typed_rational: Rational = poly.dense.evaluate([0, 1], r)
    var typed_complex: ExactComplex = poly.sparse.evaluate([(1, 1)], ExactComplex(0, 1))
    assert_equal(typed_rational, r)
    assert_equal(typed_complex, ExactComplex(0, 1))
    assert_equal(mixed.domain(), "QQi")
    assert_equal(poly.dense.evaluate(mixed, 2), ExactComplex(4, 4))
    var complex_zero = poly.dense.asarray([ExactComplex(0, 0)])
    assert_equal(complex_zero.domain(), "QQi")
    assert_equal(len(complex_zero), 0)
    var rational_zero = poly.dense.asarray([Rational(0)])
    assert_equal(rational_zero.domain(), "QQ")
    assert_equal(poly.dense.evaluate([0, 1], Rational(1, 2)), Rational(1, 2))
    assert_equal(poly.dense.add([1], [Rational(1, 2)])[0], Rational(3, 2))
    assert_equal(poly.dense.add([Rational(1, 2)], [1])[0], Rational(3, 2))
    var c = poly.dense.add([Rational(1, 2)], [ExactComplex(0, 1)])
    assert_equal(c[0], ExactComplex(Rational(1, 2), 1))


def test_modular() raises:
    var m = poly.dense.asarray([-1, 2], modulus=7)
    assert_equal(m.domain(), "GF(7)")
    assert_equal(m[0].value, UInt64(6))
    assert_equal(poly.dense.evaluate(m, -1).value, UInt64(4))
    var mixed = poly.coefficients((Int(-2), poly.modint(3, 7)))
    assert_equal(mixed[0].value, UInt64(5))
    var inferred = poly.dense.asarray([poly.modint(-1, 7), poly.modint(2, 7)])
    assert_equal(inferred.to_list(), m.to_list())
    var zero = poly.dense.add(m, [1, -2])
    assert_equal(len(zero), 0)
    assert_equal(zero.modulus(), UInt64(7))
    assert_equal(poly.dense.add(zero, [8])[0].value, UInt64(1))
    assert_equal(poly.dense.asarray(List[Int](), modulus=7).modulus(), UInt64(7))
    assert_equal(poly.dense.asarray(List[ModInt](), modulus=7).modulus(), UInt64(7))
    assert_equal(poly.modint(Integer(2) ** 130, 7).value, UInt64(2))


def test_sparse() raises:
    # Construct the endpoint at runtime so compiler literal narrowing cannot
    # turn the maximum valid power into -1 before the API receives it.
    var maximum_power = Int("4294967295")
    assert_equal(UInt64(maximum_power), UInt64(4294967295))
    assert_equal(poly.sparse.asarray([(maximum_power, 1)])[0][0], maximum_power)
    var source = [(2, 1), (0, -1), (2, 2), (1, 0)]
    var p = poly.sparse.asarray(source)
    assert_equal(p.to_list(), [(0, Integer(-1)), (2, Integer(3))])
    assert_equal(source, [(2, 1), (0, -1), (2, 2), (1, 0)])
    assert_equal(poly.sparse.evaluate(p, 2), Integer(11))
    assert_equal(poly.sparse.derivative(p).to_list(), [(1, Integer(6))])
    assert_equal(poly.sparse.add([(0, 1)], [(0, Rational(1, 2))])[0][1], Rational(3, 2))
    var b = poly.sparse.asarray([((1, 0), 2), ((0, 2), 3)])
    assert_equal(poly.sparse.evaluate(b, (2, 3)), Integer(31))
    assert_equal(poly.sparse.evaluate(b, (Rational(1, 2), 1)), Rational(4))
    assert_equal(poly.sparse.derivative(b, axis=1).to_list(), [((0, 1), Integer(6))])
    assert_equal(poly.sparse.multiply(b, [((0, 1), 1)]).to_list(), [((0, 3), Integer(3)), ((1, 1), Integer(2))])
    var empty = poly.sparse.asarray([((1, 0), 1), ((1, 0), -1)])
    assert_equal(len(empty), 0)
    assert_equal(poly.sparse.evaluate(empty, (2, 3)), Integer(0))
    assert_equal(poly.sparse.add(empty, [((0, 1), 2)])[0][0], (0, 1))
    var tri = poly.sparse.asarray([((1, 2, 3), 2)])
    assert_equal(poly.sparse.evaluate(tri, (2, 3, 4)), Integer(2304))
    var modular = poly.sparse.asarray([(0, -1), (3, 2)], modulus=7)
    assert_equal(poly.sparse.evaluate(modular, 2).value, UInt64(1))
    var zero = poly.sparse.add(modular, [(0, 1), (3, -2)])
    assert_equal(len(zero), 0)
    assert_equal(zero.modulus(), UInt64(7))
    assert_equal(poly.sparse.add(zero, [(0, 8)])[0][1].value, UInt64(1))
    var fixed: Array[Tuple[Int, Int], 2] = [(0, 1), (2, 3)]
    assert_equal(poly.sparse.evaluate(fixed, 2), Integer(13))
    assert_equal(poly.sparse.evaluate(Span(fixed), 2), Integer(13))


def test_single_pass() raises:
    var calls = ArcPointer(Counters(0, 0))
    var source = Once[Int]([1, 2, 0], calls)
    var p = poly.dense.asarray(source^)
    assert_equal(p.to_list(), [Integer(1), Integer(2)])
    assert_equal(calls[].iterations, 1)
    assert_equal(calls[].reads, 4)
    var terms_calls = ArcPointer(Counters(0, 0))
    var terms = Once[Tuple[Int, Int]]([(2, 3), (0, 1)], terms_calls)
    var s = poly.sparse.asarray(terms^)
    assert_equal(poly.sparse.evaluate(s, 2), Integer(13))
    assert_equal(terms_calls[].iterations, 1)
    assert_equal(terms_calls[].reads, 3)
    var raw = [3, 4]
    var cursor = iter(raw)
    assert_equal(poly.dense.asarray(cursor^).to_list(), [Integer(3), Integer(4)])



def expect_error[operation: def() raises thin -> None](message: String) raises:
    var failed = False
    try:
        operation()
    except error:
        failed = True
        assert_true(message in String(error), String("unexpected error: ", error))
    assert_true(failed, String("invalid input was accepted: ", message))


def reject_mixed_moduli() raises:
    _ = poly.dense.asarray([poly.modint(1, 7), poly.modint(1, 11)])


def reject_mixed_moduli_zero() raises:
    _ = poly.sparse.asarray([(0, poly.modint(0, 7)), (0, poly.modint(0, 11))])


def reject_explicit_mismatch() raises:
    _ = poly.dense.asarray([poly.modint(1, 7)], modulus=11)


def reject_composite_modulus() raises:
    _ = poly.dense.asarray([1], modulus=9)


def reject_missing_modulus() raises:
    _ = poly.dense.asarray(List[ModInt]())


def reject_negative_power() raises:
    _ = poly.sparse.asarray([(-1, 0)])


def reject_oversize_power() raises:
    _ = poly.sparse.asarray([(Int("4294967296"), 1)])


def reject_negative_multivariate_power() raises:
    _ = poly.sparse.asarray([((0, -1), 1)])


def reject_bad_axis() raises:
    _ = poly.sparse.derivative([((1, 0), 1)], axis=2)


def reject_batch_rank() raises:
    _ = poly.dense.asarray(Batch[Integer]([1, 2, 3, 4]).reshape([2, 2]))


def reject_owned_batch_rank() raises:
    var batch = Batch[Rational]([1, 2, 3, 4]).reshape([2, 2])
    _ = poly.dense.asarray(batch^)


def reject_zero_modulus() raises:
    _ = poly.dense.asarray([1], modulus=0)


def test_rejections() raises:
    expect_error[reject_mixed_moduli]("different moduli")
    expect_error[reject_mixed_moduli_zero]("different moduli")
    expect_error[reject_explicit_mismatch]("different moduli")
    expect_error[reject_composite_modulus]("prime")
    expect_error[reject_missing_modulus]("empty modular input")
    expect_error[reject_negative_power]("powers must")
    expect_error[reject_oversize_power]("powers must")
    expect_error[reject_negative_multivariate_power]("powers must")
    expect_error[reject_bad_axis]("axis is outside")
    expect_error[reject_batch_rank]("coefficient axis")
    expect_error[reject_owned_batch_rank]("coefficient axis")
    expect_error[reject_zero_modulus]("prime")



def test_mixed_sparse() raises:
    var p = poly.sparse.asarray(((0, 1), (1, Rational(1, 2)), (0, 2)))
    assert_equal(p.domain(), "QQ")
    var value: Rational = poly.sparse.evaluate(p, 2)
    assert_equal(value, Rational(4))
    assert_equal(p.to_list(), [(0, Rational(3)), (1, Rational(1, 2))])
    var z = poly.sparse.asarray((((1, 0), 2), ((0, 1), ExactComplex(1, 2)), ((1, 0), Rational(1, 2))))
    var complex_value: ExactComplex = poly.sparse.evaluate(z, (2, 3))
    assert_equal(complex_value, ExactComplex(8, 6))
    assert_equal(z.domain(), "QQi")
    var zero = poly.sparse.asarray((((2, 0), Rational(0)), ((2, 0), 0)))
    assert_equal(len(zero), 0)
    assert_equal(zero.domain(), "QQ")
    assert_equal(poly.sparse.evaluate(zero, (3, 4)), Rational(0))
    var reversed = poly.sparse.asarray(((1, Rational(1, 2)), (0, 3)))
    assert_equal(reversed.to_list(), p.to_list())
    var wide = poly.sparse.asarray(((0, UInt64.MAX), (1, Int64.MIN)))
    assert_equal(wide[0][1], Integer(UInt64.MAX))
    assert_equal(wide[1][1], Integer(Int64.MIN))
    var m = poly.sparse.asarray(((0, -1), (2, poly.modint(3, 7))))
    assert_equal(m.modulus(), UInt64(7))
    assert_equal(poly.sparse.evaluate(m, 2).value, UInt64(4))
    var explicit = poly.sparse.asarray(((0, Int8(-1)), (1, UInt64(2))), modulus=7)
    assert_equal(poly.sparse.evaluate(explicit, 2).value, UInt64(3))
    var modzero = poly.sparse.asarray(((2, 1), (2, poly.modint(6, 7))))
    assert_equal(len(modzero), 0)
    assert_equal(modzero.modulus(), UInt64(7))
    var imaginary_zero = poly.sparse.asarray(((0, ExactComplex(0, 0)), (0, 0)))
    assert_equal(imaginary_zero.domain(), "QQi")
    assert_equal(len(imaginary_zero), 0)
    var copied = z.to_list()
    copied[0] = ((0, 0), ExactComplex(99))
    assert_equal(poly.sparse.evaluate(z, (2, 3)), ExactComplex(8, 6))


def run_budgeted(scenario: Int, half: Rational, mut budget: Budget) raises PolyError:
    # Every branch is a complete public workflow, including the shared output copy.
    if scenario == 0:
        _ = poly.dense.asarray([1, 2, 0], budget).to_list(budget)
    elif scenario == 1:
        _ = poly.dense.asarray([-1, 2], budget, modulus=7).to_list(budget)
    elif scenario == 2:
        _ = poly.dense.add([1, 2], [half], budget).to_list(budget)
    elif scenario == 3:
        _ = poly.dense.multiply([1, 1], [1, -1], budget).to_list(budget)
    elif scenario == 4:
        _ = poly.dense.evaluate([1, 2, 3], half, budget)
    elif scenario == 5:
        _ = poly.dense.derivative([1, 2, 3], budget).to_list(budget)
    elif scenario == 6:
        var values = [1, 2, 3]
        _ = poly.dense.asarray(iter(values), budget).to_list(budget)
    elif scenario == 7:
        var values = [-1, 2]
        _ = poly.dense.asarray(iter(values), budget, modulus=7).to_list(budget)
    elif scenario == 8:
        _ = poly.sparse.asarray([(2, 3), (0, 1), (2, -3)], budget).to_list(budget)
    elif scenario == 9:
        _ = poly.sparse.asarray([(2, -1)], budget, modulus=7).to_list(budget)
    elif scenario == 10:
        _ = poly.sparse.asarray(((0, 1), (1, half)), budget).to_list(budget)
    elif scenario == 11:
        _ = poly.sparse.asarray(((0, 1), (1, Int8(-1))), budget, modulus=7).to_list(budget)
    elif scenario == 12:
        _ = poly.sparse.add([(1, 1)], [(0, half)], budget).to_list(budget)
    elif scenario == 13:
        _ = poly.sparse.multiply([(0, 1), (1, 1)], [(0, 1), (1, -1)], budget).to_list(budget)
    elif scenario == 14:
        _ = poly.sparse.evaluate([(0, 1), (1, 2)], half, budget)
    elif scenario == 15:
        _ = poly.sparse.evaluate([((1, 2), 3)], (1, half), budget)
    elif scenario == 16:
        _ = poly.sparse.derivative([((1, 2), 3)], budget, axis=1).to_list(budget)
    elif scenario == 17:
        var values = [(1, 2), (0, 3)]
        _ = poly.sparse.asarray(iter(values), budget).to_list(budget)
    elif scenario == 18:
        var values = [(1, -2), (0, 3)]
        _ = poly.sparse.asarray(iter(values), budget, modulus=7).to_list(budget)
    elif scenario == 19:
        _ = poly.coefficients((1, half, ExactComplex(0, 1)), budget).to_list(budget)
    else:
        _ = poly.modint(Integer(-99), 7, budget)


def test_budgets() raises:
    var half = Rational(1, 2)
    for scenario in range(21):
        var measured = Budget()
        run_budgeted(scenario, half, measured)
        var spent = UInt64.MAX - measured.remaining()
        assert_true(spent > 0)
        var exact = Budget(spent)
        run_budgeted(scenario, half, exact)
        assert_equal(exact.remaining(), UInt64(0))
        var allowances: List[UInt64] = [0, spent - 1]
        for allowance in allowances:
            var limited = Budget(allowance)
            var failed = False
            try:
                run_budgeted(scenario, half, limited)
            except error:
                failed = True
                assert_equal(error.kind, PolyErrorKind.LIMIT)
                assert_true(limited.remaining() <= allowance)
            assert_true(failed, String("missing budget charge in scenario ", scenario))

    var p = poly.dense.asarray([1, 2, 3])
    var sparse = poly.sparse.asarray([((1, 2), 3), ((0, 0), 1)])
    var zero_budget = Budget(0)
    var shared = poly.dense.asarray(p, zero_budget)
    var sparse_shared = poly.sparse.asarray(sparse, zero_budget)
    assert_equal(shared.to_list(), p.to_list())
    assert_equal(sparse_shared.to_list(), sparse.to_list())
    assert_true(Pointer(to=shared._core[]) == Pointer(to=p._core[]))
    assert_true(Pointer(to=sparse_shared._core[]) == Pointer(to=sparse._core[]))
    var failed = False
    try:
        _ = poly.dense.multiply(p, p, max_steps=1)
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.LIMIT)
    assert_true(failed)
    assert_equal(poly.dense.evaluate(p, 2), Integer(17))
    failed = False
    try:
        _ = sparse.to_list(max_steps=0)
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.LIMIT)
    assert_true(failed)
    assert_equal(poly.sparse.evaluate(sparse, (2, 3)), Integer(55))

    # A zero allowance must stop before calling next(), even for an owned source.
    var calls = ArcPointer(Counters(0, 0))
    var source = Once([1, 2, 3], calls)
    failed = False
    try:
        _ = poly.dense.asarray(source^, max_steps=0)
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.LIMIT)
    assert_true(failed)
    assert_equal(calls[].reads, 0)
    var term_calls = ArcPointer(Counters(0, 0))
    var term_source = Once([(1, 2), (0, 3)], term_calls)
    failed = False
    try:
        _ = poly.sparse.asarray(term_source^, max_steps=1)
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.LIMIT)
    assert_true(failed)
    assert_equal(term_calls[].reads, 1)

    # Invalid input is a typed domain failure, including a malformed zero term.
    failed = False
    try:
        _ = poly.sparse.asarray(((0, half), (-1, 0)))
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.DOMAIN)
    assert_true(failed)
    var seven = poly.modint(0, 7)
    var five = poly.modint(0, 5)
    failed = False
    try:
        _ = poly.sparse.asarray(((0, seven.copy()), (1, five.copy())))
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.DOMAIN)
    assert_true(failed)

    # Byte-size overflow is rejected before reserving any element storage.
    var untouched = List[SIMD[DType.uint64, 4]]()
    var capacity_budget = Budget()
    failed = False
    try:
        _allocation[SIMD[DType.uint64, 4]](Int.MAX, capacity_budget)
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.CAPACITY)
    assert_true(failed)
    assert_equal(capacity_budget.remaining(), UInt64.MAX)
    assert_equal(len(untouched), 0)



def test_modular_context_reuse() raises:
    # A large prime makes accidental per-coefficient primality checks visible.
    var value = poly.modint(1, UInt64(2305843009213693951))
    var once: List[ModInt] = [value.copy()]
    var repeated = List[ModInt](length=20, fill=value)
    var first_budget = Budget()
    var many_budget = Budget()
    _ = poly.dense.asarray(once, first_budget)
    var prepared = poly.dense.asarray(repeated, many_budget)
    var extra = first_budget.remaining() - many_budget.remaining()
    assert_true(extra < 512, "field validation repeated for every coefficient")
    var no_work = Budget(0)
    assert_equal(poly.dense.asarray(prepared, no_work).modulus(), value.modulus)
    assert_equal(poly.dense.asarray(prepared, no_work, modulus=value.modulus).modulus(), value.modulus)

    var invalid = value.copy()
    invalid.modulus = 0
    var input: List[ModInt] = [value.copy(), invalid.copy()]
    var failed = False
    try:
        _ = poly.dense.asarray(input)
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.DOMAIN)
    assert_true(failed)


def test_prepared_arithmetic() raises:
    var p = poly.dense.asarray([1, 2, 3])
    var q = poly.dense.asarray([3, -2])
    var reused = Budget()
    var raw = Budget()
    assert_equal(poly.dense.add(p, q, reused).to_list(), [Integer(4), Integer(0), Integer(3)])
    _ = poly.dense.add([1, 2, 3], [3, -2], raw)
    assert_true(reused.remaining() > raw.remaining())
    var s = poly.sparse.asarray([(0, 1), (2, 3)])
    var t = poly.sparse.asarray([(0, -1), (1, 2)])
    reused = Budget()
    raw = Budget()
    assert_equal(poly.sparse.add(s, t, reused).to_list(), [(1, Integer(2)), (2, Integer(3))])
    _ = poly.sparse.add([(0, 1), (2, 3)], [(0, -1), (1, 2)], raw)
    assert_true(reused.remaining() > raw.remaining())



def test_compaction() raises:
    var coefficients = List[Int](length=1000, fill=0)
    coefficients[0] = 1
    coefficients[1] = 2
    coefficients[2] = 3
    coefficients[3] = 4
    coefficients[4] = 5
    var p = poly.dense.asarray(coefficients)
    assert_equal(p.to_list(), [Integer(1), Integer(2), Integer(3), Integer(4), Integer(5)])
    assert_true(p._core[].coefficients.capacity() <= 2 * len(p))
    assert_equal(len(coefficients), 1000)
    assert_equal(coefficients[4], 5)

    var full = poly.dense.asarray(List[Int](length=1000, fill=1))
    var negative = poly.dense.asarray(List[Int](length=1000, fill=-1))
    var zero = poly.dense.add(full, negative)
    assert_equal(len(zero), 0)
    assert_equal(zero._core[].coefficients.capacity(), 0)
    assert_equal(len(full), 1000)

    var left = List[Tuple[Int, Int]]()
    var right = List[Tuple[Int, Int]]()
    for i in range(1000):
        left.append((i, i + 1))
        right.append((i, -i - 1))
    var sparse_zero = poly.sparse.add(left, right)
    assert_equal(sparse_zero._core[].terms.capacity(), 0)
    # Leave five coefficients to exercise compaction's order-preserving moves.
    right[0] = (0, 0)
    right[1] = (1, 0)
    right[2] = (2, 0)
    right[3] = (3, 0)
    right[4] = (4, 0)
    var sparse_small = poly.sparse.add(left, right)
    assert_equal(sparse_small.to_list(), [(0, Integer(1)), (1, Integer(2)), (2, Integer(3)), (3, Integer(4)), (4, Integer(5))])
    assert_true(sparse_small._core[].terms.capacity() <= 2 * len(sparse_small))

    # A compaction failure must leave borrowed source coefficients unchanged.
    var measured = Budget()
    _ = poly.dense.asarray(coefficients, measured)
    var almost = Budget(UInt64.MAX - measured.remaining() - 2)
    var failed = False
    try:
        _ = poly.dense.asarray(coefficients, almost)
    except error:
        failed = True
        assert_equal(error.kind, PolyErrorKind.LIMIT)
    assert_true(failed)
    assert_equal(coefficients[0], 1)
    assert_equal(coefficients[4], 5)


def main() raises:
    test_dense()
    test_promotion()
    test_modular()
    test_sparse()
    test_single_pass()
    test_rejections()
    test_mixed_sparse()
    test_budgets()
    test_modular_context_reuse()
    test_prepared_arithmetic()
    test_compaction()
    print("array boundary checks passed")
