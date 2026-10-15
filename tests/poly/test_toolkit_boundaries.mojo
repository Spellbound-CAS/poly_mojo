"""Precision, shape, single-pass preparation, work accounting, and bad evidence."""
from poly import series, basis, rational, dense, sparse
from apn_mojo import Rational
from polynomials import Budget, PolyError, PolyErrorKind, QQ
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from std.testing import assert_equal, assert_true, assert_false


@fieldwise_init
struct Cursor(Iterator):
    comptime Element = Int
    var position: ArcPointer[Int]
    var reads: ArcPointer[Int]
    var values: ArcPointer[List[Int]]

    def __next__(mut self) raises StopIteration -> Int:
        self.reads[] += 1
        if self.position[] == len(self.values[]):
            raise StopIteration()
        var result = self.values[][self.position[]]
        self.position[] += 1
        return result


@fieldwise_init
struct Once(Copyable, Movable, Iterable):
    var position: ArcPointer[Int]
    var reads: ArcPointer[Int]
    var values: ArcPointer[List[Int]]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = Cursor

    def __iter__(ref self) -> Cursor:
        return Cursor(self.position.copy(), self.reads.copy(), self.values.copy())


def source() -> Once:
    var values: List[Int] = [1, 0, 0]
    return Once(ArcPointer(Int(0)), ArcPointer(Int(0)), ArcPointer(values^))


def workflow(s: series.Array[1], b: basis.Array[1], r: rational.Function[Int], mut budget: Budget) raises PolyError:
    _ = series.coefficients(series.reciprocal(s, budget), budget)
    _ = basis.coefficients(basis.subdivide(b, 1, budget)[0], budget)
    _ = rational.check(rational.partial_fractions(r, budget), budget)
    _ = rational.check(rational.laurent(r, 0, budget, terms=3), budget)


def main() raises:
    var first = source()
    var s = series.asarray(first, order=5)
    assert_equal(first.reads[], 4)
    assert_equal(s[4], Rational(0))
    var second = source()
    var b = basis.asarray(second, basis="bernstein")
    assert_equal(second.reads[], 4)
    assert_equal(len(b), 3)
    var third = source()
    var r = rational.asarray(third, [0, 1, -1])
    assert_equal(third.reads[], 4)
    var sites = source()
    var node_basis = basis.asarray([1, 2, 3, 4], basis="newton", nodes=sites)
    assert_equal(sites.reads[], 4)
    assert_equal(len(basis.nodes(node_basis)), 3)
    # Factorization populates the core's process-wide prime cache on first use.
    # Compare exact allowances only after that separate initialization work.
    assert_true(rational.check(rational.partial_fractions(r)))
    var work = Budget()
    workflow(s, b, r, work)
    var cost = UInt64.MAX - work.remaining()
    var exact = Budget(cost)
    workflow(s, b, r, exact)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(cost - 1)
    var failed = False
    try:
        workflow(s, b, r, short)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    assert_equal(s[0], Rational(1))
    assert_equal(b[0], Rational(1))
    assert_equal(rational.evaluate(r, 2), Rational(-1, 2))
    var mod5 = series.asarray(dense.asarray([1], modulus=5), order=3)
    var mod7 = series.asarray(dense.asarray([1], modulus=7), order=3)
    var char2 = series.asarray(dense.asarray([0, 1], modulus=2), order=3)
    for scenario in range(26):
        failed = False
        try:
            if scenario == 0:
                _ = series.asarray([1])
            elif scenario == 1:
                _ = series.asarray([1], order=-2)
            elif scenario == 2:
                _ = series.asarray(s, order=6)
            elif scenario == 3:
                _ = series.add(s, series.truncate(s, 2))
            elif scenario == 4:
                _ = series.compose(s, [1])
            elif scenario == 5:
                _ = series.reciprocal([0, 1], order=4)
            elif scenario == 6:
                _ = series.log([2], order=3)
            elif scenario == 7:
                _ = series.exp([1], order=3)
            elif scenario == 8:
                _ = series.revert([0, 0, 1], order=4)
            elif scenario == 9:
                _ = series.add(mod5, mod7)
            elif scenario == 10:
                _ = series.integral(char2)
            elif scenario == 11:
                _ = basis.asarray([1, 2], basis="newton", nodes=[0, 1])
            elif scenario == 12:
                _ = basis.asarray(List[Int](), basis="bernstein")
            elif scenario == 13:
                _ = basis.asarray([1], basis="unsupported")
            elif scenario == 14:
                _ = basis.from_polynomial([0, 0, 1], basis="bernstein", degree=1)
            elif scenario == 15:
                _ = basis.from_polynomial([1], basis="chebyshev", nodes=[0])
            elif scenario == 16:
                _ = basis.subdivide(basis.asarray([1], basis="chebyshev"), 0)
            elif scenario == 17:
                _ = basis.multiply(basis.asarray([0, 1], basis="newton", nodes=[0]), basis.asarray([0, 1], basis="newton", nodes=[0]))
            elif scenario == 18:
                _ = basis.add(b, basis.asarray([1], basis="chebyshev"))
            elif scenario == 19:
                _ = rational.asarray([1], [0])
            elif scenario == 20:
                _ = rational.divide(r, rational.asarray([0]))
            elif scenario == 21:
                _ = rational.evaluate(r, 0)
            elif scenario == 22:
                _ = rational.laurent(r, 0, terms=-1)
            elif scenario == 23:
                _ = rational.derivative(r, axis=1)
            elif scenario == 24:
                _ = rational.substitute(r, [sparse.asarray([(0, 0)])])
            else:
                _ = rational.poles(rational.from_terms([((0, 1), 1)], [((1, 0), 1)]))
        except e:
            failed = e.kind == PolyErrorKind.DOMAIN
        assert_true(failed, String("scenario ", scenario))
    failed = False
    try:
        _ = series.coefficients(series.asarray([1], order=Int.MAX), max_steps=10)
    except e:
        failed = e.kind == PolyErrorKind.CAPACITY or e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var raw = series.to_core(s)
    raw.order = 0
    failed = False
    try:
        _ = series.from_core(raw)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var invalid_basis = basis.to_core(b)
    invalid_basis.kind = 9
    failed = False
    try:
        _ = basis.from_core(invalid_basis)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var altered = rational.asarray([1], [0, 1])
    altered._data[].value.numerator.terms[0].coefficient = QQ(Rational(2))
    assert_false(rational.check(altered))
    print("array toolkit boundary checks passed")
