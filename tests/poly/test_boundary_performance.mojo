"""Independent results and failure boundaries for borrowed arithmetic and ingestion."""
from poly import dense, sparse
from apn_mojo import Integer, Rational, ExactComplex
from polynomials import Budget, PolyError, PolyErrorKind
from std.iter import iter
from std.testing import assert_equal, assert_true


def test_dense_lists() raises:
    for n in range(12):
        for m in range(12):
            var a = List[Rational]()
            var b = List[Int]()
            for i in range(n):
                a.append(Rational(i % 5 - 2, i + 1))
            for i in range(m):
                b.append(2 - i % 5)
            var sum = dense.add(a, b)
            var difference = dense.subtract(a, b)
            var expected_sum = List[Rational]()
            var expected_difference = List[Rational]()
            for i in range(max(n, m)):
                var x = Rational(i % 5 - 2, i + 1) if i < n else Rational(0)
                var y = Rational(2 - i % 5) if i < m else Rational(0)
                expected_sum.append(x + y)
                expected_difference.append(x - y)
            while len(expected_sum) and expected_sum[len(expected_sum) - 1] == Rational(0):
                _ = expected_sum.pop()
            while len(expected_difference) and expected_difference[len(expected_difference) - 1] == Rational(0):
                _ = expected_difference.pop()
            assert_equal(sum.to_list(), expected_sum)
            assert_equal(difference.to_list(), expected_difference)
            assert_equal(len(a), n)
            assert_equal(len(b), m)
    var values: List[Rational] = [Rational(1, 3), Rational(2, 3)]
    var before = dense.add(values, values)
    values[0] = Rational(4, 3)
    assert_equal(before[0], Rational(2, 3))
    assert_equal(dense.add(values, values)[0], Rational(8, 3))
    assert_equal(dense.add([UInt64.MAX], [Int64.MIN])[0], Integer(UInt64.MAX) + Integer(Int64.MIN))
    var unsigned: List[UInt64] = [UInt64.MAX, UInt64(0)]
    var signed: List[Int64] = [Int64.MIN, Int64(0)]
    assert_equal(dense.asarray(unsigned).to_list(), [Integer(UInt64.MAX)])
    assert_equal(dense.asarray(signed).to_list(), [Integer(Int64.MIN)])
    var fixed = [Rational(1, 3), Rational(0)]
    var converted = dense.asarray(fixed)
    fixed[0] = Rational(2)
    assert_equal(converted.to_list(), [Rational(1, 3)])
    assert_equal(dense.asarray(fixed).to_list(), [Rational(2)])
    var imaginary: List[ExactComplex] = [ExactComplex(1, 2), ExactComplex(-1, 3)]
    assert_equal(dense.subtract(imaginary, [Rational(1, 2)]).to_list(), [ExactComplex(Rational(1, 2), 2), ExactComplex(-1, 3)])


def test_sparse_lists() raises:
    var a: List[Tuple[Int, Rational]] = [(0, Rational(1, 2)), (2, Rational(3)), (9, Rational(0))]
    var b: List[Tuple[Int, Int]] = [(0, -1), (1, 2), (2, -3)]
    assert_equal(sparse.add(a, b).to_list(), [(0, Rational(-1, 2)), (1, Rational(2))])
    assert_equal(sparse.subtract(a, b).to_list(), [(0, Rational(3, 2)), (1, Rational(-2)), (2, Rational(6))])
    assert_equal(sparse.add(List[Tuple[Int, Rational]](), b).to_list(), [(0, Rational(-1)), (1, Rational(2)), (2, Rational(-3))])
    # The fallback must normalize repeated powers even when they cancel.
    var unordered: List[Tuple[Int, Int]] = [(2, 4), (0, 1), (2, -4), (1, 2), (0, -1)]
    assert_equal(sparse.add(unordered, [(1, -2)]).to_list(), List[Tuple[Int, Integer]]())
    assert_equal(sparse.subtract(unordered, [(1, 2)]).to_list(), List[Tuple[Int, Integer]]())
    var multi: List[Tuple[Tuple[Int, Int], Rational]] = [((0, 2), Rational(1, 2)), ((1, 0), Rational(3))]
    var other: List[Tuple[Tuple[Int, Int], Int]] = [((0, 1), 2), ((0, 2), 1), ((1, 0), -3)]
    assert_equal(sparse.add(multi, other).to_list(), [((0, 1), Rational(2)), ((0, 2), Rational(3, 2))])
    var kept = sparse.add(a, b)
    b[0] = (0, 1)
    assert_equal(kept[0][1], Rational(-1, 2))
    assert_equal(sparse.add(a, b)[0][1], Rational(3, 2))
    for power in [-1, Int(4294967296)]:
        var failed = False
        try:
            # Invalid powers must not be hidden by a zero coefficient.
            _ = sparse.add([(power, 0)], [(0, 1)])
        except error:
            failed = error.kind == PolyErrorKind.DOMAIN
        assert_true(failed)
    var failed = False
    try:
        _ = sparse.subtract([((0, -1), 0)], [((0, 0), 1)])
    except error:
        failed = error.kind == PolyErrorKind.DOMAIN
    assert_true(failed)


def workflow(scenario: Int, a: List[Rational], b: List[Tuple[Int, Rational]], mut budget: Budget) raises PolyError:
    if scenario == 0:
        _ = dense.asarray(a, budget)
    elif scenario == 1:
        _ = dense.asarray(iter(a), budget)
    elif scenario == 2:
        _ = sparse.asarray(b, budget)
    elif scenario == 3:
        _ = sparse.asarray(iter(b), budget)
    elif scenario == 4:
        _ = dense.add(a, a, budget)
    elif scenario == 5:
        _ = dense.subtract(a, a, budget)
    elif scenario == 6:
        _ = sparse.add(b, b, budget)
    else:
        _ = sparse.subtract(b, b, budget)


def test_failure_and_ownership() raises:
    var a: List[Rational] = [Rational(1, 3), Rational(0), Rational(2)]
    var b: List[Tuple[Int, Rational]] = [(0, Rational(1, 3)), (2, Rational(2))]
    for scenario in range(8):
        var measured = Budget()
        workflow(scenario, a, b, measured)
        var spent = UInt64.MAX - measured.remaining()
        var exact = Budget(spent)
        workflow(scenario, a, b, exact)
        assert_equal(exact.remaining(), UInt64(0))
        # Fail at every step, including growth, conversion, arithmetic and wrap.
        for allowance in range(Int(spent)):
            var budget = Budget(UInt64(allowance))
            var failed = False
            try:
                workflow(scenario, a, b, budget)
            except error:
                failed = error.kind == PolyErrorKind.LIMIT
            assert_true(failed)
            assert_equal(a, [Rational(1, 3), Rational(0), Rational(2)])
            assert_equal(b, [(0, Rational(1, 3)), (2, Rational(2))])
    var huge = Rational(Integer(1) << 2048, 3)
    var source: List[Rational] = [huge.copy(), Rational(0)]
    var result = dense.asarray(iter(source))
    source[0] = Rational(7)
    assert_equal(result.to_list(), [huge.copy()])
    var imaginary: List[ExactComplex] = [ExactComplex(1, 2), ExactComplex(3, 4)]
    assert_equal(dense.asarray(iter(imaginary)).to_list(), imaginary)


def main() raises:
    test_dense_lists()
    test_sparse_lists()
    test_failure_and_ownership()
    print("array boundary performance checks passed")
