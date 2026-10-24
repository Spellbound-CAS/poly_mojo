"""Sparse array products, powers and quotients agree with sparse-born arrays, including at the packed kernels' key-range limit."""
import poly
from apn_mojo import Integer
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_equal, assert_true, assert_false


def base() raises PolyError -> poly.sparse.Array[Tuple[Int, Int, Int], 0]:
    """1 + x + y + z."""
    return poly.sparse.asarray([((0, 0, 0), 1), ((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), 1)])


def binomial(n: Int, k: Int) -> Int:
    var result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


def test_packed_results() raises:
    var f = base()
    var f7 = poly.sparse.power(f, 7)
    assert_equal(len(f7), binomial(10, 3))
    assert_equal(poly.sparse.evaluate(f7, (2, 3, 5)), Integer(11) ** 7)
    # Chained packed operands give the same terms as a sparse-born round trip.
    var f3 = poly.sparse.power(f, 3)
    var f4 = poly.sparse.power(f, 4)
    var product = poly.sparse.multiply(f3, f4)
    assert_equal(product.to_list(), f7.to_list())
    var rebuilt = poly.sparse.asarray(f7.to_list())
    assert_equal(poly.sparse.dumps(rebuilt), poly.sparse.dumps(product))
    var quotient = poly.sparse.divide_exact(product, f3)
    assert_true(Bool(quotient))
    assert_equal(quotient.value().to_list(), f4.to_list())
    var shifted = poly.sparse.add(product, [((0, 0, 0), 1)])
    assert_false(Bool(poly.sparse.divide_exact(shifted, f3)))
    # Indexing and iteration materialize the same terms in order.
    var listed = f7.to_list()
    var position = 0
    for entry in f7:
        assert_equal(entry[0], listed[position][0])
        assert_equal(entry[1], listed[position][1])
        position += 1
    assert_equal(position, len(listed))
    assert_equal(f7[len(f7) - 1][0], listed[len(listed) - 1][0])


def test_sparse_born_operands() raises:
    # Large sparse-born operands take the packed kernels; the result matches the power.
    var f = base()
    var f8 = poly.sparse.asarray(poly.sparse.power(f, 8).to_list())
    var square = poly.sparse.multiply(f8, f8)
    assert_equal(square.to_list(), poly.sparse.power(f, 16).to_list())
    var back = poly.sparse.divide_exact(poly.sparse.asarray(square.to_list()), f8)
    assert_equal(back.value().to_list(), f8.to_list())


def test_modular() raises:
    var f = poly.sparse.asarray([((0, 0, 0), 3), ((2, 0, 0), 5), ((0, 1, 3), -1)], modulus=101)
    var f9 = poly.sparse.power(f, 9)
    assert_equal(f9.domain(), "GF(101)")
    var f4 = poly.sparse.power(f, 4)
    var f5 = poly.sparse.power(f, 5)
    assert_equal(poly.sparse.multiply(f4, f5).to_list(), f9.to_list())
    assert_equal(poly.sparse.divide_exact(f9, f4).value().to_list(), f5.to_list())
    assert_equal(poly.sparse.evaluate(f9, (2, 3, 5)), poly.sparse.evaluate(poly.sparse.asarray(f9.to_list(), modulus=101), (2, 3, 5)))


def test_wide_exponents() raises:
    # A product whose packed key range would reach 2**62 stays sparse.
    var wide = poly.sparse.asarray([((1 << 30, 1 << 30), 1), ((0, 0), 1)])
    var square = poly.sparse.power(wide, 2)
    assert_equal(len(square), 3)
    assert_equal(square[2][0], (1 << 31, 1 << 31))
    assert_equal(poly.sparse.multiply(wide, wide).to_list(), square.to_list())


def chain(mut budget: Budget) raises PolyError:
    var f = base()
    var f6 = poly.sparse.power(f, 6, budget)
    var g = poly.sparse.multiply(f6, f6, budget)
    _ = poly.sparse.divide_exact(g, f6, budget)
    _ = g.to_list(budget)


def test_budget() raises:
    var measured = Budget()
    chain(measured)
    var used = UInt64.MAX - measured.remaining()
    var again = Budget()
    chain(again)
    assert_equal(UInt64.MAX - again.remaining(), used)
    var exact = Budget(used)
    chain(exact)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(used - 1)
    var failed = False
    try:
        chain(short)
    except error:
        assert_equal(error.kind, PolyErrorKind.LIMIT)
        failed = True
    assert_true(failed)


def test_charges_ignore_materialization() raises:
    # A packed result charges the same later work whether or not it was materialized.
    var f6 = poly.sparse.power(base(), 6)
    var first = Budget()
    _ = poly.sparse.multiply(f6, f6, first)
    _ = f6.to_list()
    var second = Budget()
    _ = poly.sparse.multiply(f6, f6, second)
    assert_equal(first.remaining(), second.remaining())


def main() raises:
    test_packed_results()
    test_sparse_born_operands()
    test_modular()
    test_wide_exponents()
    test_budget()
    test_charges_ignore_materialization()
    print("sparse packed array checks passed")
