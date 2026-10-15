"""The complete call shares one budget, including nested reconstruction and proofs."""
from std.testing import assert_equal, assert_true
from poly import io, series, basis, rational, batch
from polynomials import Budget, PolyErrorKind


def boundaries[T: Copyable & Deinitable](value: T) raises:
    var work = Budget()
    var text = io.dumps(value, work)
    var spent = UInt64.MAX - work.remaining()
    var exact = Budget(spent)
    assert_equal(io.dumps(value, exact), text)
    assert_equal(exact.remaining(), UInt64(0))
    var failed = False
    try:
        _ = io.dumps(value, max_steps=spent - 1)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    work = Budget()
    _ = io.loads[T](text, work)
    spent = UInt64.MAX - work.remaining()
    exact = Budget(spent)
    _ = io.loads[T](text, exact)
    assert_equal(exact.remaining(), UInt64(0))
    failed = False
    try:
        _ = io.loads[T](text, max_steps=spent - 1)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    assert_equal(io.dumps(value), text)


def main() raises:
    boundaries(series.asarray([1, 2], order=8))
    boundaries(basis.asarray([1, 2, 0], basis="bernstein"))
    boundaries(rational.asarray([1], [1, 1]))
    boundaries(batch.asarray([series.asarray([1], order=5), series.asarray([2], order=10)]))
    var failed = False
    try:
        _ = io.loads[series.Array[1]]("", max_field_bytes=0)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    failed = False
    try:
        _ = io.loads[series.Array[1]]("POLYNOMIALS-VALUE 1\nSERIES(QQ)\n9999999999999999999999999999999\n")
    except e:
        failed = e.kind == PolyErrorKind.CAPACITY
    assert_true(failed)
    # A forged count never reserves its claimed capacity before fields arrive.
    failed = False
    try:
        _ = io.loads[basis.Array[1]]("POLYNOMIALS-VALUE 1\nBASIS(QQ)\nchebyshev\n9223372036854775807\n", max_steps=100000)
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    failed = False
    try:
        _ = io.loads[batch.Array[Int]]("POLYNOMIALS-VALUE 1\nBATCH(Int)\n2\n9223372036854775807\n2\n")
    except e:
        failed = e.kind == PolyErrorKind.CAPACITY
    assert_true(failed)
    print("value serialization budget checks passed")
