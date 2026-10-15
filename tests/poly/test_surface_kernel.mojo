"""Independent kernel bytes, mapped axes, malformed inputs, and resource limits."""
from poly import dense, sparse
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def main() raises:
    # One term, one pair (ID 9, exponent 1), coefficient 2/1.
    var expected: List[UInt8] = [1,0,0,0, 1,0,0,0, 9,0,0,0,0,0,0,0, 1,0,0,0, 7,0,0,0, 2,0,0,0, 0,2,1]
    var encoded = dense.to_kernel_bytes([0, 2], variable=9)
    assert_equal(encoded, expected)
    assert_equal(dense.from_kernel_bytes(Span(expected), variable=9).to_list(), List[Rational]([Rational(0), Rational(2)]))
    var p = sparse.asarray([((2, 0, 0), Rational(1, 3)), ((0, 1, 0), Rational(-7, 2)), ((0, 0, 0), Rational(5))])
    var ids: List[Int64] = [19, -7, 400]
    var payload = sparse.to_kernel_bytes(p, variables=ids)
    assert_equal(sparse.from_kernel_bytes[ndim=3](Span(payload), variables=ids).to_list(), p.to_list())
    assert_equal(dense.from_kernel_bytes(Span(dense.to_kernel_bytes([2])))[0], Rational(2))
    for mode in range(9):
        var rejected = False
        try:
            if mode == 0: _ = dense.to_kernel_bytes([0])
            elif mode == 1: _ = sparse.to_kernel_bytes([(0, 1)])
            elif mode == 2: _ = dense.from_kernel_bytes(Span(expected))
            elif mode == 3: _ = sparse.from_kernel_bytes[ndim=3](Span(payload), variables=List[Int64]([Int64(1), Int64(1), Int64(2)]))
            elif mode == 4: _ = dense.from_kernel_bytes(Span(payload), variable=19)
            else:
                var bad = expected.copy()
                if mode == 5: bad.append(0)
                elif mode == 6: _ = bad.pop()
                elif mode == 7: bad[16] = 0
                else: bad[29] = 0
                _ = dense.from_kernel_bytes(Span(bad), variable=9)
        except e:
            rejected = e.kind == (PolyErrorKind.DOMAIN if mode >= 2 and mode <= 4 else PolyErrorKind.ENCODING)
        assert_true(rejected)
    var meter = Budget()
    _ = sparse.from_kernel_bytes[ndim=3](Span(payload), meter, variables=ids)
    var used = UInt64.MAX - meter.remaining()
    var exact = Budget(used)
    _ = sparse.from_kernel_bytes[ndim=3](Span(payload), exact, variables=ids)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(used - 1)
    var limited = False
    try: _ = sparse.from_kernel_bytes[ndim=3](Span(payload), short, variables=ids)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    var huge = sparse.to_kernel_bytes([(1000000000, 2)])
    limited = False
    try: _ = dense.from_kernel_bytes(Span(huge), max_steps=1000000)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    print("surface kernel adapter checks passed")
