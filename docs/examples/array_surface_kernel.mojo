"""Exchange a rational polynomial using the kernel literal wire format."""
from poly import dense, sparse
from apn_mojo import Rational
from std.testing import assert_equal


def main() raises:
    var payload = dense.to_kernel_bytes([0, 2], variable=9)
    var restored = dense.from_kernel_bytes(Span(payload), variable=9)
    assert_equal(restored[1], Rational(2))
    print("Kernel payload bytes:", len(payload))
    var p = sparse.asarray([((2, 0, 0), 3), ((0, 1, 0), -1)])
    var ids: List[Int64] = [19, -7, 400]
    var data = sparse.to_kernel_bytes(p, variables=ids)
    var q = sparse.from_kernel_bytes[ndim=3](Span(data), variables=ids)
    assert_equal(sparse.evaluate(q, (2, 3, 0)), Rational(9))
    print("Restored", q.domain(), "coefficients with three axes")
