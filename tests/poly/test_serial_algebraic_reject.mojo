"""Reject malformed metadata, noncanonical coefficients and forged root evidence."""
from poly import dense, sparse, algebraic
from polynomials import PolyErrorKind, QQ
from apn_mojo import Rational
from std.testing import assert_true, assert_equal


def bad(text: String) raises:
    var rejected = False
    try:
        _ = dense.loads[algebraic.Real](text)
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(rejected)


def main() raises:
    var prefix = "POLYNOMIALS-ALGEBRAIC 1\nDENSE\nAA\n0\n1\n"
    var body = "3\n2\n2\n2\n-2\n1\n"
    var fixture = String(prefix, body)
    for length in range(fixture.byte_length()):
        bad(String(fixture[byte=0:length]))
    bad(String(fixture, "\n"))
    bad(String("POLYNOMIALS-ALGEBRAIC 2\nDENSE\nAA\n0\n1\n", body))
    bad(String("POLYNOMIALS-ALGEBRAIC 1\nDENSE\nQQbar\n0\n1\n", body))
    bad(String("POLYNOMIALS-ALGEBRAIC 1\nDENSE\nAA\n-0\n1\n", body))
    bad(String("POLYNOMIALS-ALGEBRAIC 1\nDENSE\nAA\n0\n01\n", body))
    bad(String(prefix, "3\n2/1\n2\n2\n-2\n1\n"))
    bad(String(prefix, "3\n2\n3\n2\n-2\n1\n"))
    bad(String(prefix, "3\n2\n2\n2\n-3\n1\n"))
    bad(String(prefix, "3\n0\n0\n2\n0\n1\n"))
    bad(String(prefix, "3\n2\n2\n3\n-2\n1\n0\n"))
    # Tags 0 and 1 carried Sturm sequences and are no longer accepted.
    bad(String(prefix, "1\n2\n2\n2\n-2\n1\n0\n0\n"))
    # x^2-2 on (1,2): one sign change and Descartes bound 1.
    var irrational = "2\n1\n2\n3\n-2\n0\n1\n"
    var root = dense.loads[algebraic.Real](String(prefix, irrational))[0]
    assert_true(root == dense.root_between([-2, 0, 1], 1, 2))
    bad(String(prefix, "2\n-2\n2\n3\n-2\n0\n1\n"))
    bad(String(prefix, "2\n2\n1\n3\n-2\n0\n1\n"))
    var invalid = root.copy()
    invalid._data[].plain.low = QQ(Rational(-2))
    var rejected = False
    try:
        _ = dense.dumps([invalid.copy()])
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    assert_true(rejected)
    rejected = False
    try:
        _ = sparse.loads[algebraic.Real, ndim=2](fixture)
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(rejected)
    var sparse_prefix = "POLYNOMIALS-ALGEBRAIC 1\nSPARSE\n1\n2\n0\n1\n"
    # Domain is mandatory; unknown or mismatched layouts never get repaired.
    rejected = False
    try:
        _ = sparse.loads[algebraic.Real, ndim=2](String(sparse_prefix, "1\n1\n0\n1\n", body))
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(rejected)
    # Valid sparse metadata with repeated monomials is still noncanonical.
    var sparse_head = "POLYNOMIALS-ALGEBRAIC 1\nSPARSE\nAA\n1\n2\n0\n1\n"
    for payload in [String(sparse_head, "2\n1\n0\n1\n", body, "1\n0\n1\n", body), String(sparse_head, "1\n1\n0\n0\n", body), String(sparse_head, "1\n1\n3\n1\n", body)]:
        rejected = False
        try:
            _ = sparse.loads[algebraic.Real, ndim=2](payload)
        except e:
            rejected = e.kind == PolyErrorKind.ENCODING
        assert_true(rejected)
    rejected = False
    try:
        _ = dense.loads[algebraic.Real]("POLYNOMIALS-ALGEBRAIC 1\nDENSE\nAA\n0\n18446744073709551616\n")
    except e:
        rejected = e.kind == PolyErrorKind.CAPACITY
    assert_true(rejected)
    rejected = False
    try:
        _ = dense.loads[algebraic.Complex](String("POLYNOMIALS-ALGEBRAIC 1\nDENSE\nQQbar\n0\n1\n", body, "3\n2\n3\n2\n-2\n1\n"))
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(rejected)
    print("array algebraic serialization rejection checks passed")
