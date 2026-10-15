"""Handle work exhaustion separately from an invalid polynomial expression."""
from std.testing import assert_true
from polynomials import Budget, PolyErrorKind, QQ, NoContext, PolynomialRing


def main() raises:
    var exhausted = Budget(0)
    var saw_limit = False
    try:
        exhausted.charge(1)
    except error:
        saw_limit = error.kind == PolyErrorKind.LIMIT
    assert_true(saw_limit)

    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var saw_domain = False
    try:
        _ = ring.parse("y + 1", work)
    except error:
        saw_domain = error.kind == PolyErrorKind.DOMAIN
    assert_true(saw_domain)
    print("Work exhaustion handled:", saw_limit)
    print("Unknown variable handled:", saw_domain)
