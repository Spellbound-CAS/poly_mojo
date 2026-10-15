"""Handle work exhaustion separately from division by the zero polynomial."""
import poly
from polynomials import PolyErrorKind
from std.testing import assert_true


def main() raises:
    var saw_limit = False
    try:
        _ = poly.dense.power([1, 1], 100, max_steps=1000)
    except error:
        saw_limit = error.kind == PolyErrorKind.LIMIT
    assert_true(saw_limit)

    var saw_domain = False
    try:
        _ = poly.dense.divmod([1, 1], [0])
    except error:
        saw_domain = error.kind == PolyErrorKind.DOMAIN
    assert_true(saw_domain)
    print("Work exhaustion handled:", saw_limit)
    print("Zero divisor handled:", saw_domain)
