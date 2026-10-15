"""Checked exact results retain their operands and proof across archives."""
from poly import dense, sparse, exact, io, algebraic, coefficients, modint
from polynomials import PolyError, PolyErrorKind
from std.testing import assert_true, assert_equal

def roundtrip[k: Int, P: Copyable & Deinitable, sparse: Bool, E: Copyable & Deinitable](value: exact.Result[k, P, sparse, E]) raises:
    comptime T = exact.Result[k, P, sparse, E]
    var text = io.dumps(value)
    var restored = io.loads[T](text)
    assert_true(exact.check(restored))
    assert_equal(io.dumps(restored), text)
    var failed = False
    try:
        _ = io.loads[T](String(text[byte=0:text.byte_length()-2]))
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)

def main() raises:
    roundtrip(dense.divide_checked([1, 2, 1], [1, 1]))
    roundtrip(dense.gcd_checked([1, 2, 1], [1, 1]))
    roundtrip(dense.resultant_certified([1, 0, 1], [1, 1]))
    roundtrip(dense.resultant_certified([modint(1, 5), modint(0, 5), modint(1, 5)], [modint(1, 5), modint(1, 5)]))
    print("completion exact archive checks passed")
