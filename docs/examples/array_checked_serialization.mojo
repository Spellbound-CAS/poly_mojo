from poly import dense, exact, io
from std.testing import assert_true

def main() raises:
    var result = dense.divide_checked([1, 2, 1], [1, 1])
    var saved = io.dumps(result)
    var restored = io.loads[type_of(result)](saved)
    assert_true(exact.check(restored))
    print("Restored division evidence checked")
