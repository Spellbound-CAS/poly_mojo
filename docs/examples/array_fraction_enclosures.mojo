from poly import rational, io
from apn_mojo import BallContext
from std.testing import assert_true

def main() raises:
    var f = rational.asarray([1], [1, 1])
    var context = BallContext(precision=96)
    var result = rational.enclose(f, 1, context=context)
    assert_true(rational.check(result))
    var restored = io.loads[type_of(result)](io.dumps(result))
    assert_true(rational.check(restored))
    print("Certified fraction value and restored evidence checked")
