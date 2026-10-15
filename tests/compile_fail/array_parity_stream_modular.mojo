# expect: enclosure evaluation requires characteristic-zero exact coefficients
from poly import dense
from apn_mojo import BallContext

def main() raises:
    _ = dense.enclose_stream(dense.asarray([1, 2], modulus=7), 3, context=BallContext(precision=80))
