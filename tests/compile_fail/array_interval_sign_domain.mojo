# expect: real roots require exact integer or rational coefficients
from poly import dense

def main() raises:
    _ = dense.sign_on_interval(dense.asarray([1, 2], modulus=5), 0, 1)
