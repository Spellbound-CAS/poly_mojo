# expect: real roots require exact integer or rational coefficients
from poly import dense

def main() raises:
    _ = dense.real_roots(dense.asarray([1, 2], modulus=5))
