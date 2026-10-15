# expect: real roots require exact integer or rational coefficients
from poly import dense

def main() raises:
    _ = dense.real_roots([1.0, 2.0])
