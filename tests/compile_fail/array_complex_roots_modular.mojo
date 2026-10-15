# expect: exact complex roots require characteristic-zero exact coefficients
from poly import dense, modint

def main() raises:
    _ = dense.complex_roots([modint(1, 7), modint(1, 7)])
