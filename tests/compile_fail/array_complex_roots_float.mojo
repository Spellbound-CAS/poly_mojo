# expect: exact complex roots require characteristic-zero exact coefficients
from poly import dense

def main() raises:
    _ = dense.complex_roots([1.0, 0.0, 1.0])
