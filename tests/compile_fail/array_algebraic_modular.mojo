# expect: modular coefficients cannot mix implicitly with rational or complex coefficients
from poly import dense, algebraic, coefficients, modint

def main() raises:
    _ = coefficients((algebraic.asreal(1), modint(1, 7)))
