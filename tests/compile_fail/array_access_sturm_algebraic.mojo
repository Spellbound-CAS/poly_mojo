# expect: Sturm preparation requires integer or rational coefficients
from poly import dense, algebraic

def main() raises:
    _ = dense.sturm_sequence([algebraic.asreal(1)])
