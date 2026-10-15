# expect: Sturm preparation requires integer or rational coefficients
from poly import sparse, modint

def main() raises:
    _ = sparse.sturm_sequence([(1, modint(1, 7))])
