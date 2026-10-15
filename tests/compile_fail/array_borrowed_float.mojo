# expect: exact coefficients must be integers
from poly import dense

def main() raises:
    var values: List[Float64] = [1.5, 2.5]
    _ = dense.asarray(values)
