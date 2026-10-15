# expect: basis queries cannot widen the prepared coefficient field
from poly import sparse
from apn_mojo import ExactComplex

def main() raises:
    var b = sparse.groebner([[(1, 1)]])
    _ = sparse.reduce([(1, ExactComplex(0, 1))], b)
