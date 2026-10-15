# expect: evaluation point has the wrong number of coordinates
from poly import sparse
from apn_mojo import ArithmeticContext

def main() raises:
    var context = ArithmeticContext()
    var p = sparse.asarray([((1, 0), 1)], context=context)
    _ = sparse.evaluate(p, (1,), context=context)
