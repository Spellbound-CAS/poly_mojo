# expect: evaluation point has the wrong number of coordinates
from poly import sparse
from apn_mojo import BallContext

def main() raises:
    _ = sparse.enclose([((1, 0), 1)], (1,), context=BallContext(precision=80))
