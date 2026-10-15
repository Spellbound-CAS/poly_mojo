# expect: sparse tuples need numeric pairs with one power shape
from poly import sparse
from apn_mojo import Ball, BallContext

def main() raises:
    _ = sparse.astype[Ball](((1, 2), ((1, 0), 3)), context=BallContext(precision=80))
