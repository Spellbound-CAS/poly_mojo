# expect: context must match the requested dtype
from poly import sparse
from apn_mojo import Ball, ArithmeticContext

def main() raises:
    _ = sparse.astype[Ball]([(0, 1)], context=ArithmeticContext())
