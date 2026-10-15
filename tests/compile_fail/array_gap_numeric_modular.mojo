# expect: numeric coefficients require real/complex numbers
from poly import dense, modint
from apn_mojo import ArithmeticContext

def main() raises:
    _ = dense.asarray([modint(1, 7)], context=ArithmeticContext())
