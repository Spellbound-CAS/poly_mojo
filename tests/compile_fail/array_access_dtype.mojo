# expect: dtype must be Float, Complex, Ball, or ComplexBall
from poly import dense
from apn_mojo import ArithmeticContext

def main() raises:
    _ = dense.astype[Int]([1], context=ArithmeticContext())
