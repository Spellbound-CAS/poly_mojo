"""Choose a numerical destination and re-round prepared coefficients."""
from apn_mojo import Float, ComplexBall, Rational, ArithmeticContext, FloatFormat, BallContext, contains
from poly import dense, sparse
from std.testing import assert_equal, assert_true


def main() raises:
    var high = ArithmeticContext(format=FloatFormat(precision=100))
    var low = ArithmeticContext(format=FloatFormat(precision=30))
    var p = dense.asarray([Rational(1, 3), Rational(2)], context=high)
    var rounded = dense.astype[Float](p, context=low)
    assert_equal(rounded[0].precision(), 30)
    assert_equal(p[0].precision(), 100)
    var enclosed = sparse.astype[ComplexBall]([(0, Rational(1, 3)), (2, Rational(1))], context=BallContext(precision=80))
    assert_true(contains(enclosed[0][1].real(), Rational(1, 3)))
    assert_true(contains(enclosed[0][1].imag(), 0))
    print("Rounded precision:", rounded[0].precision())
    print("Real coefficients enclosed in complex balls")
