"""Unknown constants and rounded underflow cannot justify formal composition."""
from apn_mojo import Rational, Float, FloatFormat, ArithmeticContext, ComplexContext, BallContext
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, TruncatedSeries, NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, NumericDensePoly, NumericSeries
from support.check import Checker


def numeric[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    var unknown = NumericSeries[D].from_polynomial(NumericDensePoly[D](0), 0, context, work)
    try:
        _ = unknown.compose(unknown, context, work)
        c.check(False, "order-zero numerical composition accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "order-zero numerical constant is unknown")


def main() raises:
    var c = Checker("series boundaries")
    var work = Budget()
    var unknown = TruncatedSeries[QQ].from_polynomial(DensePoly[QQ](0, NoContext()), 0, work)
    try:
        _ = unknown.compose(unknown, work)
        c.check(False, "order-zero exact composition accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "order-zero exact constant is unknown")
    numeric[FloatDomain](ArithmeticContext(), c, work)
    numeric[ComplexDomain](ComplexContext(ArithmeticContext()), c, work)
    numeric[BallDomain](BallContext(128), c, work)
    numeric[ComplexBallDomain](BallContext(128), c, work)
    var poly = NumericDensePoly[FloatDomain].from_coefficients(0, [Float(Rational(1, 1024)), Float(1)], work)
    var narrow = ArithmeticContext(format=FloatFormat(2, emin=-2, emax=2))
    var series = NumericSeries[FloatDomain].from_polynomial(poly, 2, narrow, work)
    try:
        _ = series.exp(narrow, work)
        c.check(False, "underflow erased a nonzero formal constant")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "formal constants checked before rounding")
    c.finish()
