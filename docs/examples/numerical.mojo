"""Rounded arithmetic, rigorous evaluation, and validated real/complex roots."""
from apn_mojo import Rational, Float, FloatFormat, ArithmeticContext, Ball, BallContext
from polynomials import Budget, QQ, NoContext, DensePoly, FloatDomain, BallDomain, NumericDensePoly, NumericOptions, EnclosureStream, EnclosureStreamChecker, enclose_dense, check_enclose_dense, real_roots_exact, check_real_roots_exact, complex_roots_exact, check_complex_roots_exact


def main() raises:
    var budget = Budget()
    var context = ArithmeticContext(format=FloatFormat(128))
    var numerical = NumericDensePoly[FloatDomain].from_coefficients(0, [Float(1), Float(2), Float(1)], budget)
    print("Rounded value at 3:", numerical.evaluate(Float(3), context, budget))
    var exact = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], budget)
    var argument = Ball(Rational(3, 2), Rational(1, 16))
    var enclosure = enclose_dense[BallDomain](exact, argument, BallContext(128), budget)
    print("Enclosure checked:", check_enclose_dense[BallDomain](exact, argument, enclosure, budget))
    var stream = EnclosureStream[BallDomain].from_dense(exact, argument, BallContext(128), budget)
    var checker = EnclosureStreamChecker[BallDomain].from_dense(exact, argument, budget)
    while True:
        var bounds = stream.next(budget)
        if not bounds:
            break
        if not checker.accept(bounds.take(), budget):
            raise Error("streamed bound rejected")
    if not checker.finish(stream.finish(budget), budget):
        raise Error("streamed evaluation is incomplete or invalid")
    var options = NumericOptions(absolute_bits=4)
    var real = real_roots_exact(exact, options, budget, complete_only=True)
    print("Real root cover checked:", check_real_roots_exact(exact, real, options, budget))
    var complex = complex_roots_exact(exact, options, budget, complete_only=True)
    print("Complex roots:", len(complex.regions))
    print("Complex root cover checked:", check_complex_roots_exact(exact, complex, options, budget))
