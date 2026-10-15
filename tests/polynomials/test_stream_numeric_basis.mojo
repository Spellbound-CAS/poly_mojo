"""Closed-form monomial coefficients in the Chebyshev basis."""
from apn_mojo import Integer, Rational, FloatFormat, ArithmeticContext, ComplexContext, BallContext
from polynomials import *
from support.check import Checker


def domains[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    for n in [2, 16, 64]:
        var zero = numeric_value[D](NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), work), context, work)
        var input = NumericBox.point(QQ(Rational(1)), QQ(Rational(0)), work)
        comptime if D == BallDomain:
            input = NumericBox.real_interval(QQ(Rational(-1)), QQ(Rational(2)), work)
        elif D == ComplexBallDomain:
            input = NumericBox(QQ(Rational(-1)), QQ(Rational(2)), QQ(Rational(-3)), QQ(Rational(4)))
        var leading = numeric_value[D](input, context, work)
        var values = List[D.Scalar](length=n + 1, fill=zero)
        values[n] = leading^
        var p = NumericDensePoly[D].from_coefficients(0, values^, work)
        var basis = NumericBasisPolynomial[D].chebyshev(p, context, work)
        var choose = QQ(Rational(1))
        var scale = QQ(Rational(Integer(1), Integer(1) << (n - 1)))
        # x^n = 2^(1-n) sum_{j<n/2} binomial(n,j) T_(n-2j)
        #       + 2^(-n) binomial(n,n/2) T_0, for positive even n.
        for j in range(n // 2 + 1):
            var factor = choose.mul(scale, work)
            if j == n // 2:
                factor = factor.mul(QQ(Rational(1, 2)), work)
            var expected = NumericBox(input.real_low.mul(factor, work), input.real_high.mul(factor, work), input.imag_low.mul(factor, work), input.imag_high.mul(factor, work))
            c.check(numeric_bounds[D](basis.coefficients[n - 2 * j], work).contains(expected, work), "closed-form Chebyshev coefficient enclosure")
            choose = choose.mul(QQ(Rational(n - j)), work).mul(QQ(Rational(j + 1)).inverse(work), work)
        for j in range(1, n, 2):
            c.check(numeric_bounds[D](basis.coefficients[j], work).contains(NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), work), work), "Chebyshev parity zero")


def main() raises:
    var c = Checker("stream numeric basis")
    var work = Budget()
    domains[FloatDomain](ArithmeticContext(format=FloatFormat(128)), c, work)
    domains[ComplexDomain](ComplexContext(ArithmeticContext(format=FloatFormat(128))), c, work)
    domains[BallDomain](BallContext(128), c, work)
    domains[ComplexBallDomain](BallContext(128), c, work)
    c.finish()
