"""Exact binary endpoint extraction against public rational bounds and identities."""
from apn_mojo import Integer, Rational, Ball, ComplexBall
from polynomials import *
from polynomials.numeric_domains import _numeric_cost
from polynomials.numeric_bounds import _evaluation_bounds
from support.check import Checker


def exercise(value: Ball, mut check: Checker, mut work: Budget) raises:
    var box = numeric_bounds[BallDomain](value, work)
    check.check(_evaluation_bounds[BallDomain](value, work) == box, "producer endpoints equal independent rational endpoints")
    check.check(box.real_low.value == value.lower_rational() and box.real_high.value == value.upper_rational(), "binary endpoints equal the public rational construction")
    check.check(box.imag_low.is_zero() and box.imag_high.is_zero(), "real endpoint extraction has zero imaginary bounds")
    var z = ComplexBall(value, Ball(-2, precision=value.midpoint().precision()))
    var rectangle = numeric_bounds[ComplexBallDomain](z, work)
    check.check(_evaluation_bounds[ComplexBallDomain](z, work) == rectangle, "complex producer endpoints equal independent rational endpoints")
    check.check(rectangle.real_low == box.real_low and rectangle.real_high == box.real_high and rectangle.imag_low.value == Rational(-2) and rectangle.imag_high.value == Rational(-2), "complex endpoint extraction preserves both coordinates")
    var tiny = Budget(0)
    try:
        _ = numeric_bounds[BallDomain](value, tiny)
        check.check(False, "endpoint extraction ignored its budget")
    except e:
        check.check(e.kind == PolyErrorKind.LIMIT, "endpoint exhaustion preserves LIMIT")


def main() raises:
    var check = Checker("exact numerical endpoints")
    var work = Budget()
    for bits in [2, 17, 64, 256, 4096]:
        for numerator in [-7, -1, 0, 1, 7]:
            for denominator in [1, 2, 3, 17]:
                for radius in [0, 1, 2, 7]:
                    exercise(Ball(Rational(numerator, denominator), Rational(radius, 19), precision=bits), check, work)
    for exponent in [-10000, -511, -1, 0, 1, 511, 10000]:
        var power = Integer(1) << abs(exponent)
        var scale = Rational(power) if exponent >= 0 else Rational(Integer(1), power)
        for negative in [False, True]:
            exercise(Ball(-scale if negative else scale, Rational(1, 3), precision=128), check, work)
            exercise(Ball(Rational(-1 if negative else 1, 3), scale, precision=128), check, work)
    for midpoint in [-2, -1, 0, 1, 2]:
        var value = Ball(midpoint, 1, precision=32)
        var box = numeric_bounds[BallDomain](value, work)
        check.check(box.real_low.value == Rational(midpoint - 1) and box.real_high.value == Rational(midpoint + 1), "hand-computed endpoints include exact cancellation and radius dominance")
    var nonexact = Ball(Rational(1, 3), precision=256)
    var limited = Budget(UInt64(_numeric_cost[BallDomain](nonexact)))
    try:
        _ = _evaluation_bounds[BallDomain](nonexact, limited)
        check.check(False, "endpoint alignment ignored its budget")
    except e:
        check.check(e.kind == PolyErrorKind.LIMIT, "alignment exhaustion is not reclassified as CAPACITY")
    check.finish()
