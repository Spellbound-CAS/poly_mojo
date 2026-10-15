"""Array numeric."""
from apn_mojo import ArithmeticContext, FloatFormat, BallContext, Rational
from poly import dense, sparse


def main() raises:
    var context = ArithmeticContext(format=FloatFormat(precision=80))
    var p = dense.asarray([1, 2, 3], context=context)
    print("Value:", dense.evaluate(p, 2, context=context).to_rational_exact())
    var balls = dense.asarray([Rational(1, 3), Rational(2)], context=BallContext(precision=80))
    var enclosure = dense.evaluate(balls, 1, context=BallContext(precision=80))
    print("Finite enclosure:", enclosure.is_finite())
    var q = sparse.asarray([((1, 0), 2), ((0, 1), 3)], context=context)
    print("Sparse value:", sparse.evaluate(q, (2, 1), context=context).to_rational_exact())
