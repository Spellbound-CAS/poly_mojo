"""Canonical text roundtrips exact domains, numeric precision and sparse shape."""
from std.testing import assert_equal
from apn_mojo import Integer, Rational, Float, Ball, ArithmeticContext, FloatFormat, BallContext
from poly import dense, sparse, modint
from polynomials import ModInt


def main() raises:
    var text = dense.dumps([1, -2, 3])
    var p = dense.loads[Integer](text)
    assert_equal(p[2], 3)
    assert_equal(dense.dumps(p), text)
    var rational = dense.asarray([Rational(1, 3), Rational(2)])
    assert_equal(dense.dumps(dense.loads[Rational](dense.dumps(rational))), dense.dumps(rational))
    var modular = dense.asarray([modint(3, 7), modint(1, 7)])
    assert_equal(dense.loads[ModInt](dense.dumps(modular), modulus=7).modulus(), UInt64(7))
    var sp = sparse.asarray([((1, 0), 2), ((0, 2), 3)])
    assert_equal(sparse.dumps(sparse.loads[Integer, ndim=2](sparse.dumps(sp))), sparse.dumps(sp))
    var numeric = dense.asarray([Rational(1, 3), Rational(2)], context=ArithmeticContext(format=FloatFormat(precision=80)))
    assert_equal(dense.dumps(dense.loads_numeric[Float](dense.dumps(numeric))), dense.dumps(numeric))
    var balls = sparse.asarray([((1, 0), Rational(1, 3))], context=BallContext(precision=80))
    assert_equal(sparse.dumps(sparse.loads_numeric[Ball, ndim=2](sparse.dumps(balls))), sparse.dumps(balls))
    var constant = sparse.asarray([((), 7)], context=ArithmeticContext(format=FloatFormat(precision=80)))
    var zero_axes = sparse.loads_numeric[Float, ndim=0](sparse.dumps(constant))
    assert_equal(sparse.evaluate(zero_axes, (), context=ArithmeticContext(format=FloatFormat(precision=80))).to_rational_exact(), Rational(7))
    print("array gap serialization checks passed")
