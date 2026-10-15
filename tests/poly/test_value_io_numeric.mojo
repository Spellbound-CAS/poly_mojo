"""Numerical records preserve formats, zeros, basis extent, and enclosures."""
from std.testing import assert_equal, assert_true
from apn_mojo import Float, Complex, Ball, ComplexBall, Rational, ArithmeticContext, ComplexContext, BallContext, FloatFormat, contains
from poly import io, series, basis, dense, sparse


def same[T: Copyable & Deinitable](value: T) raises:
    var text = io.dumps(value)
    assert_equal(io.dumps(io.loads[T](text)), text)


def main() raises:
    var ctx = ArithmeticContext(format=FloatFormat(precision=96))
    var bc = BallContext(precision=96)
    var cc = ComplexContext(ctx)
    var s = series.asarray([Rational(1, 3)], order=6, context=ctx)
    var copy = io.loads[type_of(s)](io.dumps(s))
    assert_equal(copy.order(), 6)
    assert_equal(copy[0].to_json(), s[0].to_json())
    assert_equal(copy[5].to_json(), s[5].to_json())
    assert_equal(copy[5].precision(), 96)
    same(s)
    same(series.asarray(List[Int](), order=0, context=ctx))
    same(series.asarray(List[Int](), order=10, context=bc))
    same(series.asarray([Complex(1, 1)], order=3, context=cc))
    same(series.asarray([Complex(1, 1)], order=3, context=bc))
    for name in [String("newton"), String("chebyshev"), String("bernstein")]:
        var b = basis.asarray([Rational(1, 3), Rational(0)], basis=name, nodes=List[Int]([2]) if name == "newton" else List[Int](), context=bc)
        var saved = io.loads[type_of(b)](io.dumps(b))
        assert_equal(len(saved), 2)
        assert_equal(saved[1].to_json(), b[1].to_json())
        assert_true(contains(basis.evaluate(saved, 0, context=bc), Rational(1, 3)))
        same(b)
    same(basis.asarray([Complex(1, 1), Complex(0, 0)], basis="newton", nodes=[Complex(2, 1)], context=cc))
    same(basis.asarray([Complex(1, 1), Complex(0, 0)], basis="bernstein", context=bc))
    same(basis.asarray(List[Int](), basis="chebyshev", context=ctx))
    same(dense.asarray([1, 2], context=ctx))
    same(sparse.asarray([((0, 2, 0), Rational(1, 3))], context=bc))
    same(sparse.asarray([((), 0)], context=ctx))
    same(Float(0, context=ctx))
    same(Complex(0, 0, context=cc))
    same(Ball(Rational(1, 3), precision=96))
    same(ComplexBall(Complex(1, 1), precision=96))
    print("value serialization numerical checks passed")
