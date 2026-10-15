"""Malformed records fail closed instead of silently normalizing saved data."""
from std.testing import assert_true
from apn_mojo import Integer, Float, ArithmeticContext, FloatFormat
from poly import io, series, basis, dense, sparse, rational, batch, algebraic, modint
from polynomials import PolyErrorKind


def bad[T: Copyable & Deinitable](text: String) raises:
    var rejected = False
    try:
        _ = io.loads[T](text)
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(rejected)


def main() raises:
    var good = String("POLYNOMIALS-VALUE 1\nSERIES(QQ)\n5\n2\n1\n2\n")
    for length in range(good.byte_length()):
        bad[series.Array[1]](String(good[byte=0:length]))
    bad[series.Array[1]](String(good, "\n"))
    bad[series.Array[1]](good.replace("VALUE 1", "VALUE 2"))
    bad[series.Array[1]](good.replace("5\n2", "1\n2"))
    bad[series.Array[1]](good.replace("5\n", "05\n"))
    bad[series.Array[1]](good.replace("5\n", "-1\n"))
    bad[series.Array[1]](good.replace("1\n2\n", "1/1\n2\n"))
    bad[series.Array[1]](good.replace("2\n1\n2\n", "2\n1\n0\n"))
    bad[dense.Array[1]](good)
    bad[series.Array[2]](good)
    bad[Integer]("POLYNOMIALS-VALUE 1\nZZ\n1/2\n")
    bad[basis.Array[1]]("POLYNOMIALS-VALUE 1\nBASIS(QQ)\nlagrange\n0\n0\n")
    bad[basis.Array[1]]("POLYNOMIALS-VALUE 1\nBASIS(QQ)\nbernstein\n0\n0\n")
    bad[basis.Array[1]]("POLYNOMIALS-VALUE 1\nBASIS(QQ)\nchebyshev\n1\n1\n1\n2\n")
    bad[basis.Array[1]]("POLYNOMIALS-VALUE 1\nBASIS(QQ)\nnewton\n2\n1\n2\n0\n")
    bad[batch.Array[Int]]("POLYNOMIALS-VALUE 1\nBATCH(Int)\n1\n2\n1\nInt\n0\n")
    bad[batch.Array[Integer]]("POLYNOMIALS-VALUE 1\nBATCH(Int)\n1\n0\n0\n")
    bad[batch.Array[Int]]("POLYNOMIALS-VALUE 1\nBATCH(Int)\n0\n0\n")
    bad[batch.Array[Int]]("POLYNOMIALS-VALUE 1\nBATCH(Int)\n1\n-1\n0\n")
    var sp = String("POLYNOMIALS-VALUE 1\nSPARSE(1,ZZ)\n")
    bad[sparse.Array[Int, 0]](String(sp, "1\n1\n1\n1\n2\n"))
    bad[sparse.Array[Int, 0]](String(sp, "1\n1\n0\n0\n2\n"))
    bad[sparse.Array[Int, 0]](String(sp, "1\n0\n0\n"))
    bad[sparse.Array[Int, 0]](String(sp, "2\n0\n1\n0\n2\n"))
    var mod = io.dumps(series.asarray([modint(1, 7)], order=2))
    bad[series.Array[3]](mod.replace("\n7\n", "\n9\n"))
    bad[series.Array[3]](mod.replace("2\n1\n1\n", "2\n1\n8\n"))
    var frac = String("POLYNOMIALS-VALUE 1\nRATIONAL(1,QQ)\n")
    var constant = String("SPARSE(1,QQ)\n1\n0\n")
    bad[rational.Function[Int]](String(frac, constant, "1\nSPARSE(1,QQ)\n0\n"))
    bad[rational.Function[Int]](String(frac, constant, "1\n", constant, "2\n"))
    var factor = String("SPARSE(1,QQ)\n1\n1\n0\n1\n1\n")
    bad[rational.Function[Int]](String(frac, factor, factor))
    # The defining polynomial and Sturm quotient must agree with the root.
    var irrational = String("0\n1\n2\n3\n-2\n0\n1\n3\n3\n-2\n0\n1\n2\n0\n2\n1\n2\n2\n2\n0\n1/2\n2\n0\n1\n")
    bad[algebraic.Real](String("POLYNOMIALS-VALUE 1\nAA\n", irrational.replace("1/2", "1/3")))
    bad[algebraic.Complex](String("POLYNOMIALS-VALUE 1\nQQbar\n", irrational, irrational.replace("1/2", "1/3")))
    var ctx = ArithmeticContext(format=FloatFormat(precision=80))
    var numeric = series.asarray([1], order=3, context=ctx)
    var zero_json = Float(0, context=ctx).to_json()
    bad[type_of(numeric)](io.dumps(numeric).replace(zero_json, Float(1, context=ctx).to_json()))
    print("value serialization rejection checks passed")
