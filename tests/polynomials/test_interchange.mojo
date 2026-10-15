"""Versioned exact and numerical interchange preserves representation metadata."""
from apn_mojo import Rational, Integer, ExactComplex, ArithmeticContext, ComplexContext, BallContext
from polynomials import Budget, PolyErrorKind, CoefficientRing, QQ, QQi, ZZ, ModInt, Modulus, NoContext, DensePoly, MonomialOrder, NumericDomain, FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain, NumericDensePoly, NumericSparsePoly, NumericBox, numeric_value, numeric_bounds, serialize_dense, deserialize_dense, serialize_sparse, deserialize_sparse, serialize_numeric_dense, deserialize_numeric_dense, serialize_numeric_sparse, deserialize_numeric_sparse
from support.check import Checker


def exact[C: CoefficientRing](ring: C.Context, coefficients: List[C], mut c: Checker, mut work: Budget) raises:
    var poly = DensePoly[C].from_coefficients(7, ring, coefficients.copy(), work)
    var text = serialize_dense(poly, work)
    c.check(deserialize_dense[C](text, ring, work) == poly, "exact dense interchange round trip")
    var order = MonomialOrder(MonomialOrder.GREVLEX, [9, 7], work)
    var sparse = poly.to_sparse(order, work)
    c.check(deserialize_sparse[C](serialize_sparse(sparse, work), ring, work) == sparse, "exact sparse interchange round trip")
    for n in [63, 64, 65, 127, 129]:
        var repeated = DensePoly[C].from_coefficients(7, ring, List[C](length=n, fill=coefficients[0]), work).to_sparse(order, work)
        c.check(deserialize_sparse[C](serialize_sparse(repeated, work), ring, work) == repeated, "bounded read batches preserve coefficients and ordering")
    var zero = DensePoly[C](7, ring.copy())
    c.check(deserialize_dense[C](serialize_dense(zero, work), ring, work) == zero, "zero retains dense variable and domain")
    var bad = String(text, "trailing")
    try:
        _ = deserialize_dense[C](bad, ring, work)
        c.check(False, "trailing interchange data accepted")
    except e:
        c.check(e.kind == PolyErrorKind.ENCODING, "trailing interchange data rejected")
    var tiny = Budget(0)
    try:
        _ = deserialize_dense[C](text, ring, tiny)
        c.check(False, "interchange ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "interchange budget")


def numeric[D: NumericDomain](context: D.Context, mut c: Checker, mut work: Budget) raises:
    var zero = QQ(Rational(0))
    var scalar = numeric_value[D](NumericBox.point(QQ(Rational(3, 2)), zero, work), context, work)
    var poly = NumericDensePoly[D].from_coefficients(11, [scalar.copy(), scalar.copy()], work)
    var text = serialize_numeric_dense(poly, work)
    var decoded = deserialize_numeric_dense[D](text, work)
    c.check(serialize_numeric_dense(decoded, work) == text and decoded.variable == 11, "numeric scalar formats preserved exactly")
    var order = MonomialOrder(MonomialOrder.LEX, [11, 12], work)
    var sparse = NumericSparsePoly[D].from_dense(poly, order, context, work)
    var sparse_text = serialize_numeric_sparse(sparse, work)
    var restored = deserialize_numeric_sparse[D](sparse_text, work)
    c.check(serialize_numeric_sparse(restored, work) == sparse_text, "numerical sparse metadata and scalars preserved")


def main() raises:
    var c = Checker("interchange")
    var work = Budget()
    var qq: List[QQ] = [QQ(Rational(-3, 7)), QQ(Rational(0)), QQ(Rational(5))]
    exact[QQ](NoContext(), qq, c, work)
    var zz: List[ZZ] = [ZZ(Integer(-9)), ZZ(Integer(1) << 256)]
    exact[ZZ](NoContext(), zz, c, work)
    var qi: List[QQi] = [QQi(ExactComplex(Rational(1, 3), Rational(-2, 5))), QQi(ExactComplex(0, 1))]
    exact[QQi](NoContext(), qi, c, work)
    var prime = Modulus(101, work)
    var mod: List[ModInt] = [ModInt(17, prime), ModInt(1, prime)]
    exact[ModInt](prime, mod, c, work)
    numeric[FloatDomain](ArithmeticContext(), c, work)
    numeric[ComplexDomain](ComplexContext(ArithmeticContext()), c, work)
    numeric[BallDomain](BallContext(128), c, work)
    numeric[ComplexBallDomain](BallContext(128), c, work)
    var negative = DensePoly[QQ].from_coefficients(-17, NoContext(), [QQ(Rational(1))], work)
    c.check(deserialize_dense[QQ](serialize_dense(negative, work), NoContext(), work) == negative, "signed variable IDs retained")
    var malformed: List[String] = ["POLYNOMIALS 2\nDENSE\nQQ\n0\n0\n", "POLYNOMIALS 1\nDENSE\nQQ\n00\n0\n", "POLYNOMIALS 1\nDENSE\nQQ\n0\n1\n2/4\n", "POLYNOMIALS 1\nDENSE\nQQ\n0\n1\n0\n"]
    for text in malformed:
        try:
            _ = deserialize_dense[QQ](text, NoContext(), work)
            c.check(False, "noncanonical interchange accepted")
        except e:
            c.check(e.kind == PolyErrorKind.ENCODING, "noncanonical interchange rejected")
    for corrupt in [0, 1]:
        var text = String("POLYNOMIALS 1\nSPARSE\nQQ\n1\n1\n0\n65\n")
        for i in range(65):
            text.write("1\n0\n", i + 1 - Int(corrupt == 0 and i == 64), "\n", "0" if corrupt == 1 and i == 64 else "1", "\n")
        try:
            _ = deserialize_sparse[QQ](text, NoContext(), work)
            c.check(False, "malformed term at batch boundary accepted")
        except e:
            c.check(e.kind == PolyErrorKind.ENCODING, "batch boundary rejects duplicate monomials and zero coefficients")
    c.finish()
