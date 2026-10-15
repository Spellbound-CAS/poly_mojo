"""Exact calculus, composition, Taylor coefficients and simultaneous substitution."""
from apn_mojo import Integer, Rational, ExactComplex
from polynomials import Budget, PolyErrorKind, QQ, ZZ, QQi, ModInt, Modulus, NoContext, DensePoly, SparsePoly, MonomialOrder, derivative_dense, integrate_dense, integrate_integer_dense, definite_integral, compose_dense, affine_transform, taylor_shift, jet_evaluate, integrate_sparse, substitute_many, evaluate_dense, derivative
from support.check import Checker


def main() raises:
    var c = Checker("calculus")
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1)), QQ(Rational(2)), QQ(Rational(3))], work)
    var d = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(2)), QQ(Rational(6))], work)
    c.check(derivative_dense(f, work) == d, "dense derivative")
    var integral = integrate_dense(f, work)
    c.check(derivative_dense(integral, work) == f and integral.coefficients[0].is_zero(), "zero-constant antiderivative")
    c.check(definite_integral(f, QQ(Rational(0)), QQ(Rational(2)), work) == QQ(Rational(14)), "definite integral")
    var ints = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(0)), ZZ(Integer(1))], work)
    c.check(integrate_integer_dense(ints, work).coefficients[2] == QQ(Rational(1, 2)), "explicit integer promotion")
    var inner = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1)), QQ(Rational(2))], work)
    var expected = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(6)), QQ(Rational(16)), QQ(Rational(12))], work)
    c.check(compose_dense(f, inner, work) == expected, "hand-derived composition")
    c.check(affine_transform(f, QQ(Rational(2)), QQ(Rational(1)), work) == expected, "affine composition")
    var shifted = taylor_shift(f, QQ(Rational(1)), work)
    var expected_shift: List[QQ] = [QQ(Rational(6)), QQ(Rational(8)), QQ(Rational(3))]
    c.check(shifted.coefficients.to_list(work) == expected_shift, "Taylor shift")
    var jet = jet_evaluate(f, QQ(Rational(1)), 5, work)
    c.check(len(jet) == 5 and jet[0] == QQ(Rational(6)) and jet[1] == QQ(Rational(8)) and jet[2] == QQ(Rational(3)) and jet[3].is_zero() and jet[4].is_zero(), "Taylor jet with zero padding")
    c.check(len(jet_evaluate(f, QQ(Rational(1)), 0, work)) == 0, "empty jet")
    var order = MonomialOrder(MonomialOrder.LEX, [0, 1], work)
    var x = SparsePoly[QQ].variable(0, order, NoContext(), work)
    var y = SparsePoly[QQ].variable(1, order, NoContext(), work)
    var p = x.pow(2, work).add(y.scale(QQ(Rational(3)), work), work)
    var assignments: List[Tuple[Int64, SparsePoly[QQ]]] = [(Int64(0), y.copy()), (Int64(1), x.copy())]
    c.check(substitute_many(p, assignments, work) == y.pow(2, work).add(x.scale(QQ(Rational(3)), work), work), "simultaneous variable swap")
    c.check(derivative(integrate_sparse(p, 1, work), 1, work) == p, "sparse integration")
    var prime = Modulus(3, work)
    var modular = DensePoly[ModInt].from_coefficients(0, prime, [ModInt(0, prime), ModInt(0, prime), ModInt(1, prime)], work)
    try:
        _ = integrate_dense(modular, work)
        c.check(False, "characteristic-three obstruction ignored")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "nonintegrable finite-field coefficient")
    var z = QQi(ExactComplex(Rational(1, 2), Rational(2, 3)))
    var complex = DensePoly[QQi].from_coefficients(0, NoContext(), [z.copy(), QQi.one(NoContext())], work)
    c.check(derivative_dense(integrate_dense(complex, work), work) == complex, "Gaussian rational calculus")
    var tiny = Budget(0)
    var original = f.copy()
    try:
        _ = compose_dense(f, inner, tiny)
        c.check(False, "composition ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and f == original, "budget preserves input")
    c.finish()
