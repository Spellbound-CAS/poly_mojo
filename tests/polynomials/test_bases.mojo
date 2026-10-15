"""Hand-derived basis coefficients, evaluations and subdivisions."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, BasisPolynomial, affine_transform, ModInt, Modulus
from support.check import Checker


def main() raises:
    var c = Checker("bases")
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1)), QQ(Rational(2)), QQ(Rational(3))], work)
    var nodes: List[QQ] = [QQ(Rational(1)), QQ(Rational(2))]
    var newton = BasisPolynomial[QQ].newton(f, nodes, work)
    var expected_newton: List[QQ] = [QQ(Rational(6)), QQ(Rational(11)), QQ(Rational(3))]
    c.check(newton.coefficients == expected_newton and newton.to_polynomial(work) == f, "Newton conversion coefficients")
    var chebyshev = BasisPolynomial[QQ].chebyshev(f, work)
    var expected_chebyshev: List[QQ] = [QQ(Rational(5, 2)), QQ(Rational(2)), QQ(Rational(3, 2))]
    c.check(chebyshev.coefficients == expected_chebyshev and chebyshev.to_polynomial(work) == f, "Chebyshev conversion coefficients")
    var bernstein = BasisPolynomial[QQ].bernstein(f, 2, work)
    var expected_bernstein: List[QQ] = [QQ(Rational(1)), QQ(Rational(2)), QQ(Rational(6))]
    c.check(bernstein.coefficients == expected_bernstein and bernstein.to_polynomial(work) == f, "Bernstein conversion coefficients")
    c.check(newton.evaluate(QQ(Rational(2)), work) == QQ(Rational(17)), "nested Newton evaluation")
    c.check(chebyshev.evaluate(QQ(Rational(2)), work) == QQ(Rational(17)), "Clenshaw evaluation")
    c.check(bernstein.evaluate(QQ(Rational(2)), work) == QQ(Rational(17)), "de Casteljau evaluation outside unit interval")
    var split = bernstein.subdivide(QQ(Rational(1, 2)), work)
    c.check(split[0].to_polynomial(work) == affine_transform(f, QQ(Rational(1, 2)), QQ(Rational(0)), work), "left subdivision identity")
    c.check(split[1].to_polynomial(work) == affine_transform(f, QQ(Rational(1, 2)), QQ(Rational(1, 2)), work), "right subdivision identity")
    var elevated = BasisPolynomial[QQ].bernstein(f, 4, work)
    var expected_elevated: List[QQ] = [QQ(Rational(1)), QQ(Rational(3, 2)), QQ(Rational(5, 2)), QQ(Rational(4)), QQ(Rational(6))]
    c.check(elevated.coefficients == expected_elevated and elevated.to_polynomial(work) == f, "Bernstein degree elevation")
    var zero = DensePoly[QQ](0, NoContext())
    c.check(len(BasisPolynomial[QQ].bernstein(zero, 5, work).coefficients) == 6, "zero retains Bernstein degree")
    var repeated: List[QQ] = [QQ(Rational(1)), QQ(Rational(1))]
    c.check(BasisPolynomial[QQ].newton(f, repeated, work).to_polynomial(work) == f, "repeated Newton nodes allowed")
    var malformed = chebyshev.copy()
    malformed.nodes.append(QQ(Rational(0)))
    try:
        _ = malformed.to_polynomial(work)
        c.check(False, "invalid basis shape accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "invalid basis shape rejected")
    var p = Modulus(2, work)
    var x2 = DensePoly[ModInt].from_coefficients(0, p, [ModInt(0, p), ModInt(0, p), ModInt(1, p)], work)
    try:
        _ = BasisPolynomial[ModInt].chebyshev(x2, work)
        c.check(False, "singular Chebyshev basis accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "characteristic-two conversion rejected")
    var half = QQ(Rational(1, 2))
    var tiny = Budget(0)
    try:
        _ = elevated.subdivide(half, tiny)
        c.check(False, "basis budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and elevated.to_polynomial(work) == f, "basis budget preserves input")
    c.finish()
