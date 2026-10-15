"""Named construction, safe grammar and exact expression round trips."""
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, PolyErrorKind, QQ, QQi, ZZ, NoContext, PolynomialRing, MonomialOrder, Modulus, ModInt
from support.check import Checker


def main() raises:
    var c = Checker("builders")
    var work = Budget()
    var names: List[String] = ["x", "y"]
    var ring = PolynomialRing[QQ].create(names, NoContext(), work)
    var x = ring.variable("x", work)
    var y = ring.variable("y", work)
    var f = ring.parse("(x + y)^2 - x^2 - y^2", work)
    c.check(f == x.mul(y, work).scale(QQ(Rational(2)), work), "binomial grammar expansion")
    var fractional = ring.parse("-x^2 + (3/2)*y - 2", work)
    c.check(ring.parse(ring.format(fractional, work), work) == fractional, "exact readable round trip")
    c.check(ring.parse("-x^2", work) == x.pow(2, work).neg(work), "power binds before unary minus")
    c.check(ring.parse("(-x)^2", work) == x.pow(2, work), "parentheses bind unary minus")
    var invalid: List[String] = ["", "x/(x+1)", "x^-1", "x^4294967296", "x +", "(x", "x)", "x y", "unknown", "x^2^3", "f(x)", "2.0*x"]
    for text in invalid:
        try:
            _ = ring.parse(text, work)
            c.check(False, "invalid expression accepted")
        except e:
            c.check(e.kind == PolyErrorKind.ENCODING or e.kind == PolyErrorKind.CAPACITY or e.kind == PolyErrorKind.DOMAIN, "bounded expression rejected")
    var gaussian = PolynomialRing[QQi].create(names, NoContext(), work)
    var complex = gaussian.parse("(1 + i)*x + (2 - 3*i)*y", work)
    c.check(gaussian.parse(gaussian.format(complex, work), work) == complex, "Gaussian rational expression round trip")
    var integer_ring = PolynomialRing[ZZ].create(names, NoContext(), work)
    c.check(integer_ring.parse("(x+x)/2", work) == integer_ring.variable("x", work), "exact integer coefficient division")
    var prime = Modulus(7, work)
    var modular = PolynomialRing[ModInt].create(names, prime, work)
    c.check(modular.parse("8*x", work) == modular.variable("x", work), "modular integer literals")
    var tiny = Budget(0)
    try:
        _ = ring.parse("x", tiny)
        c.check(False, "parser budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "parser budget")
    c.finish()
