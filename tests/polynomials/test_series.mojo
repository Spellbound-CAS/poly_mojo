"""Formal series identities, precision boundaries and finite characteristic."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, TruncatedSeries, Modulus, ModInt
from support.check import Checker


def series(values: List[Int], order: Int, mut work: Budget) raises -> TruncatedSeries[QQ]:
    var coefficients = List[QQ]()
    for value in values:
        coefficients.append(QQ(Rational(value)))
    return TruncatedSeries[QQ].from_polynomial(DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work), order, work)


def main() raises:
    var c = Checker("series")
    var work = Budget()
    var one = series([1], 7, work)
    var geometric = series([1, -1], 7, work).inverse(work)
    c.check(geometric.polynomial == series([1, 1, 1, 1, 1, 1, 1], 7, work).polynomial, "geometric series coefficients")
    c.check(geometric.mul(series([1, -1], 7, work), work).polynomial == one.polynomial, "inverse identity modulo x^7")
    var f = series([0, 1, 1], 7, work)
    var reverted = f.revert(work)
    c.check(reverted.polynomial == series([0, 1, -1, 2, -5, 14, -42], 7, work).polynomial, "Catalan reversion coefficients")
    c.check(f.compose(reverted, work).polynomial == series([0, 1], 7, work).polynomial, "compositional inverse identity")
    c.check(one.add(f, work).log(work).exp(work).polynomial == one.add(f, work).polynomial, "exp log identity")
    c.check(f.exp(work).log(work).polynomial == f.polynomial, "log exp identity")
    c.check(f.derivative(work).integrate(work).polynomial == f.polynomial, "series calculus preserves known precision")
    var tiny = Budget(0)
    try:
        _ = f.revert(tiny)
        c.check(False, "reversion ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and f.order == 7, "series budget preserves input")
    try:
        _ = f.compose(one, work)
        c.check(False, "nonzero inner constant accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "unknown series tail cannot compose at nonzero constant")
    try:
        _ = f.truncate(8, work)
        c.check(False, "series invented extra precision")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "precision cannot increase")
    var p = Modulus(3, work)
    var finite = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, p, [ModInt(0, p), ModInt(1, p)], work), 4, work)
    try:
        _ = finite.exp(work)
        c.check(False, "finite-characteristic exponential accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "noninvertible exponential denominator rejected")
    c.finish()
