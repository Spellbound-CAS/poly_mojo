"""Exact multiplicities require product, squarefreeness and root-count evidence."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, NumericOptions, PolyErrorKind, real_roots_exact, real_root_multiplicities, check_real_root_multiplicities
from support.check import Checker


def main() raises:
    var c = Checker("root_multiplicity")
    var work = Budget()
    # (x+2)^3*(x-1)^2 = x^5+4x^4+x^3-10x^2-4x+8.
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(8)), QQ(Rational(-4)), QQ(Rational(-10)), QQ(Rational(1)), QQ(Rational(4)), QQ(Rational(1))], work)
    var options = NumericOptions(absolute_bits=3, max_subdivisions=256)
    var roots = real_roots_exact(f, options, work, complete_only=True)
    var certificate = real_root_multiplicities(f, roots, options, work)
    c.check(check_real_root_multiplicities(f, roots, options, certificate, work), "complete multiplicity certificate")
    var multiplicities = List[Int]()
    for value in certificate.multiplicities:
        if value:
            multiplicities.append(value)
    var expected: List[Int] = [3, 2]
    c.check(multiplicities == expected, "known repeated-root multiplicities")
    var bad = certificate.copy()
    bad.multiplicities[0] += 1
    c.check(not check_real_root_multiplicities(f, roots, options, bad, work), "forged multiplicity rejected")
    bad = certificate.copy()
    _ = bad.sequences.pop()
    c.check(not check_real_root_multiplicities(f, roots, options, bad, work), "missing factor Sturm evidence rejected")
    bad = certificate.copy()
    var zero = DensePoly[QQ](0, NoContext())
    bad.squarefree_witnesses[0] = (zero.copy(), zero.copy())
    c.check(not check_real_root_multiplicities(f, roots, options, bad, work), "false squarefreeness witness rejected")
    var tiny = Budget(0)
    try:
        _ = real_root_multiplicities(f, roots, options, tiny)
        c.check(False, "multiplicity budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "multiplicity budget")
    c.finish()
