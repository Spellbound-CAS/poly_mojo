"""Exact multiplicities require product evidence, coprime blocks and a vanishing block per root."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, NumericOptions, PolyErrorKind, real_roots_exact, real_root_multiplicities, check_real_root_multiplicities, check_squarefree_product
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
    var expected: List[Int] = [3, 2]
    c.check(certificate.multiplicities == expected, "known repeated-root multiplicities in root order")
    var bad = certificate.copy()
    bad.multiplicities[0] += 1
    c.check(not check_real_root_multiplicities(f, roots, options, bad, work), "forged multiplicity rejected")
    bad = certificate.copy()
    bad.blocks[0] = bad.blocks[1]
    bad.multiplicities[0] = bad.multiplicities[1]
    c.check(not check_real_root_multiplicities(f, roots, options, bad, work), "block that does not vanish at the root rejected")
    bad = certificate.copy()
    _ = bad.blocks.pop()
    c.check(not check_real_root_multiplicities(f, roots, options, bad, work), "missing root block rejected")
    # (x - 1)(x + 2) squared times (x + 2) is f, but the blocks share a root,
    # so their product is not the squarefree part.
    bad = certificate.copy()
    var x_minus_1 = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-1)), QQ(Rational(1))], work)
    var x_plus_2 = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(2)), QQ(Rational(1))], work)
    bad.decomposition.factors = [(x_minus_1.mul(x_plus_2, work), 2), (x_plus_2.copy(), 1)]
    c.check(check_squarefree_product(f, bad.decomposition, work), "the overlapping blocks still multiply to f")
    bad.blocks = [1, 0]
    bad.multiplicities = [1, 2]
    c.check(not check_real_root_multiplicities(f, roots, options, bad, work), "blocks sharing a root rejected")
    var tiny = Budget(0)
    try:
        _ = real_root_multiplicities(f, roots, options, tiny)
        c.check(False, "multiplicity budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "multiplicity budget")
    c.finish()
