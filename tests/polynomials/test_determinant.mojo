"""Independent determinant values and corruption tests for CAD projection."""
from polynomials import Budget, MonomialOrder, projection_resultant, check_projection_resultant
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var work = Budget()
    var c = Checker("projection determinants")
    var order = MonomialOrder()
    var f = parse_polynomial("x1^2 - x0", order, work)
    var g = parse_polynomial("2*x1", order, work)
    var expected = parse_polynomial("-4*x0", order, work)
    var proof = projection_resultant(f, g, 1, work)
    c.check(proof.resultant == expected, "Res(y^2-x,2y)=-4x")
    c.check(check_projection_resultant(f, g, proof, work), "fraction-free determinant proof")
    var bad = proof.copy()
    bad.steps[0].trailing[0][0] = parse_polynomial("17", order, work)
    c.check(not check_projection_resultant(f, g, bad, work), "corrupted trailing matrix rejected")
    bad = proof.copy()
    bad.resultant = expected.neg(work)
    c.check(not check_projection_resultant(f, g, bad, work), "wrong determinant sign rejected")
    var zero = projection_resultant(f, f, 1, work)
    c.check(zero.resultant.is_zero() and check_projection_resultant(f, f, zero, work), "singular Sylvester matrix")
    var linear = parse_polynomial("x1 + x0", order, work)
    var higher = parse_polynomial("x1^3 + x0", order, work)
    var value = parse_polynomial("-x0^3 + x0", order, work)
    var asymmetric = projection_resultant(linear, higher, 1, work)
    c.check(asymmetric.resultant == value and check_projection_resultant(linear, higher, asymmetric, work), "linear-cubic evaluation determinant")
    c.finish()
