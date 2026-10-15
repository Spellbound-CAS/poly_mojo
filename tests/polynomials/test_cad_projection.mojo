"""Projection closure, repeated/shared-factor bases and forged-proof rejection."""
from polynomials import Budget, MonomialOrder, SparsePoly, QQ, cad_projection, check_cad_projection
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var work = Budget()
    var c = Checker("CAD projection")
    var order = MonomialOrder()
    var circle = parse_polynomial("x0^2 + x1^2 - 1", order, work)
    var diagonal = parse_polynomial("x1 - x0", order, work)
    var inputs: List[SparsePoly[QQ]] = [circle.pow(2, work).mul(diagonal, work), diagonal.copy()]
    for mode in range(2):
        var proof = cad_projection(inputs, [Int64(0), Int64(1)], work, lazard=Bool(mode))
        c.check(check_cad_projection(inputs, proof, work), "complete McCallum/Lazard closure")
        c.check(len(proof.levels[1]) == 2, "shared and repeated factors become two coprime basis elements")
        var expected = parse_polynomial("x0^2 - 1/2", order, work)
        var found = False
        for basis in proof.levels[0]:
            found = found or basis == expected
        c.check(found, "circle-line intersection resultant retained")
        var bad = proof.copy()
        _ = bad.factorizations.pop()
        c.check(not check_cad_projection(inputs, bad, work), "missing projection factorization rejected")
        bad = proof.copy()
        bad.resultants[0].resultant = parse_polynomial("0", order, work)
        c.check(not check_cad_projection(inputs, bad, work), "false resultant identity rejected")
        bad = proof.copy()
        _ = bad.coprimality.pop()
        c.check(not check_cad_projection(inputs, bad, work), "missing basis coprimality rejected")
    c.finish()
