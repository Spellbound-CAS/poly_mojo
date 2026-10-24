"""Projection closure, repeated/shared-factor bases and forged-proof rejection."""
from polynomials import Budget, MonomialOrder, SparsePoly, QQ, cad_projection, check_cad_projection
from polynomials.cad_projection import _project, _empty
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
    # Reduced McCallum (revision 13): coefficients from the leading one down to
    # the first nonzero constant. Here x1 has the constant coefficient 1, so
    # x0^2 - 2 is not added, while the leading coefficient x0 is.
    var f = parse_polynomial("x0*x1^2 + x1 + x0^2 - 2", order, work)
    var single: List[SparsePoly[QQ]] = [f.copy()]
    var variables: List[Int64] = [Int64(0), Int64(1)]
    var reduced = cad_projection(single, variables, work)
    var leading = parse_polynomial("x0", order, work)
    var omitted = parse_polynomial("x0^2 - 2", order, work)
    var has_leading = False
    var has_omitted = False
    for basis in reduced.levels[0]:
        has_leading = has_leading or basis == leading
        has_omitted = has_omitted or basis == omitted
    c.check(has_leading and not has_omitted, "reduced McCallum stops at the first nonzero constant coefficient")
    c.check(check_cad_projection(single, reduced, work), "reduced McCallum projection checks")
    var full = _project[False](single, variables, False, _empty(variables, False, work), work, reduced=False)
    var full_has_omitted = False
    for basis in full.levels[0]:
        full_has_omitted = full_has_omitted or basis == omitted
    c.check(full_has_omitted and check_cad_projection(single, full, work), "pre-revision-13 full-coefficient projection still checks")
    var forged = reduced.copy()
    var dropped = False
    for i in range(len(forged.factorizations)):
        ref proof = forged.factorizations[i]
        if len(proof.factors) == 1 and proof.factors[0][0] == leading:
            _ = forged.factorizations.pop(i)
            dropped = True
            break
    c.check(dropped and not check_cad_projection(single, forged, work), "missing leading coefficient rejected")
    var constant = parse_polynomial("x1^2 + x0*x1 + x0^2 - 1", order, work)
    var monic: List[SparsePoly[QQ]] = [constant.copy()]
    var ldcf_only = cad_projection(monic, variables, work)
    var omitted_monic = False
    for basis in ldcf_only.levels[0]:
        omitted_monic = omitted_monic or basis == parse_polynomial("x0", order, work) or basis == parse_polynomial("x0^2 - 1", order, work)
    c.check(not omitted_monic and check_cad_projection(monic, ldcf_only, work), "constant leading coefficient adds no other coefficient")
    c.finish()
