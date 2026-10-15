"""Selective overlapping coverings, early cylinders, and corrupted witnesses."""
from polynomials import Budget, QQ, SparsePoly, MonomialOrder, PolyErrorKind
from polynomials.cylindrical_covering import cylindrical_covering, check_cylindrical_covering
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var work = Budget()
    var c = Checker("selective cylindrical covering")
    var order = MonomialOrder()
    # x < 1 and x > 0 and x^2 >= 4 is impossible. The three conflict
    # intervals overlap and need fewer samples than all nine CAD cells.
    var inputs: List[SparsePoly[QQ]] = [parse_polynomial("x0 - 1", order, work), parse_polynomial("x0", order, work), parse_polynomial("x0^2 - 4", order, work)]
    var masks: List[UInt8] = [1, 4, 6]
    var variables: List[Int64] = [0]
    var proof = cylindrical_covering(inputs, masks, variables, work)
    c.check(proof.witness_cell == -1 and check_cylindrical_covering(inputs, masks, proof, work), "overlapping conflicts prove infeasibility")
    c.check(len(proof.conflicts) < len(proof.stacks[0].samples) and len(proof.cells) == 1, "conflicts skip irrelevant sections and sectors")
    var overlap = False
    for i in range(len(proof.conflicts)):
        for j in range(i):
            if max(proof.conflicts[i].lower, proof.conflicts[j].lower) <= min(proof.conflicts[i].upper, proof.conflicts[j].upper):
                overlap = True
    c.check(overlap, "overlapping intervals are retained and checked")
    var bad = proof.copy()
    _ = bad.conflicts.pop()
    c.check(not check_cylindrical_covering(inputs, masks, bad, work), "coverage gap rejected")
    bad = proof.copy()
    bad.conflicts[0].upper = len(proof.stacks[0].samples)
    c.check(not check_cylindrical_covering(inputs, masks, bad, work), "interval beyond stack rejected")
    bad = proof.copy()
    bad.conflicts[0].proof.sign.relation = 1 if bad.conflicts[0].proof.sign.relation == 0 else -bad.conflicts[0].proof.sign.relation
    c.check(not check_cylindrical_covering(inputs, masks, bad, work), "false conflict sign rejected")
    bad = proof.copy()
    bad.conflicts[0].proof.constraint = len(inputs)
    c.check(not check_cylindrical_covering(inputs, masks, bad, work), "out-of-range conflict constraint rejected")
    bad = proof.copy()
    bad.conflicts[0].lower = 0
    bad.conflicts[0].upper = len(proof.stacks[0].samples) - 1
    c.check(not check_cylindrical_covering(inputs, masks, bad, work), "expansion across a relevant root rejected")
    inputs = [parse_polynomial("x0^2 + 1", order, work), parse_polynomial("x1^2 + x2^2 - 1", order, work)]
    masks = [1, 2]
    variables = [0, 1, 2]
    proof = cylindrical_covering(inputs, masks, variables, work, lazard=True)
    c.check(proof.witness_cell == -1 and check_cylindrical_covering(inputs, masks, proof, work), "lower-dimensional conflict excludes entire cylinder")
    c.check(len(proof.stacks) == 1 and len(proof.cells) == 1, "higher coordinates never lifted after conflict")
    inputs = [parse_polynomial("x0^2 + x1^2 - 1", order, work)]
    masks = [1]
    variables = [0, 1]
    proof = cylindrical_covering(inputs, masks, variables, work)
    c.check(proof.witness_cell >= 0 and check_cylindrical_covering(inputs, masks, proof, work), "disk interior found selectively")
    c.check(len(proof.stacks) == 2, "zero-first sampling finds disk before other branches")
    inputs = [parse_polynomial("x0^2 - 2", order, work)]
    masks = [2]
    variables = [0]
    proof = cylindrical_covering(inputs, masks, variables, work)
    c.check(proof.witness_cell >= 0 and check_cylindrical_covering(inputs, masks, proof, work), "irrational satisfying section")
    bad = proof.copy()
    bad.cells[bad.witness_cell].signs[0].relation = 1
    c.check(not check_cylindrical_covering(inputs, masks, bad, work), "false satisfying sign rejected")
    inputs = [parse_polynomial("x0^2 + x1^2 - 1", order, work), parse_polynomial("x0 - 1", order, work)]
    masks = [2, 4]
    variables = [0, 1]
    proof = cylindrical_covering(inputs, masks, variables, work, lazard=True)
    c.check(proof.witness_cell == -1 and check_cylindrical_covering(inputs, masks, proof, work), "recursive coverage excludes circle beyond tangent")
    bad = proof.copy()
    if len(bad.cells) > 1:
        bad.cells[1].parent = 1
    c.check(not check_cylindrical_covering(inputs, masks, bad, work), "cyclic parent rejected")
    inputs = [parse_polynomial("1", order, work)]
    masks = [2]
    proof = cylindrical_covering(inputs, masks, variables, work)
    c.check(proof.root_conflict and len(proof.stacks) == 0 and check_cylindrical_covering(inputs, masks, proof, work), "constant contradiction needs no lifting")
    inputs = [parse_polynomial("x0", order, work)]
    masks = [0]
    proof = cylindrical_covering(inputs, masks, [Int64(0)], work)
    c.check(proof.witness_cell == -1 and len(proof.conflicts) == 1 and check_cylindrical_covering(inputs, masks, proof, work), "empty sign mask excludes the entire line")
    var invalid_mask = False
    try:
        _ = cylindrical_covering(inputs, [UInt8(8)], [Int64(0)], work)
    except e:
        invalid_mask = e.kind == PolyErrorKind.DOMAIN
    c.check(invalid_mask, "invalid mask rejected")
    inputs = [parse_polynomial("x0", order, work), parse_polynomial("x0*x1", order, work)]
    masks = [2, 4]
    proof = cylindrical_covering(inputs, masks, variables, work, lazard=True)
    c.check(proof.witness_cell == -1 and check_cylindrical_covering(inputs, masks, proof, work), "Lazard nullified input gives a valid selective conflict")
    inputs = List[SparsePoly[QQ]]()
    masks = List[UInt8]()
    proof = cylindrical_covering(inputs, masks, List[Int64](), work)
    c.check(proof.witness_cell == 0 and check_cylindrical_covering(inputs, masks, proof, work), "empty conjunction in zero variables")
    var exhausted = False
    var tiny = Budget(0)
    try:
        _ = check_cylindrical_covering(inputs, masks, proof, tiny)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "checking preserves budget failure")
    c.finish()
