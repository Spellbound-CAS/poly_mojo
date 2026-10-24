"""Section sign witnesses: a dividing basis member's root proves a terminal zero, and forgeries fail."""
from polynomials import Budget, QQ, SparsePoly, MonomialOrder, cad, check_cad, CADCertificate
from support.check import Checker
from support.corpus import parse_polynomial


def first_witness(decomposition: CADCertificate, kind: Int, input: Int) -> Int:
    for i in range(len(decomposition.cells)):
        ref cell = decomposition.cells[i]
        if len(cell.witnesses) > input and cell.witnesses[input].kind == kind:
            return i
    return -1


def sector_cell(decomposition: CADCertificate, input: Int) -> Int:
    for i in range(len(decomposition.cells)):
        ref cell = decomposition.cells[i]
        if len(cell.witnesses) > input and cell.position % 2 == 0 and cell.signs[input].relation != 0:
            return i
    return -1


def main() raises:
    var work = Budget()
    var c = Checker("CAD section witnesses")
    var order = MonomialOrder()
    var circle = parse_polynomial("x0^2 + x1^2 - 1", order, work)
    var line = parse_polynomial("x1 - x0", order, work)
    var product = circle.mul(line, work)
    var inputs: List[SparsePoly[QQ]] = [circle.copy(), line.copy(), product.copy()]
    var variables: List[Int64] = [0, 1]
    for mode in range(2):
        var decomposition = cad(inputs, variables, work, lazard=Bool(mode))
        c.check(check_cad(inputs, decomposition, work), "decomposition with section witnesses checks")
        var sections = 0
        var zeros = 0
        for cell in decomposition.cells:
            for i in range(len(cell.witnesses)):
                if cell.witnesses[i].kind == 4:
                    sections += 1
                if cell.signs[i].relation == 0:
                    zeros += 1
        c.check(sections > 0 and sections == zeros, "every terminal zero is a section witness")
        # The product vanishes on the sections of both factors.
        var top = len(decomposition.projection.levels[1])
        c.check(top == 2, "two top-level basis members")
        var on_circle = first_witness(decomposition, 4, 0)
        var on_line = first_witness(decomposition, 4, 1)
        c.check(on_circle >= 0 and on_line >= 0, "both factors carry section witnesses")
        var circle_member = decomposition.cells[on_circle].witnesses[0].primitive
        var line_member = decomposition.cells[on_line].witnesses[1].primitive
        c.check(circle_member != line_member, "the factors are distinct basis members")
        c.check(decomposition.cells[on_circle].witnesses[2].kind == 4 and decomposition.cells[on_circle].witnesses[2].primitive == circle_member, "the product's witness names the circle on its sections")

        # A basis member that does not divide the input.
        var bad = decomposition.copy()
        bad.cells[on_circle].witnesses[0].primitive = line_member
        c.check(not check_cad(inputs, bad, work), "a non-dividing basis member is rejected")
        bad = decomposition.copy()
        bad.cells[on_circle].witnesses[0].primitive = top
        c.check(not check_cad(inputs, bad, work), "an out-of-range basis member is rejected")
        bad = decomposition.copy()
        bad.cells[on_circle].witnesses[0].primitive = -1
        c.check(not check_cad(inputs, bad, work), "a negative basis member is rejected")
        # A dividing member that has no root at this section: the circle at a
        # section of the line away from their intersections.
        var line_only = -1
        for i in range(len(decomposition.cells)):
            ref cell = decomposition.cells[i]
            if len(cell.witnesses) == 3 and cell.witnesses[1].kind == 4 and cell.signs[0].relation != 0:
                line_only = i
                break
        c.check(line_only >= 0, "a section of the line alone")
        bad = decomposition.copy()
        bad.cells[line_only].witnesses[0].kind = 4
        bad.cells[line_only].witnesses[0].primitive = circle_member
        bad.cells[line_only].witnesses[0].point = List[type_of(bad.cells[line_only].sample)]()
        bad.cells[line_only].witnesses[0].representations.clear()
        bad.cells[line_only].signs[0].relation = 0
        c.check(not check_cad(inputs, bad, work), "a section witness on another member's section is rejected")
        # A sector is no section.
        var sector = sector_cell(decomposition, 0)
        c.check(sector >= 0, "a sector cell")
        bad = decomposition.copy()
        bad.cells[sector].witnesses[0].kind = 4
        bad.cells[sector].witnesses[0].primitive = circle_member
        bad.cells[sector].witnesses[0].point = List[type_of(bad.cells[sector].sample)]()
        bad.cells[sector].witnesses[0].representations.clear()
        bad.cells[sector].signs[0].relation = 0
        c.check(not check_cad(inputs, bad, work), "a section witness on a sector is rejected")
        # The witness proves zero only.
        bad = decomposition.copy()
        bad.cells[on_circle].signs[0].relation = 1
        c.check(not check_cad(inputs, bad, work), "a nonzero sign with a section witness is rejected")
        bad = decomposition.copy()
        bad.cells[on_circle].witnesses[0].point.append(bad.cells[on_circle].sample.copy())
        c.check(not check_cad(inputs, bad, work), "a section witness with refined points is rejected")
    c.finish()
