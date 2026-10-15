"""Independent existential formulas shared with the pinned QEPCAD B cross-check."""
from polynomials import Budget, QQ, SparsePoly, MonomialOrder, cad_covering, check_cad_covering
from polynomials.cylindrical_covering import cylindrical_covering, check_cylindrical_covering
from support.corpus import parse_polynomial
from support.check import Checker


def verify(inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], expected: Bool, label: String, mut c: Checker) raises:
    var budget = Budget()
    var result = cad_covering(inputs, masks, variables, budget, lazard=True)
    c.check((result.witness_cell >= 0) == expected and check_cad_covering(inputs, masks, result, budget), label)
    var selective = cylindrical_covering(inputs, masks, variables, budget, lazard=True)
    c.check((selective.witness_cell >= 0) == expected and check_cylindrical_covering(inputs, masks, selective, budget), label + " (selective)")


def main() raises:
    var work = Budget()
    var c = Checker("CAD independent formulas")
    var order = MonomialOrder()
    var positive = parse_polynomial("x0^2 + 1", order, work)
    var quadratic = parse_polynomial("x0^2 - 2", order, work)
    var circle = parse_polynomial("x0^2 + x1^2 - 1", order, work)
    var diagonal = parse_polynomial("x1 - x0", order, work)
    var tangent = parse_polynomial("x0 - 1", order, work)
    verify([positive.copy()], [UInt8(1)], [Int64(0)], False, "positive quadratic cannot be negative", c)
    verify([quadratic.copy()], [UInt8(2)], [Int64(0)], True, "irrational real zero exists", c)
    verify([circle.copy()], [UInt8(1)], [Int64(0), Int64(1)], True, "unit disk has an interior", c)
    verify([circle.copy(), diagonal.copy()], [UInt8(2), UInt8(2)], [Int64(0), Int64(1)], True, "circle intersects the diagonal", c)
    verify([circle.copy(), tangent.copy()], [UInt8(2), UInt8(4)], [Int64(0), Int64(1)], False, "circle cannot lie beyond its tangent", c)
    verify([circle.copy(), circle.copy()], [UInt8(1), UInt8(4)], [Int64(0), Int64(1)], False, "opposite strict signs cannot hold together", c)
    c.finish()
