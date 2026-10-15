"""Construct a complete certified cylindrical decomposition of the unit circle."""
from std.testing import assert_true
from polynomials import Budget, QQ, NoContext, SparsePoly, MonomialOrder, cad, check_cad, cylindrical_covering, check_cylindrical_covering


def main() raises:
    var budget = Budget()
    var order = MonomialOrder()
    var x = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var y = SparsePoly[QQ].variable(1, order, NoContext(), budget)
    var one = SparsePoly[QQ].constant(QQ.one(NoContext()), order, budget)
    var circle = x.pow(2, budget).add(y.pow(2, budget), budget).sub(one, budget)
    var inputs: List[SparsePoly[QQ]] = [circle.copy()]
    var result = cad(inputs, [Int64(0), Int64(1)], budget, lazard=True)
    assert_true(check_cad(inputs, result, budget))
    var plane_cells = 0
    for cell in result.cells:
        if cell.level == 1:
            plane_cells += 1
    assert_true(plane_cells == 13)
    var interior = cylindrical_covering(inputs, [UInt8(1)], [Int64(0), Int64(1)], budget, lazard=True)
    assert_true(interior.witness_cell >= 0 and check_cylindrical_covering(inputs, [UInt8(1)], interior, budget))
    assert_true(len(interior.stacks) == 2)
    print("Unit-circle plane cells:", plane_cells)
    print("Projection, lifting and coverage checked")
    print("Disk interior found with two selective stacks")
