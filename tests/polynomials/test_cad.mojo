"""Hand-derived CAD cells, nullification, exact samples and coverage corruption."""
from apn_mojo import Rational
from polynomials import Budget, QQ, DensePoly, SparsePoly, RealAlgebraic, MonomialOrder, real_rational, real_root, cad, check_cad, cad_cell, check_cad_cell, cad_covering, check_cad_covering, PolyErrorKind
from support.check import Checker
from support.corpus import parse_polynomial
from polynomials import cad_projection
from polynomials._cad_lifting import _stack, _check_stack, _position


def q(value: Int, mut work: Budget) raises -> RealAlgebraic:
    return real_rational(QQ(Rational(value)), work)


def main() raises:
    var work = Budget()
    var c = Checker("CAD")
    var order = MonomialOrder()
    var circle = parse_polynomial("x0^2 + x1^2 - 1", order, work)
    var inputs: List[SparsePoly[QQ]] = [circle.copy()]
    for mode in range(2):
        var decomposition = cad(inputs, [Int64(0), Int64(1)], work, lazard=Bool(mode))
        c.check(check_cad(inputs, decomposition, work), "certified unit-circle decomposition")
        c.check(len(decomposition.cells) == 19 and len(decomposition.stacks) == 6, "five base cells and thirteen plane cells")
        var negative = 0
        var zero = 0
        var positive = 0
        for cell in decomposition.cells:
            if cell.level == 1:
                if cell.signs[0].relation < 0:
                    negative += 1
                elif cell.signs[0].relation == 0:
                    zero += 1
                else:
                    positive += 1
        c.check(negative == 1 and zero == 4 and positive == 8, "hand-derived signs across all thirteen cells")
        if not mode:
            var bad = decomposition.copy()
            _ = bad.cells.pop()
            c.check(not check_cad(inputs, bad, work), "missing sector breaks coverage")
            bad = decomposition.copy()
            bad.cells[6].parent = 0
            c.check(not check_cad(inputs, bad, work), "wrong parent breaks cylindricity")
            bad = decomposition.copy()
            bad.stacks[0].sources[0] = (99, 0)
            c.check(not check_cad(inputs, bad, work), "forged root-function source rejected")
            bad = decomposition.copy()
            bad.cells[len(bad.cells) - 1].signs[0].relation = -1
            c.check(not check_cad(inputs, bad, work), "false global cell sign rejected")
    var local_point: List[RealAlgebraic] = [real_rational(QQ(Rational(1, 2)), work), q(0, work)]
    var local = cad_cell(inputs, [Int64(0), Int64(1)], local_point, work)
    c.check(check_cad_cell(inputs, local_point, local, work) and local.signs[0].relation == -1, "single cell uses root functions at the query base")
    c.check(len(local.stacks) == 2, "single-cell construction lifts only its own path")
    var wrong_local = local.copy()
    wrong_local.positions[1] = 0
    c.check(not check_cad_cell(inputs, local_point, wrong_local, work), "point outside selected sector rejected")
    var nullified = parse_polynomial("x0*x2 + x1*x3", order, work)
    var null_inputs: List[SparsePoly[QQ]] = [nullified.copy()]
    var null_point: List[RealAlgebraic] = [q(0, work), q(0, work), q(1, work), q(1, work)]
    var variables: List[Int64] = [0, 1, 2, 3]
    var rejected = False
    try:
        _ = cad_cell(null_inputs, variables, null_point, work)
    except e:
        rejected = e.kind == PolyErrorKind.UNSUPPORTED
    c.check(rejected, "positive-dimensional McCallum nullification fails closed")
    var lazard = cad_cell(null_inputs, variables, null_point, work, lazard=True)
    c.check(check_cad_cell(null_inputs, null_point, lazard, work) and lazard.signs[0].relation == 0, "Lazard handles the same nullification exactly")
    var positive_poly = parse_polynomial("x0^2 + 1", order, work)
    var impossible: List[SparsePoly[QQ]] = [positive_poly.copy()]
    var covering = cad_covering(impossible, [UInt8(1)], [Int64(0)], work)
    c.check(covering.witness_cell == -1 and check_cad_covering(impossible, [UInt8(1)], covering, work), "x^2+1<0 has a complete conflict covering")
    var bad_covering = covering.copy()
    _ = bad_covering.violations.pop()
    c.check(not check_cad_covering(impossible, [UInt8(1)], bad_covering, work), "incomplete conflict covering rejected")
    var equation = parse_polynomial("x0^2 - 2", order, work)
    var feasible: List[SparsePoly[QQ]] = [equation.copy()]
    var solution = cad_covering(feasible, [UInt8(2)], [Int64(0)], work)
    c.check(solution.witness_cell >= 0 and check_cad_covering(feasible, [UInt8(2)], solution, work), "irrational satisfying section is certified")
    var tiny = Budget(1)
    var exhausted = False
    try:
        _ = check_cad_cell(inputs, local_point, local, tiny)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "CAD verification preserves LIMIT")
    var shared_inputs = List[SparsePoly[QQ]]()
    for i in range(6):
        shared_inputs.append(parse_polynomial(String("x1^2 - 2 - ", i + 1, "*x0"), order, work))
    var projection = cad_projection(shared_inputs, [Int64(0), Int64(1)], work)
    var base: List[Tuple[Int64, RealAlgebraic]] = [(Int64(0), q(0, work))]
    var shared = _stack(projection, 0, 1, base, 0, None, work)
    c.check(len(shared.roots) == 2 and len(shared.samples) == 5, "six specialized quadratics share just two boundaries")
    c.check(_check_stack(projection, 0, 1, base, 0, shared, work), "shared-root merging has independent evidence")
    for indices in shared.mapping:
        c.check(len(indices) == 2 and indices[0] == 0 and indices[1] == 1, "every specialized root maps to its ordered boundary")
    for i in range(len(shared.samples)):
        c.check(_position(shared.roots, shared.samples[i], work) == i, "boundary and sector position")
    c.check(shared.sources[0] == (0, 0) and shared.sources[1] == (0, 1), "first source wins deterministic equality ties")
    shared.mapping[5][1] = 0
    c.check(not _check_stack(projection, 0, 1, base, 0, shared, work), "forged shared-root mapping rejected")
    var distinct_inputs = List[SparsePoly[QQ]]()
    for i in range(5):
        distinct_inputs.append(parse_polynomial(String("x1^2 - ", i + 2, " - ", i + 1, "*x0"), order, work))
    var distinct_projection = cad_projection(distinct_inputs, [Int64(0), Int64(1)], work)
    var distinct_stack = _stack(distinct_projection, 0, 1, base, 0, None, work)
    c.check(len(distinct_stack.roots) == 10 and _position(distinct_stack.roots, q(0, work), work) == 10, "five positive and five negative quadratic roots bracket zero")
    c.check(_check_stack(distinct_projection, 0, 1, base, 0, distinct_stack, work), "odd run count merges rational and irrational boundaries")
    c.finish()
