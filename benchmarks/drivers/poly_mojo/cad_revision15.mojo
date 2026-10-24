"""Full-CAD projection-policy comparison, using dynamic variable counts.

Arguments: MODE INPUT ANSWER REPETITIONS. INPUT starts with the variable count,
then one polynomial per line in the independent test corpus grammar. Each
repetition includes full certificate checking; JSON records go to stdout.
"""
from polynomials import Budget, SparsePoly, QQ, RealAlgebraic, MonomialOrder, cad, check_cad
from polynomials.cad import _cad
from polynomials._cad_common import _coordinates
from support.corpus import parse_polynomial
from std.sys import argv
from std.time import perf_counter_ns


def main() raises:
    var args = argv()
    if len(args) != 5:
        raise Error("expected MODE INPUT ANSWER REPETITIONS")
    var mode = String(args[1])
    if mode != "lazard" and mode != "mccallum" and mode != "revision15" and mode != "revision15-finite":
        raise Error("unknown projection policy")
    var lines = open(String(args[2]), "r").read().splitlines()
    var dimensions = Int(lines[0])
    var work = Budget()
    var inputs = List[SparsePoly[QQ]]()
    var variables = List[Int64]()
    for i in range(dimensions):
        variables.append(Int64(i))
    for i in range(1, len(lines)):
        if lines[i]:
            inputs.append(parse_polynomial(String(lines[i]), MonomialOrder(), work))
    for repetition in range(Int(args[4])):
        var build_work = Budget()
        var start = perf_counter_ns()
        var result = _cad(inputs, variables, List[RealAlgebraic](), build_work, mode != "mccallum", brown=mode == "revision15" or mode == "revision15-finite", brown_finite=mode == "revision15-finite")
        var build_ns = perf_counter_ns() - start
        var check_work = Budget()
        start = perf_counter_ns()
        var valid = check_cad(inputs, result, check_work)
        var check_ns = perf_counter_ns() - start
        if not valid:
            raise Error("CAD certificate rejected")
        var terminal = 0
        var answer = String()
        for i in range(len(result.cells)):
            ref cell = result.cells[i]
            if cell.level == dimensions - 1:
                terminal += 1
                if repetition == 0:
                    var point = _coordinates(result.cells, i, variables, work)
                    answer.write(point[1], " ")
                    for j in range(len(cell.signs)):
                        answer.write("," if j else "", cell.signs[j].relation)
                    for entry in point[0]:
                        ref value = entry[1]
                        if value.rational:
                            answer.write(" Q:", value.low.value)
                        else:
                            answer.write(" A:")
                            for j in range(len(value.polynomial.coefficients)):
                                answer.write("," if j else "", value.polynomial.coefficients[j].value)
                            answer.write(":", value.low.value, ":", value.high.value)
                    answer.write("\n")
        if repetition == 0:
            with open(String(args[3]), "w") as output:
                output.write(String(terminal, "\n", answer))
        var factors = String("[")
        for i in range(len(result.projection.levels)):
            factors.write("," if i else "", len(result.projection.levels[i]))
        factors.write("]")
        print(String('{"repetition":', repetition, ',"construct_ns":', build_ns,
            ',"check_ns":', check_ns, ',"cells":', terminal,
            ',"stacks":', len(result.stacks), ',"factors":', factors,
            ',"check_ok":true}'))
