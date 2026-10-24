"""Compare historical full CAD using sources frozen at the original commit.

Arguments: INPUT ANSWER. INPUT is the variable count followed by corpus-grammar
polynomials. Uses the explicit plain-Lazard reference policy, then its
complete checker. Phase records survive a timeout; only the final JSON record
denotes an accepted result. Answer aggregation is outside both timed regions.
"""
from polynomials import Budget, SparsePoly, QQ, MonomialOrder, cad, check_cad
from support.corpus import parse_polynomial
from std.sys import argv
from std.time import perf_counter_ns


def progress(path: String, stage: String, elapsed: Int, stacks: Int) raises:
    with open(path + ".progress", "a") as output:
        output.write(String('{"stage":"', stage, '","elapsed_ns":', elapsed,
            ',"stacks":', stacks, '}\n'))


def main() raises:
    var args = argv()
    if len(args) != 3:
        raise Error("expected INPUT ANSWER")
    var answer = String(args[2])
    with open(answer + ".progress", "w") as output:
        output.write("")
    var lines = open(String(args[1]), "r").read().splitlines()
    var dimensions = Int(lines[0])
    var work = Budget()
    var inputs = List[SparsePoly[QQ]]()
    var variables = List[Int64]()
    for i in range(dimensions):
        variables.append(Int64(i))
    for i in range(1, len(lines)):
        if lines[i]:
            inputs.append(parse_polynomial(String(lines[i]), MonomialOrder(), work))
    progress(answer, "construct_started", 0, 0)
    var build_work = Budget()
    var start = perf_counter_ns()
    var result = cad(inputs, variables, build_work, lazard=True)
    var construct_ns = perf_counter_ns() - start
    progress(answer, "constructed", construct_ns, len(result.stacks))
    var check_work = Budget()
    start = perf_counter_ns()
    var valid = check_cad(inputs, result, check_work)
    var check_ns = perf_counter_ns() - start
    if not valid:
        raise Error("historical CAD certificate rejected")
    progress(answer, "verified", check_ns, len(result.stacks))
    var seen = Dict[String, Bool]()
    var cells = 0
    with open(answer, "w") as output:
        for cell in result.cells:
            if cell.level != dimensions - 1:
                continue
            cells += 1
            var row = String()
            for j in range(len(cell.signs)):
                row.write("," if j else "", cell.signs[j].relation)
            if row not in seen:
                seen[row] = True
                output.write(String(row, "\n"))
    var factors = String("[")
    for i in range(len(result.projection.levels)):
        factors.write("," if i else "", len(result.projection.levels[i]))
    factors.write("]")
    print(String('{"construct_ns":', construct_ns, ',"check_ns":', check_ns,
        ',"cells":', cells, ',"stacks":', len(result.stacks),
        ',"factors":', factors, ',"check_ok":true}'))
