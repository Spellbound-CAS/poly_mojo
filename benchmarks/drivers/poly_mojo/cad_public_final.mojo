"""Time the accepted public core result, without eagerly creating output cells.

Arguments: MODE INPUT ANSWER REPETITIONS. accepted_ns includes construction,
final-data assembly and full SIMD acceptance. recheck_scalar_ns is an extra
independent check, performed once and never added to the accepted-call cost.
Private data access below is untimed benchmark instrumentation only.
"""
from polynomials import Budget, SparsePoly, QQ, MonomialOrder, cad_final, check_cad_final
from support.corpus import parse_polynomial
from std.sys import argv
from std.time import perf_counter_ns


def _json_ints(values: List[Int]) -> String:
    var text = String("[")
    for i in range(len(values)):
        text.write("," if i else "", values[i])
    text.write("]")
    return text^


def _signature(signs: List[Int8], count: Int) -> String:
    var text = String()
    for i in range(count):
        text.write("," if i else "", signs[len(signs) - count + i])
    return text^


def main() raises:
    var args = argv()
    if len(args) != 5:
        raise Error("expected MODE INPUT ANSWER REPETITIONS")
    var mode = String(args[1])
    if mode != "final" and mode != "lazard" and mode != "mccallum":
        raise Error("unknown mode")
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
    var lazard = Optional[Bool]() if mode == "final" else Optional(mode == "lazard")
    with open(String(args[3]) + ".progress", "w") as output:
        output.write("")
    for repetition in range(Int(args[4])):
        var build_work = Budget()
        var start = perf_counter_ns()
        var proof = cad_final(inputs, variables, build_work, lazard=lazard)
        var accepted_ns = perf_counter_ns() - start
        with open(String(args[3]) + ".progress", "a") as output:
            output.write(String('{"stage":"accepted","repetition":', repetition,
                ',"elapsed_ns":', accepted_ns, '}\n'))
        var recheck_ns = 0
        if repetition == 0:
            var check_work = Budget()
            start = perf_counter_ns()
            if not check_cad_final(inputs, variables, List[UInt8](), proof, check_work, vectorized=False):
                raise Error("public CAD rejected by scalar recheck")
            recheck_ns = perf_counter_ns() - start
        var cells = List[Int](length=dimensions, fill=0)
        var factors = proof.factor_counts(work)
        var signatures = Dict[String, Bool]()
        var seen = Dict[Int, Bool]()
        var encodable = True
        var possibilities = 1
        for _ in inputs:
            if possibilities > Int.MAX // 3:
                encodable = False
                break
            possibilities *= 3
        var native_nodes = 0
        var wide_nodes = 0
        var exact_fallbacks = 0
        var materializations = 0
        var field_map_hits = 0
        var max_field_degree = 1
        var retries = 0
        var reduced_projection = False
        if proof.is_compact():
            ref tree = proof._data[].compact.value()
            retries = tree.retries
            reduced_projection = tree.projection.brown
            materializations = tree.materializations
            field_map_hits = tree.field_map_hits
            for stack in tree.stacks:
                cells[stack.level] += len(stack.stack.signs)
                max_field_degree = max(max_field_degree, stack.field.modulus.degree())
                native_nodes += stack.stack.isolation.native_nodes
                wide_nodes += stack.stack.isolation.wide_nodes
                exact_fallbacks += stack.stack.isolation.exact_fallbacks
                if repetition == 0 and stack.level == dimensions - 1:
                    for row in stack.stack.signs:
                        if encodable:
                            var code = 0
                            for j in range(len(row) - len(inputs), len(row)):
                                code = 3 * code + Int(row[j]) + 1
                            if code in seen:
                                continue
                            seen[code] = True
                        signatures[_signature(row, len(inputs))] = True
        else:
            ref legacy = proof.legacy_decomposition().value()
            for cell in legacy.cells:
                if cell.level >= 0:
                    cells[cell.level] += 1
                if repetition == 0 and cell.level == dimensions - 1:
                    var text = String()
                    for i in range(len(cell.signs)):
                        text.write("," if i else "", cell.signs[i].relation)
                    signatures[text^] = True
        if repetition == 0:
            with open(String(args[3]), "w") as output:
                for row in signatures:
                    output.write(String(row, "\n"))
        print(String('{"repetition":', repetition, ',"accepted_ns":', accepted_ns,
            ',"recheck_scalar_ns":', recheck_ns, ',"cells_by_level":', _json_ints(cells),
            ',"cells":', cells[dimensions - 1], ',"stacks":', proof.stack_count(),
            ',"factors":', _json_ints(factors), ',"compact":', "true" if proof.is_compact() else "false",
            ',"materializations":', materializations, ',"field_map_hits":', field_map_hits,
            ',"max_field_degree":', max_field_degree,
            ',"construction_retries":', retries, ',"revision15":', "true" if reduced_projection else "false",
            ',"native_nodes":', native_nodes, ',"wide_nodes":', wide_nodes,
            ',"exact_fallbacks":', exact_fallbacks, ',"check_ok":true}'))
