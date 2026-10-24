"""Full private CAD: separately time construction and exact final verification.

Arguments: MODE INPUT ANSWER REPETITIONS. No terminal root is materialized for
this benchmark. ANSWER contains the realized original-input sign vectors; the
complete final checker establishes root completeness and cell-tree coverage.
"""
from polynomials import Budget, SparsePoly, QQ, MonomialOrder
from polynomials._cad_tree_final import _construct_final_cad, _check_final_cad_with_cache, _terminal_signs
from polynomials._cad_root_cache import _RootCountCache
from support.corpus import parse_polynomial
from std.sys import argv
from std.time import perf_counter_ns


def _progress(path: String, stage: String, repetition: Int, elapsed: Int, stacks: Int) raises:
    # Diagnostics survive a timeout. Only the final stdout record denotes a
    # result accepted by both complete verifiers.
    with open(path + ".progress", "a") as output:
        output.write(String('{"stage":"', stage, '","repetition":', repetition,
            ',"elapsed_ns":', elapsed, ',"stacks":', stacks, '}\n'))


def main() raises:
    var args = argv()
    if len(args) != 5:
        raise Error("expected MODE INPUT ANSWER REPETITIONS")
    var mode = String(args[1])
    if mode != "lazard" and mode != "revision15" and mode != "interval":
        raise Error("unknown projection policy")
    with open(String(args[3]) + ".progress", "w") as progress:
        progress.write("")
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
        var result = _construct_final_cad(inputs, variables, List[Int](), build_work, brown=mode != "lazard", interval=mode == "interval")
        var build_ns = perf_counter_ns() - start
        _progress(String(args[3]), "constructed", repetition, build_ns, len(result.stacks))
        var verify_scalar_ns = 0
        var verify_simd_ns = 0
        var check_stats = List[Int](length=7, fill=0)
        # Alternate order to reduce warm-cache bias. Each invocation builds its
        # own checked-fact cache and verifies the complete original claim.
        for trial in range(2):
            var vectorized = (repetition + trial) % 2 == 0
            var check_work = Budget()
            start = perf_counter_ns()
            var root_counts = _RootCountCache()
            var valid = _check_final_cad_with_cache(inputs, variables, List[Int](), result, root_counts, check_work, vectorized=vectorized)
            var elapsed = perf_counter_ns() - start
            if not valid:
                raise Error("final CAD rejected")
            _progress(String(args[3]), "verified_simd" if vectorized else "verified_scalar", repetition, elapsed, len(result.stacks))
            if vectorized:
                verify_simd_ns = elapsed
                check_stats = [root_counts.hits, root_counts.misses, root_counts.linear_checks, root_counts.flushes, root_counts.stats.word_points, root_counts.stats.big_points, root_counts.stats.vector_batches]
            else:
                verify_scalar_ns = elapsed
        var cells = 0
        var native_nodes = 0
        var wide_nodes = 0
        var exact_fallbacks = 0
        var rational_short_stacks = 0
        var rational_vector_stacks = 0
        var signatures = List[String]()
        var seen = Dict[Int, Bool]()
        var seen_text = Dict[String, Bool]()
        var encodable = True
        var possibilities = 1
        for _ in inputs:
            if possibilities > Int.MAX // 3:
                encodable = False
                break
            possibilities *= 3
        for stack in result.stacks:
            if not stack.field.embedding:
                if len(stack.stack.roots) < 4:
                    rational_short_stacks += 1
                else:
                    rational_vector_stacks += 1
            native_nodes += stack.stack.isolation.native_nodes
            wide_nodes += stack.stack.isolation.wide_nodes
            exact_fallbacks += stack.stack.isolation.exact_fallbacks
            if stack.level == dimensions - 1:
                cells += len(stack.stack.signs)
                if repetition == 0:
                    for position in range(len(stack.stack.signs)):
                        # The complete verifier has already checked every
                        # cell. Aggregate its realized sign vectors without
                        # a linear scan of all previously seen vectors per
                        # cell, or formatting duplicate vectors as strings.
                        if encodable:
                            ref stored = stack.stack.signs[position]
                            var code = 0
                            for j in range(len(stored) - len(inputs), len(stored)):
                                code = 3 * code + Int(stored[j]) + 1
                            if code in seen:
                                continue
                            seen[code] = True
                        var row = _terminal_signs(stack, len(inputs), position, work)
                        var text = String()
                        for j in range(len(row)):
                            text.write("," if j else "", row[j])
                        if not encodable:
                            if text in seen_text:
                                continue
                            seen_text[text] = True
                        signatures.append(text^)
        if repetition == 0:
            with open(String(args[3]), "w") as output:
                for row in signatures:
                    output.write(String(row, "\n"))
        var factors = String("[")
        for i in range(len(result.projection.levels)):
            factors.write("," if i else "", len(result.projection.levels[i]))
        factors.write("]")
        print(String('{"repetition":', repetition, ',"construct_ns":', build_ns,
            ',"verify_scalar_ns":', verify_scalar_ns,
            ',"verify_simd_ns":', verify_simd_ns, ',"cells":', cells,
            ',"stacks":', len(result.stacks), ',"factors":', factors,
            ',"materializations":', result.materializations,
            ',"field_map_hits":', result.field_map_hits,
            ',"verification_cache_hits":', check_stats[0],
            ',"verification_cache_misses":', check_stats[1],
            ',"verification_linear_roots":', check_stats[2],
            ',"verification_batches":', check_stats[3],
            ',"verification_word_points":', check_stats[4],
            ',"verification_big_points":', check_stats[5],
            ',"verification_vector_batches":', check_stats[6],
            ',"native_nodes":', native_nodes, ',"wide_nodes":', wide_nodes,
            ',"exact_fallbacks":', exact_fallbacks,
            ',"rational_short_stacks":', rational_short_stacks,
            ',"rational_vector_stacks":', rational_vector_stacks, ',"check_ok":true}'))
