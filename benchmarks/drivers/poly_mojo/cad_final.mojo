"""Private rational-stack prototype; timings include all exact validation work.

Print one JSON record per repetition. Construction already produces the final
data, so certificate assembly is included in construct_and_assemble_ns.
Scalar verification uses one word lane with the same exact fallback as SIMD.
This is not a multivariate CAD or QEPCAD comparison.
"""
from apn_mojo import Integer, Rational
from polynomials import Budget, DensePoly, QQ, NoContext
from polynomials._cad_final import _construct_stack, _check_final_stack, _FinalStack
from polynomials._cad_final_signs import _SignStats, _signs_at
from std.sys import argv
from std.time import perf_counter_ns


def poly(values: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var coefficients = List[QQ]()
    for value in values:
        coefficients.append(QQ(Rational(value)))
    return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work)


def data_bytes(stack: _FinalStack) -> Int:
    # Stable textual payload size; no claim that this is an archive format.
    var text = String()
    for b in stack.basis:
        text.write("basis ")
        for coefficient in b.coefficients:
            text.write(coefficient.value, " ")
        text.write("\n")
    for j in range(len(stack.leading)):
        text.write("input ", stack.leading[j].value, " ")
        for e in stack.multiplicities[j]:
            text.write(e, " ")
        text.write("\n")
    for r in stack.roots:
        text.write("root ", r.basis, " ", r.low.value, " ", r.high.value, "\n")
    for s in stack.sectors:
        text.write("sector ", s.value, "\n")
    for row in stack.signs:
        text.write("signs ")
        for s in row:
            text.write(s, " ")
        text.write("\n")
    return text.byte_length()


def input_text(inputs: List[DensePoly[QQ]]) -> String:
    var text = String()
    for f in inputs:
        text.write(f.variable, ":")
        for coefficient in f.coefficients:
            text.write(coefficient.value, ",")
        text.write(";")
    return text^


def measure(name: String, inputs: List[DensePoly[QQ]], repeats: Int) raises:
    var description = input_text(inputs)
    for repetition in range(repeats):
        var build_work = Budget()
        var start = perf_counter_ns()
        var candidate = _construct_stack(inputs, build_work)
        var construction = perf_counter_ns() - start
        var scalar_ns = 0
        var vector_ns = 0
        var vector = _SignStats()
        for turn in range(2):
            var vectorized = (repetition + turn) % 2 == 0
            var work = Budget()
            var stats = _SignStats()
            start = perf_counter_ns()
            var valid = _check_final_stack(inputs, candidate, stats, work, vectorized=vectorized)
            var elapsed = perf_counter_ns() - start
            if not valid:
                raise Error("benchmark candidate failed final verification")
            if vectorized:
                vector_ns = elapsed
                vector = stats.copy()
            else:
                scalar_ns = elapsed
        print(String('{"case":"', name, '","repetition":', repetition,
            ',"input":"', description, '"',
            ',"construct_and_assemble_ns":', construction,
            ',"verify_scalar_ns":', scalar_ns, ',"verify_simd_ns":', vector_ns,
            ',"basis":', len(candidate.basis), ',"roots":', len(candidate.roots),
            ',"cells":', len(candidate.signs), ',"final_data_text_bytes":', data_bytes(candidate),
            ',"word_points":', vector.word_points, ',"big_points":', vector.big_points,
            ',"check_ok":true}'))


def kernel(repeats: Int) raises:
    var work = Budget()
    var p = poly([6, -5, -2, 1], work)
    var points = List[QQ]()
    for i in range(4096):
        points.append(QQ(Rational(i % 31 - 15, i % 7 + 1)))
    for repetition in range(repeats):
        var scalar_ns = 0
        var vector_ns = 0
        var scalar = List[Int]()
        var vector = List[Int]()
        for turn in range(2):
            var vectorized = (repetition + turn) % 2 == 0
            var budget = Budget()
            var stats = _SignStats()
            var start = perf_counter_ns()
            var signs = _signs_at(p, points, stats, budget, vectorized=vectorized)
            var elapsed = perf_counter_ns() - start
            if vectorized:
                vector_ns = elapsed
                vector = signs^
            else:
                scalar_ns = elapsed
                scalar = signs^
        if scalar != vector:
            raise Error("endpoint kernels disagree")
        print(String('{"case":"endpoint-kernel","repetition":', repetition,
            ',"scalar_ns":', scalar_ns, ',"simd_ns":', vector_ns,
            ',"points":', len(points), ',"check_ok":true}'))


def main() raises:
    var args = argv()
    var repeats = Int(args[1]) if len(args) > 1 else 9
    if repeats < 1:
        raise Error("positive repetition count required")
    var work = Budget()
    var inputs = List[DensePoly[QQ]]()
    for k in range(2, 10):
        inputs.append(poly([-k, 0, 1], work))
    measure("eight-quadratics", inputs, repeats)
    inputs.clear()
    for k in range(-4, 5):
        var factor = poly([-k, 1], work)
        inputs.append(factor.mul(factor, work))
    measure("nine-repeated-linears", inputs, repeats)
    inputs = [poly([-1, 0, 1], work), poly([-2, 0, 1], work), poly([1, -2, 1], work)]
    inputs.append(inputs[0].mul(inputs[1], work).scale(QQ(Rational(-3, 7)), work))
    measure("shared-and-repeated", inputs, repeats)
    inputs = [DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-(Integer(1) << 90))), QQ(Rational(0)), QQ(Rational(1))], work)]
    measure("wide-coefficients", inputs, repeats)
    kernel(repeats)
