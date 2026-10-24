"""Time bivariate resultant-cofactor construction, including reconstruction.

Arguments: F G ANSWER MIN_SECONDS MIN_REPETITIONS. F and G are canonical
sparse QQ records. Parsing, repeated-answer checks, serialization, and an
extra cofactor identity check are untimed. The producer's exact checks and
internal allocation remain timed.
"""
from polynomials import Budget, SparsePoly, QQ, NoContext
from polynomials.interchange import deserialize_sparse, serialize_sparse
from polynomials._views import _variables
from polynomials._resultant_cofactors import _bivariate_resultant
from std.sys import argv
from std.time import perf_counter_ns


def main() raises:
    var args = argv()
    if len(args) != 6:
        raise Error("expected F G ANSWER MIN_SECONDS MIN_REPETITIONS")
    var limit = atof(args[4]) * 1e9
    var minimum = Int(args[5])
    if limit < 0 or minimum < 1:
        raise Error("invalid repeat limits")
    var work = Budget()
    var ft = open(args[1], "r").read()
    var gt = open(args[2], "r").read()
    var times = List[Int]()
    var total = Int(0)
    var answer: String
    var f = deserialize_sparse[QQ](ft, NoContext(), work)
    var g = deserialize_sparse[QQ](gt, NoContext(), work)
    var variables = _variables(f, g, work)
    if len(variables) != 2:
        raise Error("cofactor input must be bivariate")
    var first = True
    var expected = SparsePoly[QQ].zero(f.order, NoContext(), work)
    var eu = expected.copy()
    var ev = expected.copy()
    while len(times) < minimum or Float64(total) < limit:
        var budget = Budget()
        var start = perf_counter_ns()
        var result = _bivariate_resultant(f, g, variables[0], variables, budget)
        var elapsed = perf_counter_ns() - start
        if not result:
            raise Error("cofactor producer failed")
        ref r = result.value()
        if first:
            expected = r.resultant.copy()
            eu = r.cofactors[0].copy()
            ev = r.cofactors[1].copy()
            first = False
        elif r.resultant != expected or r.cofactors[0] != eu or r.cofactors[1] != ev:
            raise Error("repeated cofactors differ")
        times.append(elapsed)
        total += elapsed
    if eu.mul(f, work).add(ev.mul(g, work), work) != expected:
        raise Error("cofactor identity failed")
    answer = serialize_sparse(expected, work) + serialize_sparse(eu, work) + serialize_sparse(ev, work)
    with open(args[3], "w") as out:
        out.write(answer)
    var record = String('{"times_ns":[')
    for i in range(len(times)):
        if i:
            record += ","
        record += String(times[i])
    print(record + "]}")
