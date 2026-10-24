"""Complete private multivariate CAD, query binding, and final-data mutations."""
from polynomials import Budget, QQ, SparsePoly, MonomialOrder, cad, check_cad, PolyErrorKind
from polynomials._cad_tree_final import _FinalCAD, _verified_final_cad, _check_final_cad, _terminal_signs
from support.corpus import parse_polynomial
from support.check import Checker


def compare(expressions: List[String], dimensions: Int, mut c: Checker, mut budget: Budget) raises:
    var inputs = List[SparsePoly[QQ]]()
    var variables = List[Int64]()
    for expr in expressions:
        inputs.append(parse_polynomial(expr, MonomialOrder(), budget))
    for i in range(dimensions):
        variables.append(Int64(i))
    var old = cad(inputs, variables, budget, lazard=True)
    c.check(check_cad(inputs, old, budget), "historical full CAD verifies")
    var result = _verified_final_cad(inputs, variables, List[Int](), budget)
    c.check(_check_final_cad(inputs, variables, List[Int](), result, budget), "relative full CAD verifies from final data")
    c.check(_check_final_cad(inputs, variables, List[Int](), result, budget, vectorized=False), "scalar and SIMD final verification agree")
    var size = 1
    for _ in inputs:
        size *= 3
    var expected = List[Bool](length=size, fill=False)
    var actual = List[Bool](length=size, fill=False)
    for cell in old.cells:
        if cell.level == dimensions - 1:
            var index = 0
            for sign in cell.signs:
                index = 3 * index + sign.relation + 1
            expected[index] = True
    for stack in result.stacks:
        if stack.level == dimensions - 1:
            for position in range(len(stack.stack.signs)):
                var index = 0
                for sign in _terminal_signs(stack, len(inputs), position, budget):
                    index = 3 * index + sign + 1
                actual[index] = True
    c.check(actual == expected, "complete realized original-input sign sets agree")
    var interval = _verified_final_cad(inputs, variables, List[Int](), budget, interval=True)
    c.check(_check_final_cad(inputs, variables, List[Int](), interval, budget), "interval-produced tree passes the independent exact verifier")
    var same = len(interval.stacks) == len(result.stacks)
    if same:
        for i in range(len(result.stacks)):
            same = same and interval.stacks[i].stack.sectors == result.stacks[i].stack.sectors and interval.stacks[i].stack.signs == result.stacks[i].stack.signs
    c.check(same, "interval and exact producers preserve canonical sectors and all labels")


def main() raises:
    var c = Checker("final multivariate CAD")
    var work = Budget()
    compare(["x0^2 - 2"], 1, c, work)
    compare(["x1^2 - x0", "x0^2 - 2"], 2, c, work)
    compare(["x0*x1 - x0", "x1^2 - 2"], 2, c, work)
    compare(["x0^2 - 2", "x1^2 - 3", "x2 - x0 - x1"], 3, c, work)
    compare(["x0*x2^2 + x1*x2 + x0*x1^3 - x0*x1"], 3, c, work)
    compare(["x0*x3^2 + x1*x3 + x0 + x1"], 4, c, work)
    # Deterministic shifted families force coefficient varieties of dimensions
    # zero and one, rather than relying on random inputs to hit nullification.
    for shift in range(3):
        # Corpus input is an expanded sum of monomials, not a general parser.
        compare([String("x0*x1^2 - ", shift, "*x1^2 + x0*x1 - ", shift, "*x1 + x0 - ", shift)], 2, c, work)
        for dimension in range(3, 5):
            var z = String("x", dimension - 1)
            compare([String("x0*", z, "^2 - ", shift, "*", z, "^2 + x1*", z, " + ", shift + 1, "*", z,
                " + x0*x1 + ", shift + 1, "*x0 - ", shift, "*x1 - ", shift * (shift + 1))], dimension, c, work)
    var inputs: List[SparsePoly[QQ]] = [parse_polynomial("x0^2 - 2", MonomialOrder(), work), parse_polynomial("x1^2 - x0", MonomialOrder(), work)]
    var variables: List[Int64] = [0, 1]
    var masks: List[Int] = [2, 2]
    var result = _verified_final_cad(inputs, variables, masks, work)
    c.check(result.satisfiable, "query has real solutions")
    c.check(result.materializations == 2, "only the two nonterminal irrational base sections materialize")
    var first_embedding = -1
    var second_embedding = -1
    for i in range(len(result.stacks)):
        if result.stacks[i].field.embedding:
            if first_embedding < 0:
                first_embedding = i
            else:
                second_embedding = i
    c.check(first_embedding >= 0 and second_embedding >= 0, "two real embeddings exercise a shared checked modulus")
    for damage in range(15):
        var bad = result.copy()
        if damage == 0:
            bad.satisfiable = False
        elif damage == 1:
            bad.masks[0] = 1
        elif damage == 2:
            bad.inputs[0] = parse_polynomial("x0^2 - 3", MonomialOrder(), work)
        elif damage == 3:
            _ = bad.stacks.pop()
        elif damage == 4:
            bad.stacks[0].children[0] = 0
        elif damage == 5:
            bad.stacks[0].children[1] = bad.stacks[0].children[0]
        elif damage == 6:
            bad.stacks[1].parent = 1
        elif damage == 7:
            bad.stacks[1].slot += 1
        elif damage == 8:
            bad.stacks[1].level = 0
        elif damage == 9:
            bad.stacks[1].truth[0] = not bad.stacks[1].truth[0]
        elif damage == 10:
            bad.stacks[1].orders[0][0] += 1
        elif damage == 11:
            _ = bad.projection.levels[0].pop()
        elif damage == 12:
            bad.stacks[1].field.coordinates[0] = bad.stacks[1].field.coordinates[0].neg(work)
        elif damage == 13:
            # The positive root shares a checked modulus with the negative
            # root. A cached irreducibility fact must not validate its maps.
            bad.stacks[second_embedding].field.coordinates[0] = bad.stacks[second_embedding].field.coordinates[0].neg(work)
        else:
            bad.stacks[second_embedding].field.embedding = bad.stacks[first_embedding].field.embedding.copy()
        c.check(not _check_final_cad(inputs, variables, masks, bad, work), "forged tree, query, specialization or projection rejected")
        c.check(not _check_final_cad(inputs, variables, masks, bad, work, vectorized=False), "scalar verification also rejects the forged result")
    masks = [2, 1]
    result = _verified_final_cad(inputs, variables, masks, work)
    c.check(result.satisfiable and _check_final_cad(inputs, variables, masks, result, work), "inequality query binds its own truth labels")
    inputs = [parse_polynomial("x1^2 + x0^2 + 1", MonomialOrder(), work)]
    masks = [2]
    result = _verified_final_cad(inputs, variables, masks, work)
    c.check(not result.satisfiable, "unsatisfiability requires the complete checked tree")
    var limited = Budget(0)
    var exhausted = False
    try:
        _ = _check_final_cad(inputs, variables, masks, result, limited)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "global exhaustion cannot accept a full result")
    inputs = List[SparsePoly[QQ]]()
    for expr in ["x0^2 - 2", "x1^2 - x0", "x2^2 - x1", "x3^2 - x2", "x4 - x3"]:
        inputs.append(parse_polynomial(expr, MonomialOrder(), work))
    variables = [0, 1, 2, 3, 4]
    masks = List[Int]()
    result = _verified_final_cad(inputs, variables, masks, work, interval=True)
    c.check(result.field_map_hits > 0, "degree-16 nonterminal conjugates exercise the default producer cache")
    c.check(_check_final_cad(inputs, variables, masks, result, work, vectorized=False), "the complete scalar verifier accepts the cached-field tree independently")
    c.finish()
