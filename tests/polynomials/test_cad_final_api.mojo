"""Public final CAD evidence, exact output and bounded compact archives."""
from polynomials import Budget, QQ, SparsePoly, MonomialOrder, cad, algebraic_compare, cad_final, check_cad_final, CADFinalCertificate, PolyErrorKind
from polynomials._cad_final_codec import _encode_final
from support.corpus import parse_polynomial
from support.check import Checker


def compare_output(expressions: List[String], variables: List[Int64], mut c: Checker, mut budget: Budget) raises:
    var inputs = List[SparsePoly[QQ]]()
    for text in expressions:
        inputs.append(parse_polynomial(text, MonomialOrder(), budget))
    var old = cad(inputs, variables, budget, lazard=True)
    var result = cad_final(inputs, variables, budget, lazard=True)
    var cells = result.cells(budget)
    c.check(len(cells) == len(old.cells), "same-projection output has the historical cell count")
    for i in range(min(len(cells), len(old.cells))):
        ref actual = cells[i]
        ref expected = old.cells[i]
        c.check(actual.parent == expected.parent and actual.level == expected.level and actual.position == expected.position, "canonical breadth-first geometry matches legacy")
        c.check(algebraic_compare(actual.sample, expected.sample, budget).relation == 0, "canonical sample is the same exact real value")
        var signs = List[Int]()
        for sign in expected.signs:
            signs.append(sign.relation)
        c.check(actual.signs == signs, "original-input labels match legacy")


def main() raises:
    var c = Checker("final CAD public API")
    var budget = Budget()
    var inputs: List[SparsePoly[QQ]] = [parse_polynomial("x0^2 + x1^2 - 1", MonomialOrder(), budget)]
    var variables: List[Int64] = [0, 1]
    var result = cad_final(inputs, variables, budget)
    c.check(result.is_compact(), "default public CAD uses final data")
    c.check(not result.is_covering(), "decomposition claim has no query")
    c.check(check_cad_final(inputs, variables, List[UInt8](), result, budget), "complete public result verifies")
    c.check(check_cad_final(inputs, variables, List[UInt8](), result, budget, vectorized=False), "scalar reference verifies public result")
    c.check(result.variables(budget) == variables, "variable order retained")
    var cells = result.cells(budget)
    c.check(len(cells) == result.cell_count(budget), "absolute geometry accounts for every cell")
    var terminals = 0
    for cell in cells:
        terminals += Int(cell.level == 1)
    c.check(terminals == 13, "circle has thirteen terminal cells")
    var bytes = result.compact_bytes(budget)
    var restored = CADFinalCertificate.from_compact_bytes(Span(bytes), inputs, variables, List[UInt8](), False, budget)
    c.check(restored.compact_bytes(budget) == bytes, "compact payload round trip preserves exact data")
    var wrong: List[SparsePoly[QQ]] = [parse_polynomial("x0^2 + x1^2 - 2", MonomialOrder(), budget)]
    c.check(not check_cad_final(wrong, variables, List[UInt8](), restored, budget), "original polynomial substitution rejected")
    var masks: List[UInt8] = [2]
    var covering = cad_final(inputs, variables, budget, masks=masks, query=True)
    c.check(covering.is_covering() and covering.satisfiable(), "full covering retains query truth")
    c.check(Bool(covering.sample(budget)), "only requested witness path materializes")
    c.check(check_cad_final(inputs, variables, masks, covering, budget, query=True), "covering verifies its exact masks")
    c.check(not check_cad_final(inputs, variables, [UInt8(1)], covering, budget, query=True), "query substitution rejected")
    var legacy = cad_final(inputs, variables, budget, lazard=False)
    c.check(not legacy.is_compact(), "explicit McCallum retains checked legacy backend")
    c.check(check_cad_final(inputs, variables, List[UInt8](), legacy, budget), "legacy backend accepted through new API")
    var rejected = False
    bytes.append(0)
    try:
        _ = CADFinalCertificate.from_compact_bytes(Span(bytes), inputs, variables, List[UInt8](), False, budget)
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    c.check(rejected, "trailing compact bytes rejected")
    rejected = False
    var exhausted = Budget(0)
    try:
        _ = cad_final(inputs, variables, exhausted)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    c.check(rejected, "exhaustion never returns accepted evidence")
    for mutation in range(6):
        var raw = result._data[].compact.value().copy()
        if mutation == 0:
            raw.satisfiable = not raw.satisfiable
        elif mutation == 1:
            _ = raw.stacks[0].stack.roots.pop()
        elif mutation == 2:
            raw.projection.levels[1].clear()
        var damaged = _encode_final(raw, budget)
        if mutation == 3:
            damaged[0] = 0
        elif mutation == 4:
            for k in range(4, 12):
                damaged[k] = 255
        elif mutation == 5:
            _ = damaged.pop()
        rejected = False
        try:
            _ = CADFinalCertificate.from_compact_bytes(Span(damaged), inputs, variables, List[UInt8](), False, budget)
        except e:
            rejected = e.kind == PolyErrorKind.ENCODING
        c.check(rejected, "corrupt final claim or bounded payload is rejected")
    compare_output([String("x0^2 - 2"), String("x1^2 - x0")], [Int64(0), 1], c, budget)
    compare_output([String("x0^2 + x1^2 - 2")], [Int64(1), 0], c, budget)
    compare_output([String("x0^2 - 2"), String("x1 - x0"), String("x2^2 - x1")], [Int64(0), 1, 2], c, budget)
    c.finish()
