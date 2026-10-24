"""Final projection coverage without determinant or coprimality transcripts."""
from polynomials import Budget, QQ, SparsePoly, MonomialOrder, cad_projection, PolyErrorKind
from polynomials._cad_projection_final import _FinalProjection, _project_values, _check_projection_values
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var work = Budget()
    var c = Checker("final projection values")
    var order = MonomialOrder()
    var variables: List[Int64] = [0, 1, 2]
    for expr in ["x2^2 + x1*x2 + x0^2 - 1", "x0*x2^2 + x2 + x1^2 - 2", "x0*x2^2 + x1*x2 + x0*x1^3 - x0*x1", "x0*x2 - x1*x2 + x1^2 - 2", "x0*x1*x2^2 - x0*x2^3 - x0*x2^2 + x0*x2 - x1*x2^4 - x1*x2 - x1 + x2^4 + x2"]:
        var inputs: List[SparsePoly[QQ]] = [parse_polynomial(expr, order, work)]
        var historical = cad_projection(inputs, variables, work, lazard=True)
        var imported = _FinalProjection(variables.copy(), False, historical.levels.copy())
        c.check(_check_projection_values(inputs, imported, work), "old closure verifies with its histories discarded")
        for brown in range(2):
            var produced = _project_values(inputs, variables, work, brown=Bool(brown))
            c.check(_check_projection_values(inputs, produced, work), "raw projection satisfies complete final obligations")
            var bad = produced.copy()
            _ = bad.levels[2].pop()
            c.check(not _check_projection_values(inputs, bad, work), "missing input coverage rejected")
            bad = produced.copy()
            bad.levels[2].append(bad.levels[2][0].copy())
            c.check(not _check_projection_values(inputs, bad, work), "duplicate noncoprime basis rejected")
            bad = produced.copy()
            bad.levels[2][0] = bad.levels[2][0].mul(bad.levels[2][0], work)
            c.check(not _check_projection_values(inputs, bad, work), "repeated projection factor rejected")
    var inputs: List[SparsePoly[QQ]] = [parse_polynomial("x1^2 + x0*x1 + x0^2 - 1", order, work)]
    var reduced = _project_values(inputs, [Int64(0), Int64(1)], work)
    c.check(_check_projection_values(inputs, reduced, work), "checked coefficient ideal justifies omitted tc")
    reduced.brown = False
    c.check(not _check_projection_values(inputs, reduced, work), "omission cannot claim plain Lazard policy")
    reduced.brown = True
    _ = reduced.levels[0].pop()
    c.check(not _check_projection_values(inputs, reduced, work), "missing mandatory discriminant coverage rejected")
    var limited = Budget(0)
    var exhausted = False
    try:
        _ = _check_projection_values(inputs, reduced, limited)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "final projection check propagates exhaustion")
    c.finish()
