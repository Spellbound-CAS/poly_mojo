"""Realized input sign vectors must survive revision-15 projection reductions."""
from polynomials import Budget, SparsePoly, QQ, RealAlgebraic, MonomialOrder, cad, check_cad
from polynomials.cad import _cad
from polynomials.cad_types import CADCertificate
from support.corpus import parse_polynomial
from support.check import Checker


def signatures(result: CADCertificate, count: Int, dimensions: Int) -> List[Bool]:
    var size = 1
    for _ in range(count):
        size *= 3
    var found = List[Bool](length=size, fill=False)
    for cell in result.cells:
        if cell.level != dimensions - 1:
            continue
        var index = 0
        for sign in cell.signs:
            index = 3 * index + sign.relation + 1
        found[index] = True
    return found^


def compare(expressions: List[String], dimensions: Int, mut c: Checker, mut work: Budget) raises:
    var inputs = List[SparsePoly[QQ]]()
    var variables = List[Int64]()
    for i in range(dimensions):
        variables.append(Int64(i))
    for expr in expressions:
        inputs.append(parse_polynomial(expr, MonomialOrder(), work))
    var legacy = cad(inputs, variables, work, lazard=True)
    var reduced = _cad(inputs, variables, List[RealAlgebraic](), work, True, brown=True)
    c.check(check_cad(inputs, legacy, work), "plain Lazard result checks")
    c.check(check_cad(inputs, reduced, work), "reduced Lazard result checks")
    c.check(signatures(legacy, len(inputs), dimensions) == signatures(reduced, len(inputs), dimensions), "all realized original-input sign vectors agree")
    var finite = _cad(inputs, variables, List[RealAlgebraic](), work, True, brown=True, brown_finite=True)
    c.check(check_cad(inputs, finite, work), "finite-nullification result checks")
    c.check(signatures(legacy, len(inputs), dimensions) == signatures(finite, len(inputs), dimensions), "finite-nullification input signs agree")


def main() raises:
    var c = Checker("Brown CAD differential")
    var work = Budget()
    compare(["x1^2 + x0*x1 + x0^2 - 1", "x1 - x0"], 2, c, work)
    compare(["x0*x1^2 + x1 + x0^2 - 2", "x1^2 - 1"], 2, c, work)
    compare(["x0*x1^2 + x0*x1 - x1 + x0^2 - 2", "x1 - 1"], 2, c, work)
    compare(["x0*x2^2 + x1*x2 + x0 + x1", "x2 - 1"], 3, c, work)
    compare(["x0*x3^2 + x1*x3 + x0 + x1"], 4, c, work)
    compare(["x0*x2^2 + x1*x2 + x0*x1^3 - x0*x1"], 3, c, work)
    compare(["x0*x2 - x1*x2 + x1^2 - 2"], 3, c, work)
    # Deterministic generated monic, shared-factor and repeated-factor inputs.
    for a in range(-2, 3):
        compare([String("x1^2 + ", a, "*x0*x1 + x0^2 - 2"), "x0*x1 - x0"], 2, c, work)
    c.finish()
