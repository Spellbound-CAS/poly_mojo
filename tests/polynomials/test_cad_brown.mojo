"""Revision-15 omissions need exact coefficient-ideal identities and budgets."""
from polynomials import Budget, PolyErrorKind, QQ, SparsePoly, MonomialOrder, cad_projection, check_cad_projection, groebner, check_groebner_basis
from polynomials.cad_projection import _cad_projection_brown, _prefer_eliminants
from polynomials._cad_brown import _omit_trailing, _coefficient_basis, _check_coefficient_basis, _finite_eliminants
from polynomials.certificates import _equal
from support.corpus import parse_polynomial
from support.check import Checker


def main() raises:
    var c = Checker("Brown Lazard projection")
    var work = Budget()
    var order = MonomialOrder()
    var variables: List[Int64] = [0, 1]
    var lower: List[Int64] = [0]
    var monic = parse_polynomial("x1^2 + x0*x1 + x0^2 - 1", order, work)
    c.check(_omit_trailing(monic, 1, lower, work, cap=0), "constant leading coefficient needs no ideal search")
    var middle = parse_polynomial("x0*x1^2 + x1 + x0^2 - 2", order, work)
    c.check(_omit_trailing(middle, 1, lower, work, cap=0), "constant interior coefficient proves non-nullification")
    var empty = parse_polynomial("x0*x1^2 + x0*x1 - x1 + x0^2 - 2", order, work)
    c.check(_omit_trailing(empty, 1, lower, work), "unit coefficient ideal without a constant coefficient")
    c.check(not _omit_trailing(empty, 1, lower, work, cap=1), "local analysis cap conservatively keeps the trailing coefficient")
    var sparse_degree = parse_polynomial("x0*x1^4294967295 + x0*x1 - x1 + x0^2 - 2", order, work)
    c.check(_omit_trailing(sparse_degree, 1, lower, work), "coefficient analysis depends on sparse terms rather than degree")
    var finite = parse_polynomial("x0*x2^2 + x1*x2 + x0 + x1", order, work)
    c.check(not _omit_trailing(finite, 2, [Int64(0), Int64(1)], work), "finite nonempty nullification is not the empty case")
    var positive_dim = parse_polynomial("x0*x3^2 + x1*x3 + x0 + x1", order, work)
    c.check(not _omit_trailing(positive_dim, 3, [Int64(0), Int64(1), Int64(2)], work), "positive-dimensional nullification retains the factor")
    c.check(not _omit_trailing(parse_polynomial("0", order, work), 1, lower, work), "zero polynomial is never certified non-nullifying")

    for f in [monic.copy(), middle.copy(), empty.copy()]:
        var inputs: List[SparsePoly[QQ]] = [f.copy()]
        var legacy = cad_projection(inputs, variables, work, lazard=True)
        var reduced = _cad_projection_brown(inputs, variables, work)
        c.check(len(reduced.levels[0]) < len(legacy.levels[0]), "reduction removes projection factors")
        c.check(check_cad_projection(inputs, reduced, work), "reduced certificate checks independently")
        c.check(check_cad_projection(inputs, legacy, work), "historical Lazard certificate still checks")
        var bad = reduced.copy()
        _ = bad.factorizations.pop()
        c.check(not check_cad_projection(inputs, bad, work), "missing mandatory factorization rejected")
        bad = reduced.copy()
        _ = bad.levels[0].pop()
        c.check(not check_cad_projection(inputs, bad, work), "missing mandatory projection factor rejected")
        bad = reduced.copy()
        bad.resultants[0].resultant = parse_polynomial("0", order, work)
        c.check(not check_cad_projection(inputs, bad, work), "forged discriminant rejected")

    # The one-sided Groebner predicate accepts {1} for <x0>. This helper must
    # reject it: a zero remainder for the input does not prove 1 belongs to I.
    var lex = MonomialOrder(MonomialOrder.LEX, lower.copy(), work)
    var generators: List[SparsePoly[QQ]] = [parse_polynomial("x0", lex, work)]
    var basis = groebner(generators, lex, work, track_cofactors=True)
    c.check(_check_coefficient_basis(generators, basis, work), "valid coefficient basis identities")
    basis.basis[0] = parse_polynomial("1", lex, work)
    c.check(check_groebner_basis(generators, basis, work), "one-sided predicate alone is insufficient")
    c.check(not _check_coefficient_basis(generators, basis, work), "forged unit basis fails the reverse containment")
    basis.cofactors = None
    c.check(not _check_coefficient_basis(generators, basis, work), "missing generator identities rejected")

    var limited = Budget(0)
    var exhausted = False
    try:
        _ = _coefficient_basis(generators, lower, limited, 1000)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "global budget exhaustion propagates")
    limited = Budget(1000)
    var declined = _coefficient_basis(generators, lower, limited, 10)
    c.check(not declined and limited.remaining() < 1000, "failed bounded analysis charges actual work")

    var axes: List[Int64] = [0, 1]
    var replacement = _finite_eliminants(finite, 2, axes, work)
    c.check(Bool(replacement) and len(replacement.value()) == 2, "finite variety has one eliminant per ambient axis")
    if replacement:
        c.check(_equal(replacement.value()[0], parse_polynomial("x0", order, work), work), "first coordinate of finite nullification")
        c.check(_equal(replacement.value()[1], parse_polynomial("x1", order, work), work), "second coordinate of finite nullification")
    replacement = _finite_eliminants(positive_dim, 3, [Int64(0), Int64(1), Int64(2)], work)
    c.check(not replacement, "unused ambient axis prevents finite reduction")
    var diagonal = parse_polynomial("x0*x2 - x1*x2 + x1^2 - 2", order, work)
    replacement = _finite_eliminants(diagonal, 2, axes, work)
    c.check(Bool(replacement), "lex rotation finds both irrational coordinates")
    if replacement:
        c.check(_equal(replacement.value()[0], parse_polynomial("x0^2 - 2", order, work), work), "rotated lex coordinate eliminant")
        c.check(_equal(replacement.value()[1], parse_polynomial("x1^2 - 2", order, work), work), "original lex coordinate eliminant")
    replacement = _finite_eliminants(empty, 1, lower, work)
    c.check(Bool(replacement) and len(replacement.value()) == 0, "empty variety distinguished from unknown")
    replacement = _finite_eliminants(diagonal, 2, axes, work, cap=1)
    c.check(not replacement, "finite analysis cap keeps tc")

    var axes3: List[Int64] = [0, 1, 2]
    var useful = parse_polynomial("x0*x2^2 + x1*x2 + x0*x1^3 - x0*x1", order, work)
    var inputs: List[SparsePoly[QQ]] = [useful.copy()]
    var first = _cad_projection_brown(inputs, axes3, work)
    var complete = _cad_projection_brown(inputs, axes3, work, finite=True)
    var first_count = 0
    var complete_count = 0
    for level in first.levels:
        first_count += len(level)
    for level in complete.levels:
        complete_count += len(level)
    c.check(complete_count < first_count, "finite replacement wins on full normalized closure")
    c.check(check_cad_projection(inputs, complete, work), "finite replacement certificate checks")
    c.check(check_cad_projection(inputs, first, work), "empty-only historical mode still checks")
    var damaged = complete.copy()
    _ = damaged.factorizations.pop()
    c.check(not check_cad_projection(inputs, damaged, work), "missing finite-mode obligation rejected")
    var known = List[List[SparsePoly[QQ]]](length=2, fill=List[SparsePoly[QQ]]())
    var tc = parse_polynomial("x1^2 - 2", order, work)
    var same: List[SparsePoly[QQ]] = [tc.copy()]
    c.check(not _prefer_eliminants(tc, same, known, axes, work), "equal closure costs retain tc")
    limited = Budget(10000)
    c.check(not _prefer_eliminants(tc, same, known, axes, limited, cap=1), "cost-analysis cap retains tc")
    c.check(limited.remaining() < 10000, "capped cost analysis charges caller")
    c.finish()
