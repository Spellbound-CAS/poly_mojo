"""Assignment order, repeated powers and partial simultaneous substitutions."""
from apn_mojo import Rational, ArithmeticContext, ComplexContext, BallContext
from polynomials import *
from support.check import Checker


def numeric[D: NumericDomain](ctx: D.Context, f: SparsePoly[QQ], x: SparsePoly[QQ], y: SparsePoly[QQ], expected: SparsePoly[QQ], mut check: Checker, mut work: Budget) raises:
    var p = numeric_from_sparse[D, QQ](f, ctx, work)
    var a = numeric_from_sparse[D, QQ](x, ctx, work)
    var b = numeric_from_sparse[D, QQ](y, ctx, work)
    var forward: List[Tuple[Int64, NumericSparsePoly[D]]] = [(Int64(-3), b.copy()), (Int64(8), a.copy())]
    var reverse: List[Tuple[Int64, NumericSparsePoly[D]]] = [(Int64(8), a.copy()), (Int64(-3), b.copy())]
    var first = numeric_substitute_many(p, forward, ctx, work)
    var second = numeric_substitute_many(p, reverse, ctx, work)
    var converted = numeric_from_sparse[D, QQ](expected, ctx, work)
    var valid = len(first.terms) == len(second.terms) and len(first.terms) == len(converted.terms)
    for k in range(len(first.terms)):
        valid = valid and first.terms[k].monomial == second.terms[k].monomial and first.terms[k].monomial == converted.terms[k].monomial
        var box = numeric_bounds[D](first.terms[k].coefficient, work)
        valid = valid and box == numeric_bounds[D](second.terms[k].coefficient, work)
        valid = valid and box.contains(numeric_bounds[D](converted.terms[k].coefficient, work), work)
    check.check(valid, "lookup order leaves numerical arithmetic and simultaneous identity unchanged")
    forward.append((Int64(-3), a.copy()))
    try:
        _ = numeric_substitute_many(p, forward, ctx, work)
        check.check(False, "duplicate numerical assignment accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "duplicate numerical assignment rejected")
    var tiny = Budget(80)
    try:
        _ = numeric_substitute_many(p, reverse, ctx, tiny)
        check.check(False, "substitution budget ignored")
    except e:
        check.check(e.kind == PolyErrorKind.LIMIT and len(p.terms) == len(f.terms), "numerical substitution failure retains source")


def main() raises:
    var check = Checker("substitution_plan")
    var work = Budget()
    var order = MonomialOrder()
    var x = SparsePoly[QQ].variable(-3, order, NoContext(), work)
    var y = SparsePoly[QQ].variable(8, order, NoContext(), work)
    var z = SparsePoly[QQ].variable(99, order, NoContext(), work)
    var f = x.pow(3, work).add(x.pow(2, work).mul(y, work), work).add(x.pow(2, work).mul(z, work), work).add(y.pow(3, work), work)
    var expected = y.pow(3, work).add(y.pow(2, work).mul(x, work), work).add(y.pow(2, work).mul(z, work), work).add(x.pow(3, work), work)
    var replacements: List[Tuple[Int64, SparsePoly[QQ]]] = [(Int64(8), x.copy()), (Int64(-3), y.copy())]
    check.check(substitute_many(f, replacements, work) == expected, "partial substitution keeps unassigned variable and repeated powers")
    replacements.append((Int64(8), y.copy()))
    try:
        _ = substitute_many(f, replacements, work)
        check.check(False, "duplicate exact assignment accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "duplicate exact assignment rejected")
    numeric[FloatDomain](ArithmeticContext(), f, x, y, expected, check, work)
    numeric[ComplexDomain](ComplexContext(ArithmeticContext()), f, x, y, expected, check, work)
    numeric[BallDomain](BallContext(64), f, x, y, expected, check, work)
    numeric[ComplexBallDomain](BallContext(64), f, x, y, expected, check, work)
    check.finish()
